#!/usr/bin/env python3
import argparse
import logging
import pathlib
import sys
assert sys.version_info.major >= 3, 'Python 3 required'

BINS = 10
NULL_STR = '.'
DESCRIPTION = """Analyze the output of overlaps.py --details and print summary statistics."""


def make_argparser():
  parser = argparse.ArgumentParser(add_help=False, description=DESCRIPTION)
  options = parser.add_argument_group('Options')
  options.add_argument(
    '-e', '--errors', action='append', nargs=2, metavar=('name', 'errors.tsv'), required=True,
    help='The errors output by overlaps.py --details. Give a sample name and the path to the '
      'errors file. Use --errors once for each sample.')
  options.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
    help='Write output to this file instead of stdout.')
  options.add_argument('-t', '--tsv', action='store_const', dest='format', const='tsv', default='human',
    help='Output the counts (total and binned) in tab-delimited format.')
  options.add_argument('-m', '--min-errors', type=int, default=1000,
    help='Minimum number of errors that must be present in a bin in order to output an error rate. '
      'Default: %(default)s')
  options.add_argument('-h', '--help', action='help',
    help='Print this argument help text and exit.')
  logs = parser.add_argument_group('Logging')
  logs.add_argument('-l', '--log', type=argparse.FileType('w'), default=sys.stderr,
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  volume = logs.add_mutually_exclusive_group()
  volume.add_argument('-q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
    default=logging.WARNING)
  volume.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  volume.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)
  return parser


def main(argv):

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')

  samples = set([e[0] for e in args.errors])
  if len(samples) < len(args.errors):
    fail('Error: Sample names must all be unique.')

  metastats = {}
  for sample, errors_path_str in args.errors:
    errors_path = pathlib.Path(errors_path_str)
    logging.warning(f'Info: Processing sample {sample}..')
    metastats[sample] = compile_stats(parse_errors(errors_path))

  if args.format == 'tsv':
    print(format_tsv(make_tsv_data(metastats)), file=args.output)
  elif args.format == 'human':
    table = make_table_data(metastats, min_errors=args.min_errors)
    print(format_table(table, spacing=2), file=args.output)


def parse_errors(path):
  pair = None
  with path.open() as file:
    for line in file:
      fields = line.rstrip('\r\n').split('\t')
      if fields[0] == 'pair':
        if pair is not None:
          assert pair['num_errors'] == len(pair['errors']), (pair['num_errors'], len(pair['errors']), pair['name'])
          yield pair
        pair = parse_pair_line(fields)
      elif fields[0] == 'error':
        error = parse_error_line(fields)
        pair['errors'].append(error)
      else:
        raise ValueError(f'Invalid value in column 1: {fields[0]!r}')
  if pair is not None:
    yield pair


def parse_pair_line(fields):
  return {
    'name': fields[1],
    'mapped': parse_value(fields[2]),
    'len1': parse_value(fields[3]),
    'len2': parse_value(fields[4]),
    'overlap': parse_value(fields[5]),
    'num_errors': parse_value(fields[6]),
    'errors': [],
  }


def parse_error_line(fields):
  return {
    'type': fields[1],
    'ref_coord': parse_value(fields[2]),
    'coord1': parse_value(fields[3]),
    'coord2': parse_value(fields[4]),
    'alt1': fields[5],
    'alt2': fields[6],
  }


def parse_value(value_str):
  if value_str.lower() == 'true':
    return True
  elif value_str.lower() == 'false':
    return False
  elif value_str == NULL_STR:
    return None
  try:
    return int(value_str)
  except ValueError as error:
    error.args = (f'Invalid value string {value_str!r}',)
    raise


def compile_stats(pairs, total_bins=BINS):
  stats = {
    'overlap': 0,
    'errors': 0,
    'overlap_binned': [0]*total_bins,
    'errors_binned': [0]*total_bins,
  }
  for pair in pairs:
    if pair['mapped']:
      overlap_binned1 = bin_overlap(pair['overlap'], pair['len1'])
      overlap_binned2 = bin_overlap(pair['overlap'], pair['len1'])
      stats['overlap'] += pair['overlap']
      add_overlap_bins(stats['overlap_binned'], overlap_binned1, overlap_binned2)
      for error in pair['errors']:
        if error['alt1'] == 'N' or error['alt2'] == 'N':
          continue
        stats['errors'] += 1
        bin1, *rest = get_bin(error['coord1'], pair['len1'], total_bins=total_bins)
        bin2, *rest = get_bin(error['coord2'], pair['len2'], total_bins=total_bins)
        bin = max(bin1, bin2)
        stats['errors_binned'][bin] += 1
  return stats


def add_overlap_bins(totals, overlaps1, overlaps2):
  """Add the binned overlap counts to the running total.
  Do this by averaging the bp in each bin for the two read mates."""
  for bin, (overlap1, overlap2) in enumerate(zip(overlaps1, overlaps2)):
    totals[bin] += (overlap1+overlap2)/2


def get_bin(num, denom, total_bins=BINS):
  """Break the read into `total_bins` bins, and figure out which bin this base is in.
  `num`:   The (1-based) coordinate of the base.
  `denom`: The read length.
  Returns a 2-tuple: the (0-based) bin, and the remainder (the number of bases past the end of the last bin).
  So, if base 41 is in bin 5 and base 42 is in bin 6, then the remainder for 42 will be 1, the remainder for
  43 will be 2, etc."""
  assert num > 0, num
  bin_size = denom/total_bins
  # Subtract 1 to prevent `bin == total_bins` when `num == denom`.
  # That is, the last `num` should be in bin `total_bins-1`.
  # Since `bin` is 0-based, if we didn't do this, there'd be `total_bins+1` total bins.
  bin_float = (num-1)/bin_size
  bin = int(bin_float)
  # The `left_boundary` is where the last bin ended.
  # That is, the last `num` that was in `bin-1`.
  left_boundary_float = (bin*bin_size)+1
  left_boundary = int(left_boundary_float)
  right_boundary_float = ((bin+1)*bin_size)+1
  right_boundary = int(right_boundary_float)
  # When the left_boundary_float is a whole number (11.000..), the left_boundary is actually the previous integer.
  if left_boundary == left_boundary_float:
    left_boundary -= 1
  if right_boundary == right_boundary_float:
    right_boundary -= 1
  return bin, left_boundary, right_boundary


def bin_overlap(overlap, readlen, total_bins=BINS):
  """Add up how many bases of overlap are in each bin.
  `overlap`: The length of the overlap, in bases. Assumed to be at the end of the read(!)
  `readlen`: The read length.
  Returns a list `total_bins` long where each element is the number of overlap bases present in each bin."""
  overlap_binned = []
  non_overlap = readlen - overlap
  boundary_bin, *rest = get_bin(max(non_overlap, 1), readlen, total_bins=total_bins)
  coord = 1
  last_bin = None
  in_overlap = False
  while coord < readlen:
    bin, left, right = get_bin(coord, readlen, total_bins=total_bins)
    if last_bin is not None:
      assert bin - last_bin == 1, (bin, last_bin, coord, readlen)
    bin_size = right - left
    if bin == boundary_bin:
      overlap_in_bin = right - non_overlap
      overlap_binned.append(overlap_in_bin)
      in_overlap = True
    elif in_overlap:
      overlap_binned.append(bin_size)
    else:
      overlap_binned.append(0)
    coord += bin_size
    last_bin = bin
  return overlap_binned


def make_tsv_data(metastats):
  rows = []
  for sample, stats in metastats.items():
    overlap_row = [sample, 'overlaps', stats['overlap']] + [int(o) for o in stats['overlap_binned']]
    rows.append(overlap_row)
    errors_row = [sample, 'errors', stats['errors']] + stats['errors_binned']
    rows.append(errors_row)
    rates = []
    for errors,  overlaps in zip(stats['errors_binned'], stats['overlap_binned']):
      if overlaps > 0:
        rates.append(100*errors/overlaps)
      else:
        rates.append(None)
    rates_row = [sample, 'rates', 100*stats['errors']/stats['overlap']] + rates
    rows.append(rates_row)
  return rows


def format_tsv(rows):
  lines = []
  for row in rows:
    lines.append('\t'.join(map(format_value, row)))
  return '\n'.join(lines)


def format_value(value):
  if value is None:
    return NULL_STR
  elif isinstance(value, float):
    return f'{value:0.5f}'
  else:
    return str(value)


def make_table_data(metastats, min_errors=1000):
  rows = []
  total_bins = len(list(metastats.values())[0]['errors_binned'])
  header1 = ['',      'Total',  'Total',   'Error']
  header2 = ['Sample', 'errors', 'overlap', 'rate']
  for i in range(1, total_bins+1):
    header1.append('')
    header2.append(f'Bin{i}')
  rows.append(header1)
  rows.append(header2)
  for sample, stats in metastats.items():
    values = []
    overlaps = stats['overlap']
    errors = stats['errors']
    for bin in range(total_bins):
      overlap_binned = stats['overlap_binned'][bin]
      errors_binned = stats['errors_binned'][bin]
      if overlap_binned <= 0:
        values.append('.')
      elif errors_binned <= min_errors:
        values.append('?')
      else:
        rate = errors_binned/overlap_binned
        values.append(f'{100*rate:0.2f}%')
    err_rate = f'{100*errors/overlaps:0.2f}%'
    row = [sample, errors, get_human_bp(overlaps), err_rate] + values
    rows.append(row)
  return rows


def get_human_bp(bp):
  units = ((0, 'bp'), (1, 'kb'), (2, 'Mb'), (3, 'Gb'), (4, 'Pb'))
  for exp, unit in reversed(units):
    bp_adj = bp/(1000**exp)
    if bp_adj >= 1:
      if unit == 'bp':
        return f'{bp} {unit}'
      else:
        return f'{bp_adj:0.1f} {unit}'


def format_table(rows, spacing=1):
  max_widths = get_max_widths(rows)
  spacer = ' '*spacing
  lines = []
  for row in rows:
    value_strs = []
    for i, value in enumerate(row):
      width = max_widths[i]
      format_str = '{{:{}s}}'.format(width)
      value_str = format_str.format(str(value))
      value_strs.append(value_str)
    lines.append(spacer.join(value_strs))
  return '\n'.join(lines)


def get_max_widths(rows):
  max_widths = []
  for row in rows:
    for i, value in enumerate(row):
      str_len = len(str(value))
      if len(max_widths) <= i:
        max_widths.append(str_len)
      else:
        max_widths[i] = max(max_widths[i], str_len)
  return max_widths


def fail(message):
  logging.critical('Error: '+str(message))
  if __name__ == '__main__':
    sys.exit(1)
  else:
    raise Exception(message)


if __name__ == '__main__':
  try:
    sys.exit(main(sys.argv))
  except BrokenPipeError:
    pass
