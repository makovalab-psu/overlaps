#!/usr/bin/env python3
import argparse
import logging
import pathlib
import sys
from utillib import simplewrap
from bfx import intervallib
import read_formats
assert sys.version_info.major >= 3, 'Python 3 required'

BINS = 10
NULL_STR = '.'
DESCRIPTION = simplewrap.wrap(
f"""Analyze the output of overlaps.py --details and print summary statistics.
The output gives the number of errors overlap bases, and the error rates for
each sample. The error rates are broken down into {BINS} bins, one for each
portion of the reads' length.
If --tsv is given, the output is tab-delimited. For each sample, three lines are
produced. They give the number of errors, number of overlap bases, and the error
rates, respectively. Each line has {BINS+3} fields:
1. The sample name.
2. What the numbers measure ('overlaps', 'errors', or 'rates').
3. The value for the sample overall.
4-{BINS+3}. The value for each bin of the read length."""
)
EPILOG = "Note: Differences where one base is N are not counted as errors."


def make_argparser():
  parser = argparse.ArgumentParser(add_help=False, description=DESCRIPTION, epilog=EPILOG,
    formatter_class=argparse.RawDescriptionHelpFormatter)
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
  options.add_argument('-b', '--num-bins', type=int, default=BINS,
    help='The number of bins to use instead of the default. Default: %(default)s')
  options.add_argument('-i', '--intervals', type=argparse.FileType('r'),
    help='A file containing a list of regions to use to filter the data. Only errors and overlaps '
      'contained in these intervals will be counted. This script will make sure the denominator '
      'for error rates only counts bases contained in these intervals (AND overlaps). The format '
      'is one interval per line, in two tab-delimited columns: the start and end coordinates of '
      'the interval (1-based).')
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

  intervals = None
  if args.intervals:
    intervals = intervallib.simplify_intervals(intervallib.read_intervals(args.intervals))

  metastats = {}
  for sample, errors_path_str in args.errors:
    errors_path = pathlib.Path(errors_path_str)
    logging.warning(f'Info: Processing sample {sample}..')
    with errors_path.open() as errors_file:
      pairs = read_formats.read_errors_file(errors_file)
      metastats[sample] = compile_stats(pairs, total_bins=args.num_bins, intervals=intervals)

  if args.format == 'tsv':
    print(format_tsv(make_tsv_data(metastats)), file=args.output)
  elif args.format == 'human':
    table = make_table_data(metastats, min_errors=args.min_errors)
    print(format_table(table, spacing=2), file=args.output)


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


def compile_stats(pairs, total_bins=BINS, intervals=None):
  stats = {
    'overlap': 0,
    'errors': 0,
    'overlap_binned': [0] * total_bins,
    'errors_binned': [0] * total_bins,
  }
  for pair in pairs:
    if not pair.mapped:
      continue
    if intervals is None:
      pair_intervals = None
    else:
      # Get the subset of intervals that cover the overlap, so we don't have to search through all
      # intervals for every calculation in this loop iteration.
      pair_intervals = intervallib.find_intervals(pair.overlap.start, pair.overlap.end, intervals)
    kwargs = dict(total_bins=total_bins, intervals=pair_intervals)
    overlap_binned1 = bin_overlap(pair.overlap.len, pair.len1, **kwargs)
    overlap_binned2 = bin_overlap(pair.overlap.len, pair.len2, **kwargs)
    stats['overlap'] += pair.overlap.len
    add_overlap_bins(stats['overlap_binned'], overlap_binned1, overlap_binned2)
    for error in pair.errors:
      if error.base1 == 'N' or error.base2 == 'N':
        continue
      # Skip errors outside the filter intervals (if given).
      if pair_intervals is not None and intervallib.find_interval(error.pos, pair_intervals) is None:
        continue
      stats['errors'] += 1
      bin1 = get_bin(error.read_coord1, pair.len1, total_bins=total_bins)[0]
      bin2 = get_bin(error.read_coord2, pair.len2, total_bins=total_bins)[0]
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
  Returns a 3-tuple: the (0-based) bin, and the left and right boundaries of that bin, in read
  coordinates. The left and right boundaries are the (1-based) coordinates of the first and last
  bases in the bin."""
  assert num > 0, num
  bin_size = denom/total_bins
  # Subtract 1 to prevent `bin == total_bins` when `num == denom`.
  # That is, the last `num` should be in bin `total_bins-1`.
  # Since `bin` is 0-based, if we didn't do this, there'd be `total_bins+1` total bins.
  bin_float = (num-1)/bin_size
  bin = int(bin_float)
  # The `left_boundary` is where the last bin ended, so it's actually just outside the bin.
  # That is, the last `num` that was in `bin-1`.
  left_boundary_float = (bin*bin_size)+1
  left_boundary = int(left_boundary_float)
  right_boundary_float = ((bin+1)*bin_size)+1
  right_boundary = int(right_boundary_float)
  # When the left_boundary_float is a whole number (11.000..), the left_boundary is actually
  # the previous integer.
  if left_boundary == left_boundary_float:
    left_boundary -= 1
  if right_boundary == right_boundary_float:
    right_boundary -= 1
  # Add 1 to `left_boundary` to get the first base of the bin.
  return bin, left_boundary+1, right_boundary


def bin_overlap(overlap, readlen, total_bins=BINS, intervals=None):
  """Add up how many bases of overlap are in each bin.
  `overlap`: The length of the overlap, in bases. Assumed to be at the end of the read(!)
  `readlen`: The read length.
  `intervals`: A list of intervals to filter by. Only count bases which appear in one of these
    intervals. See `find_interval()` for the format and constraints.
  Returns a list `total_bins` long where each element is the number of overlap bases present in each
  bin."""
  bins_bases = []
  # Get the list of bins as intervals.
  # Walk through the bins with an example coordinate. Start at 1, then increment by the computed
  # size of the current bin so that `coord` is always the first base in the bin.
  coord = 1
  bins = []
  while coord < readlen:
    start, end = get_bin(coord, readlen, total_bins=total_bins)[1:]
    bins.append((start, end))
    bin_size = end - start + 1
    coord += bin_size
  # Intersect the list of bins with the overlap interval. Note that `intersect_intervals()` only
  # returns intervals in the intersection, so any bin completely outside the overlap region does
  # not appear in the list. That is, `len(bins) == total_bins` but `len(overlap_bins) <= total_bins`
  overlap_start = readlen - overlap + 1
  overlap_end = readlen
  overlap_bins = intervallib.intersect_intervals((overlap_start, overlap_end), bins)
  # Add up the number of overlap bases in each bin.
  for overlap_bin in overlap_bins:
    if intervals is None:
      start, end = overlap_bin
      bin_bases = end - start + 1
    else:
      # Intersect the bin with the filter intervals so we get a list of intervals which are only the
      # bases in this bin AND in the overlap AND in one of the intervals.
      bin_intersection = intervallib.intersect_intervals(overlap_bin, intervals)
      bin_bases = 0
      for start, end in bin_intersection:
        bin_bases += end - start + 1
    bins_bases.append(bin_bases)
  # Bins completely outside the overlap region aren't included (see comment above), so we need to
  # left pad the list with zeros (bins with no bases in the overlap).
  bins_bases = [0] * (total_bins-len(bins_bases)) + bins_bases
  return bins_bases


def make_tsv_data(metastats):
  rows = []
  for sample, stats in metastats.items():
    overlap_row = [sample, 'overlaps', stats['overlap']] + [int(o) for o in stats['overlap_binned']]
    rows.append(overlap_row)
    errors_row = [sample, 'errors', stats['errors']] + stats['errors_binned']
    rows.append(errors_row)
    rates = []
    for errors,  overlaps in zip(stats['errors_binned'], stats['overlap_binned']):
      rates.append(divide_or_null(100*errors, overlaps))
    rates_row = [sample, 'rates', divide_or_null(100*stats['errors'], stats['overlap'])] + rates
    rows.append(rates_row)
  return rows


def divide_or_null(numerator, denominator, null=None):
  try:
    return numerator/denominator
  except ZeroDivisionError:
    return null


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
