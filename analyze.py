#!/usr/bin/env python3
import argparse
import logging
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


class LogAccumulator:
  def __init__(self, logger=logging, accumulate=True, max_lines=25):
    self.logger = logger
    self.max_lines = max_lines
    self.accumulate = accumulate
    self._lines = []
  def log(self, level, message, *args, **kwargs):
    if self.accumulate:
      self._lines.append((level, message, args, kwargs))
      while len(self._lines) > self.max_lines:
        self._lines.pop(0)
    else:
      self.logger.log(level, message, *args, **kwargs)
  def debug(self, message, *args, **kwargs):
    self.log(logging.DEBUG, message, *args, **kwargs)
  def info(self, message, *args, **kwargs):
    self.log(logging.INFO, message, *args, **kwargs)
  def warning(self, message, *args, **kwargs):
    self.log(logging.WARNING, message, *args, **kwargs)
  def error(self, message, *args, **kwargs):
    self.log(logging.ERROR, message, *args, **kwargs)
  def critical(self, message, *args, **kwargs):
    self.log(logging.CRITICAL, message, *args, **kwargs)
  def dump_all(self):
    for level, message, args, kwargs in self._lines:
      self.logger.log(level, message, *args, *kwargs)


logger = LogAccumulator()


def make_argparser():
  parser = argparse.ArgumentParser(add_help=False, description=DESCRIPTION, epilog=EPILOG,
    formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument('errors', type=argparse.FileType('r'),
    help='The errors output by overlaps.py --details.')
  options = parser.add_argument_group('Options')
  options.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
    help='Write output to this file instead of stdout.')
  options.add_argument('-t', '--tsv', action='store_const', dest='format', const='tsv', default='human',
    help='Output the counts (total and binned) in tab-delimited format.')
  options.add_argument('-m', '--min-errors', type=int, default=1000,
    help='Minimum number of errors that must be present in a bin in order to output an error rate. '
      'Default: %(default)s')
  options.add_argument('-b', '--num-bins', type=int, default=BINS,
    help='The number of bins to use instead of the default. Default: %(default)s')
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

  pairs = read_formats.read_errors_file(args.errors)
  stats = compile_stats(pairs, total_bins=args.num_bins)

  if args.format == 'tsv':
    print(format_tsv(make_tsv_data(stats)), file=args.output)
  elif args.format == 'human':
    table = make_table_data(stats, min_errors=args.min_errors)
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


def init_analysis(total_bins=BINS):
  return {
    'overlap': 0,
    'errors': 0,
    'overlap_binned': [0] * total_bins,
    'errors_binned': [0] * total_bins,
  }


def compile_stats(pairs, total_bins=BINS):
  stats = init_analysis(total_bins=total_bins)
  for pair in pairs:
    if not pair.mapped:
      continue
    overlap_binned1 = bin_overlap(pair.overlap.len, pair.len1, total_bins=total_bins)
    overlap_binned2 = bin_overlap(pair.overlap.len, pair.len2, total_bins=total_bins)
    stats['overlap'] += pair.overlap.len
    add_overlap_bins(stats['overlap_binned'], overlap_binned1, overlap_binned2)
    for error in pair.errors:
      if error.base1 == 'N' or error.base2 == 'N':
        continue
      stats['errors'] += 1
      bin1 = get_bin(error.read_coord1, pair.len1, total_bins=total_bins)[0]
      bin2 = get_bin(error.read_coord2, pair.len2, total_bins=total_bins)[0]
      bin = max(bin1, bin2)
      stats['errors_binned'][bin] += 1
  return stats


def analyze_pair(analysis, pair, errors, overlap, total_bins=BINS, intervals=None):
  readlen1, readlen2 = len(pair[0].seq), len(pair[1].seq)
  #TODO: Use an overlap length measure that's the number of bases in the overlap shared by both
  #      reads. It's currently the number of overlap bases in one of the reads, but that can
  #      include positions in that read that don't exist in the other (indels). Instead, the
  #      error rate denominator should be the number of opportunities for finding an error, which
  #      means every time a base in read 1 aligns with a base in read 2.
  logger.debug('Analyzing pair %s:', pair.first and pair.first.qname)
  analysis['overlap'] += overlap.length
  if intervals is None:
    logger.debug('  No interval filters.')
    filtered_intervals = pair_intervals1 = pair_intervals2 = None
  elif overlap.length > 0:
    # Get the subset of intervals that cover the overlap, so we don't have to search through all
    # intervals for every calculation in this loop iteration.
    logger.debug('  Overlap: %sbp.', overlap.length)
    filtered_intervals = intervallib.find_intervals(overlap.start, overlap.end, intervals)
    logger.debug('  Filtered intervals: %s', filtered_intervals)
    pair_intervals1 = get_intervals_read_coords(pair[0], filtered_intervals)
    logger.debug('    -> read 1 coords: %s', pair_intervals1)
    pair_intervals2 = get_intervals_read_coords(pair[1], filtered_intervals)
    logger.debug('    -> read 2 coords: %s', pair_intervals2)
  else:
    # There's no overlap. Set `pair_intervals` to empty set, which is technically correct (there
    # are no intervals that cover the overlap, since it doesn't exist). That won't pose problems
    # for the later calculations, since we will want to count no overlap bases or errors.
    logger.debug('  No overlap.')
    filtered_intervals = pair_intervals1 = pair_intervals2 = ()
  kwargs = dict(total_bins=total_bins)
  overlap_binned1 = bin_overlap(overlap.length, readlen1, intervals=pair_intervals1, **kwargs)
  overlap_binned2 = bin_overlap(overlap.length, readlen2, intervals=pair_intervals2, **kwargs)
  if overlap.length > 0:
    logger.debug('  Overlap bases in read 1 bins: %s', overlap_binned1)
    logger.debug('  Overlap bases in read 2 bins: %s', overlap_binned2)
  add_overlap_bins(analysis['overlap_binned'], overlap_binned1, overlap_binned2)
  for error in errors:
    logger.debug(
      '  Error @ ref %(ref_coord)s/read1 %(read_coord1)s/read2 %(read_coord2)s: %(base1)s -> '
      '%(base2)s', error._asdict()
    )
    if error.base1 == 'N' or error.base2 == 'N':
      logger.debug('    Skipping (contains N).')
      continue
    if (
      filtered_intervals is not None and
      intervallib.find_interval(error.ref_coord, filtered_intervals) is None
    ):
      logger.debug('    Skipping (not in an interval).')
      continue
    analysis['errors'] += 1
    bin1 = get_bin(error.read_coord1, readlen1, total_bins=total_bins)[0]
    bin2 = get_bin(error.read_coord2, readlen2, total_bins=total_bins)[0]
    bin = max(bin1, bin2)
    logger.debug('    Read1 bin: %s, read2 bin: %s -> final bin: %s', bin1+1, bin2+1, bin+1)
    analysis['errors_binned'][bin] += 1


def get_intervals_read_coords(read, ref_intervals):
  """Convert a list of intervals in reference coordinates to a list of intervals in read
  coordinates."""
  read_intervals = []
  readlen = len(read.seq)
  read_start_ref = read.pos
  read_end_ref = read.get_end_position()
  for interval_start_ref, interval_end_ref in ref_intervals:
    # Convert the start coord.
    interval_start_read = read.to_read_coord(interval_start_ref)
    logger.debug('      %s -> %s', interval_start_ref, interval_start_read)
    if interval_start_read is None:
      if interval_start_ref < read_start_ref:
        logger.debug('        %s < %s', interval_start_ref, read_start_ref)
        if read.forward:
          interval_start_read = -1
        else:
          interval_start_read = readlen + 1
      else:
        # We're not before the start, so we must be in an indel. Try walking away from this position
        # and find the first base outside the indel.
        interval_start_read = find_closest_read_base(read, interval_start_ref)
      logger.debug('          start = %s', interval_start_read)
    # Convert the end coord.
    interval_end_read = read.to_read_coord(interval_end_ref)
    logger.debug('      %s -> %s', interval_end_ref, interval_end_read)
    if interval_end_read is None:
      logger.debug('        interval_end_read is None')
      if interval_end_ref > read_end_ref:
        logger.debug('        %s > %s', interval_end_ref, read_end_ref)
        if read.forward:
          interval_end_read = readlen + 1
        else:
          interval_end_read = -1
      else:
        # In an indel (see above).
        interval_end_read = find_closest_read_base(read, interval_end_ref)
        logger.debug('        found interval_end_read: %s', interval_end_read)
      logger.debug('          end = %s', interval_end_read)
    # Make sure they're in the right order. A read in reverse orientation could switch this.
    try:
      if interval_start_read > interval_end_read:
        interval_start_read, interval_end_read = interval_end_read, interval_start_read
    except TypeError:
      logger.dump_all()
      print(f'interval_start_read: {interval_start_read}, interval_end_read: {interval_end_read}')
      raise
    # Add it to the list, if it overlaps the read.
    if interval_start_read < 0 and interval_end_read < 0:
      pass # The interval is to the left of the read and doesn't overlap it.
    elif interval_start_read > readlen and interval_end_read > readlen:
      pass # The interval is to the right of the read and doesn't overlap it.
    else:
      read_intervals.append((interval_start_read, interval_end_read))
  # Make sure they're sorted (by start coordinate).
  read_intervals.sort()
  return read_intervals


def find_closest_read_base(read, ref_coord):
  logger.debug('  -- finding closest read base to %d --', ref_coord)
  read_coord = None
  distance = 0
  while read_coord is None:
    distance += 1
    read_coord = read.to_read_coord(ref_coord-distance)
    if read_coord is not None:
      return read_coord
    read_coord = read.to_read_coord(ref_coord+distance)
  logger.debug('    -- went %d bases away (to %d - %d) --', distance, ref_coord-distance, ref_coord+distance)
  logger.debug('    -- result: %s --', read_coord)
  return read_coord


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
      try:
        bin_intersection = intervallib.intersect_intervals(overlap_bin, intervals)
      except AssertionError:
        print(overlap, readlen, overlap_bin, intervals, sep='\n')
        raise
      if len(bin_intersection) == 0 or bin_intersection[0] != overlap_bin:
        logger.debug('    %ss & intervals -> %s', overlap_bin, bin_intersection)
      bin_bases = 0
      for start, end in bin_intersection:
        bin_bases += end - start + 1
    bins_bases.append(bin_bases)
  # Bins completely outside the overlap region aren't included (see comment above), so we need to
  # left pad the list with zeros (bins with no bases in the overlap).
  bins_bases = [0] * (total_bins-len(bins_bases)) + bins_bases
  return bins_bases


def make_tsv_data(stats):
  rows = []
  overlap_row = ['overlaps', stats['overlap']] + [round(o) for o in stats['overlap_binned']]
  rows.append(overlap_row)
  errors_row = ['errors', stats['errors']] + stats['errors_binned']
  rows.append(errors_row)
  rates = []
  for errors, overlaps in zip(stats['errors_binned'], stats['overlap_binned']):
    rates.append(divide_or_null(100*errors, overlaps))
  rates_row = ['rates', divide_or_null(100*stats['errors'], stats['overlap'])] + rates
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
  elif value == 0:
    return '0'
  elif isinstance(value, float):
    return f'{value:0.5f}'
  else:
    return str(value)


def make_table_data(stats, min_errors=1000):
  rows = []
  total_bins = len(stats['errors_binned'])
  header1 = ['Total',  'Total',   'Error']
  header2 = ['errors', 'overlap', 'rate']
  for i in range(1, total_bins+1):
    header1.append('')
    header2.append(f'Bin{i}')
  rows.append(header1)
  rows.append(header2)
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
  row = [errors, get_human_bp(overlaps), err_rate] + values
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
  logger.critical('Error: '+str(message))
  logger.dump_all()
  if __name__ == '__main__':
    sys.exit(1)
  else:
    raise Exception(message)


if __name__ == '__main__':
  try:
    sys.exit(main(sys.argv))
  except BrokenPipeError:
    pass
