#!/usr/bin/env python3
import argparse
import collections
import logging
import pathlib
import subprocess
import sys
import time
import typing
import analyze
from bfx import cigarlib
from bfx import samreader
from bfx import intervallib
from utillib import simplewrap
assert sys.version_info.major >= 3, 'Python 3 required'

BINS = 10
VALUES_TO_STRS = {None:'.'}
ERROR_FIELDS = (
  'type', 'ref', 'ref_coord', 'read_coord1', 'read_coord2', 'base1', 'base2', 'qual1', 'qual2'
)
DESCRIPTION = simplewrap.wrap(
f"""Use the overlap between paired-end reads to find sequencing errors.
Currently only detects SNVs.
The output is tab-delimited, unless --human is given.
Null values are indicated by '{VALUES_TO_STRS.get(None, None)}', and boolean values by \
'{VALUES_TO_STRS.get(True, True)}' and '{VALUES_TO_STRS.get(False, False)}'.

The "summary" output is three lines:
The first two give summary statistics on the read lengths and overlap lengths, respectively.
Each of these first two lines contain four tab-delimited fields:
minimum, average, median, and maximum.
The third (and last) line contains 8 tab-delimited fields:
1. The total number of errors.
2. The total number of overlap bases.
3. The total number of well-mapped read pairs.
4. The total number of read pairs.
5. The total number of input bases.
6. The number of errors per base of overlap (column 1 / column 2).
7. The proportion of reads that were in well-mapped pairs (column 3 / column 4).
8. The proportion of input bases that were in overlaps (2 * column 2 / column 5).

The "details" output contains two types of lines:
Each line's type is indicated in the first column.
'pair' lines encode information on each pair of reads. The columns are:
1. 'pair'
2. The name of the first read in the pair.
3. Whether it contains two well-mapped reads (boolean).
4. Length of the first read. Null if no read 1.
5. Length of the second read. Null if no read 2.
6. Length of the overlap between the two reads. This is the number of read bases in the overlap,
   but only one of the reads, not both.
7. How many errors were detected in the overlap.
8. Whether the pair is oriented with mate 1 forward and mate 2 reverse.
9. The first base of the overlap in reference coordinates.
10. The last base of the overlap in reference coordinates.
'error' lines encode information on each error detected in the overlap of
the pair described in the preceding 'pair' line. The columns are:
1. 'error'
2. The type of error: 'snv', 'ins', or 'del'.
3. The name of the reference sequence the error is in.
4. The reference coordinate of the error.
5. The coordinate of the error in read 1.
6. The coordinate of the error in read 2.
7. The allele present in read 1.
8. The allele present in read 2.
9. The quality score for the error base in read 1.
10 The quality score for the error base in read 2."""
)


def make_argparser():
  parser = argparse.ArgumentParser(
    description=DESCRIPTION, add_help=False, formatter_class=argparse.RawDescriptionHelpFormatter
  )
  io = parser.add_argument_group('I/O')
  io.add_argument('align', type=pathlib.Path, nargs='?', default=sys.stdin,
    help='Name-sorted BAM or SAM file. Omit to read from stdin.')
  io.add_argument('-H', '--human', dest='format', action='store_const', const='human', default='tsv',
    help='Print human-readable text instead of tab-delimited output.')
  io.add_argument('-o', '--outfile', type=argparse.FileType('w'), default=sys.stdout,
    help='Print output to this file instead of stdout. Warning: If the file already exists, it '
      'will be overwritten.')
  io.add_argument('-s', '--summary', type=argparse.FileType('w'),
    help='Print summary stats for the entire input to this file.')
  io.add_argument('-a', '--analysis', type=argparse.FileType('w'),
    help='Print an analysis to this file. This is the same summary analyze.py prints. Currently '
      'always prints in tsv format.')
  filters = parser.add_argument_group('Filters')
  filters.add_argument('-m', '--mapq', type=int,
    help='Mapping quality threshold. Alignments worse than this MAPQ will be ignored.')
  #TODO:
  # filters.add_argument('-q', '--phred', type=int,
  #   help="Don't count bases with PHRED scores below this as errors.")
  # filters.add_argument('-N', '--ignore-ns', action='store_true',
  #   help="Don't count N's as errors.")
  filters.add_argument('-i', '--intervals', type=argparse.FileType('r'),
    help='A file containing a list of regions to use to filter the errors. Only errors and '
      'overlaps contained in these intervals will be counted. This script will make sure the '
      'denominator for error rates only counts bases contained in these intervals (AND overlaps). '
      'The format is one interval per line, in two tab-delimited columns: the start and end '
      'coordinates of the interval (1-based). NOTE: The stats in the "summary" output which '
      'involve the number of overlap bases will use the total overlap, *not* filtered by these '
      'intervals.')
  options = parser.add_argument_group('Options')
  options.add_argument('-p', '--progress', type=int, default=0,
    help='Print periodic updates to stderr while processing large files. Give a number X and this '
      'will print human-readable stats every X minutes (as well as an initial update X/10 minutes '
      'in). Off by default.')
  options.add_argument('-S', '--check-order', action='store_true',
    help="Validate that the input is name-sorted, according to Python's version of lexicographic "
      'ordering. This will cause the script to fail if a read is discovered out of order.')
  options.add_argument('-P', '--print', action='store_true',
    help='Debug print the alignment as cigarlib sees it.')
  options.add_argument('-R', '--detailed-read',
    help='When using --print, single out this read and print every reference coordinate of every '
      'base. Give the full read name, ending in a /1 or /2 for the mate.')
  options.add_argument('-h', '--help', action='help',
    help='Print this argument help text and exit.')
  logs = parser.add_argument_group('Logging')
  logs.add_argument('-l', '--log', type=argparse.FileType('w'), default=sys.stderr,
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  volume = logs.add_mutually_exclusive_group()
  volume.add_argument('-Q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
    default=logging.WARNING)
  volume.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  volume.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)
  return parser


def main(argv):

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')

  # Process arguments.
  intervals = None
  if args.intervals:
    raw_intervals = intervallib.read_intervals_bed(args.intervals)
    intervals = intervallib.simplify_intervals_chrom(raw_intervals)
  align_file = open_input(args.align)

  if args.print:
    print_alignment(align_file, args.detailed_read)
    return

  overlapper = Overlapper(align_file, check_order=args.check_order)
  process_file(overlapper, args.mapq, intervals, args.format, args.progress, args.outfile)

  if args.summary:
    if args.progress and args.format == 'human':
      print('Finished!', file=sys.stderr)
    print(
      format_summary_stats(overlapper.stats, overlapper.counters, args.format), file=args.summary
    )
  if args.analysis:
    tsv_data = analyze.make_tsv_data(overlapper.analysis)
    print(analyze.format_tsv(tsv_data), file=args.analysis)


def process_file(overlapper, mapq_thres, intervals, format, progress, errors_file):
  progress_sec = progress * 60
  start = int(time.time())
  last = None
  for errors, pair, overlap in overlapper.analyze_overlaps(mapq_thres, intervals):
    print(format_read_stats(errors, pair, overlap, format=format), end='', file=errors_file)
    is_progress_time, last = is_it_progress_time(progress_sec, start, last)
    if is_progress_time:
      now = int(time.time())
      output = f'\tProcessed {overlapper.stats["reads"]} reads after {human_time(now-start)}:\n'
      output += format_summary_stats(overlapper.stats, overlapper.counters, format='human')
      print(output, file=sys.stderr)
      sys.stderr.flush()


def open_input(align_path):
  if align_path is sys.stdin:
    return align_path
  elif not align_path.is_file():
    fail(f'Input file must be a regular file. {str(align_path)!r} is not.')
  if align_path is None or str(align_path) == '-':
    return sys.stdin
  else:
    if align_path.suffix == '.bam':
      return open_bam(align_path)
    else:
      return align_path.open()


def open_bam(bam_path):
  cmd = ('samtools', 'view', bam_path)
  process = subprocess.Popen(cmd, stdout=subprocess.PIPE, encoding='utf8')
  for line in process.stdout:
    yield line


def filter_errors(raw_errors, intervals):
  if intervals is None:
    return list(raw_errors)
  filtered_errors = []
  for error in raw_errors:
    if intervallib.find_interval(error.ref_coord, intervals) is not None:
      filtered_errors.append(error)
  return filtered_errors


def get_chrom_intervals(intervals, pair):
  if intervals is None:
    return None
  chrom = pair.get('rname', throw=False)
  if chrom is None:
    return None
  try:
    return intervals[chrom]
  except KeyError:
    return None


class Overlapper:

  def __init__(self, align_file, check_order=False, total_bins=BINS):
    self.align_file = align_file
    self.check_order = check_order
    self.stats = {'reads':0, 'pairs':0, 'pair_bases':0, 'overlap_bp':0, 'errors':0}
    self.counters = {'overlap_lens':collections.Counter(), 'read_lens':collections.Counter()}
    self.analysis = analyze.init_analysis(total_bins=total_bins)

  def analyze_overlaps(self, mapq_thres=None, intervals=None):
    for pair in self.get_pairs(incomplete=True):
      self.stats['reads'] += pair.num_reads
      for read in pair:
        if read is not None:
          self.counters['read_lens'][len(read.seq)] += 1
      if not (pair.is_full and pair.is_well_mapped(mapq_thres)):
        yield [], pair, None
        continue
      chrom_intervals = get_chrom_intervals(intervals, pair)
      if intervals is not None and chrom_intervals is None:
        # We're filtering by intervals, and this chromosome doesn't have any intervals listed.
        continue
      raw_errors, overlap = get_mismatches(pair)
      #TODO: Errors are getting filtered twice here: Once in analyze_pair() and once in this
      #      function. It'd be great to eliminate the duplication of effort, though it might not
      #      be a signficant performance impact. And it's conceptually more intuitive for
      #      analyze_pair() to filter both the overlap bases and the errors by the intervals, rather
      #      than just the overlap bases.
      analyze.analyze_pair(self.analysis, pair, raw_errors, overlap, intervals=chrom_intervals)
      errors = filter_errors(raw_errors, chrom_intervals)
      self.stats['pair_bases'] += len(pair[0].seq) + len(pair[1].seq)
      self.stats['overlap_bp'] += overlap.length
      self.stats['errors'] += len(errors)
      self.stats['pairs'] += 1
      self.counters['overlap_lens'][overlap.length] += 1
      yield errors, pair, overlap

  def get_pairs(self, incomplete=False):
    pair = Pair()
    last_name = None
    for read in samreader.read(self.align_file):
      # logging.debug(mated_name(read))
      name = read.qname
      if self.check_order and last_name is not None and name < last_name:
        fail(f'Reads must be sorted by name! Failed on:\n  {last_name}\n  {name}')
      last_name = name
      assert not (name.endswith('/1') or name.endswith('/2')), name
      assert read.mate in (1, 2), read.flag
      # logging.debug(f'pair: {format_pair(pair)}')
      try:
        pair.add(read)
      except InvalidState:
        fail(
          'Invalid state. Encountered a full pair too early. Pair:\n'
          f'  {pair[0].qname}\n  {pair[1].qname}'
        )
      if not pair.is_full:
        continue
      # Check the full pairs.
      if pair[0].qname != pair[1].qname:
        # If the names don't match, the first read can't have a mate in the file, since this is
        # name-sorted. Discard it, shift the second read to the first position, and loop again.
        singleton = pair.shift()
        if incomplete:
          yield Pair(singleton)
      else:
        # Seems like a properly matched pair. Yield it, and start the next pair.
        yield pair
        pair = Pair()


Overlap = collections.namedtuple('Overlap', ('length', 'start', 'end'))


def get_mismatches(pair):
  errors = []
  overlap_len = 0
  read1, read2 = pair
  base_map = get_base_map(read2)
  started = False
  read_coords, ref_coords = get_coords(read1)
  max_ref_coord = -1
  min_ref_coord = sys.maxsize
  # Iterate through every read coordinate and its corresponding reference coordinate (if any).
  for read1_coord, ref_coord, base1 in zip(read_coords, ref_coords, read1.seq):
    if ref_coord is None:
      # logging.debug(f'None: {base1} -> <INS>')
      continue
    try:
      base2, read2_coord = base_map[ref_coord]
    except KeyError:
      #TODO: Insertions don't have a reference position. How to detect them in read2?
      # if started:
      #   logging.debug(f'{ref_coord:4d}: {base1} -> <DEL>')
      continue
    started = True
    max_ref_coord = max(max_ref_coord, ref_coord)
    min_ref_coord = min(min_ref_coord, ref_coord)
    overlap_len += 1
    if base1 != base2:
      # logging.debug(f'{ref_coord:4d}: {base1} -> {base2}')
      qual1 = read1.qual[read1_coord-1]
      qual2 = read2.qual[read2_coord-1]
      error = Error(
        type='snv',
        rname=read1.qname,
        ref=read1.rname,
        read_coord1=read1_coord,
        read_coord2=read2_coord,
        ref_coord=ref_coord,
        base1=base1,
        base2=base2,
        qual1=qual1,
        qual2=qual2,
      )
      errors.append(error)
  if max_ref_coord < 0:
    max_ref_coord = None
  if min_ref_coord >= sys.maxsize:
    min_ref_coord = None
  return errors, Overlap(overlap_len, min_ref_coord, max_ref_coord)


def get_base_map(read):
  """Return a mapping of reference coordinates to read bases and coordinates.
  This returns a dict where each key is `ref_coord` and each value is a tuple of
  `(read_coord, read_base)`."""
  base_map = {}
  read_coords, ref_coords = get_coords(read)
  for read_coord, ref_coord, read_base in zip(read_coords, ref_coords, read.seq):
    base_map[ref_coord] = (read_base, read_coord)
  return base_map


def get_coords(read):
  read_coords = range(1, len(read.seq)+1)
  if read.reversed:
    read_coords = reversed(read_coords)
  ref_coords = get_reference_positions(read)
  return read_coords, ref_coords


def is_it_progress_time(progress_sec, start, last):
  if progress_sec == 0:
    return False, 0
  now = time.time()
  if last is None:
    # We haven't had a progress time yet.
    # Is it the first progress time?
    elapsed = now - start
    if elapsed > progress_sec/10:
      return True, start
  else:
    # Is it the next progress time?
    elapsed = now - last
    if elapsed > progress_sec:
      return True, now
  return False, last


def format_read_stats(errors, pair, overlap, format='tsv'):
  if format == 'tsv':
    return format_read_stats_tsv(errors, pair, overlap)
  elif format == 'human':
    return format_read_stats_human(errors, pair, overlap)


def format_read_stats_tsv(errors, pair, overlap):
  lines = []
  if pair.is_full:
    if overlap is None:
      olen = ostart = oend = None
    else:
      olen = overlap.length
      ostart = overlap.start
      oend = overlap.end
    read_line = (
      'pair', pair[0].qname, pair.is_well_mapped(cached=True), len(pair[0].seq), len(pair[1].seq),
      olen, len(errors), pair.forward, ostart, oend
    )
  else:
    read_line = ('pair', pair[0].qname, False, len(pair[0].seq), None, None, 0, None, None, None)
  lines.append(read_line)
  for error in errors:
    fields = [getattr(error, field) for field in ERROR_FIELDS]
    lines.append(['error']+fields)
  output = values_to_tsv(lines)
  if output:
    return output+'\n'
  else:
    return ''


def format_read_stats_human(errors, pair, overlap):
  lines = []
  if pair.is_full:
    if pair.is_well_mapped(cached=True):
      first_line = f'{pair[0].qname}: {len(errors)} errors in {overlap.length}bp'
    else:
      first_line = f'{pair[0].qname}: Not well-mapped.'
  else:
    first_line = f'{pair[0].qname}: Singleton.'
  if errors:
    first_line += ':'
  lines.append(first_line)
  for error in errors:
    lines.append(f'{error.ref_coord:4d} {error.type.upper()}: {error.base1} -> {error.base2}')
  return '\n'.join(lines)+'\n'


def format_summary_stats(stats, counters, format='tsv'):
  add_computed_stats(stats, counters)
  if format == 'tsv':
    return format_summary_stats_tsv(stats)
  elif format == 'human':
    return format_summary_stats_human(stats)


def format_summary_stats_tsv(stats):
  output_stats = (
    ('min_rlen', 'avg_rlen', 'med_rlen', 'max_rlen'),
    ('min_overlap', 'avg_overlap', 'med_overlap', 'max_overlap'),
    ('errors', 'overlap_bp', 'pairs', 'reads', 'pair_bases', 'error_rate', 'paired_read_frac', 'overlap_rate'),
  )
  lines = []
  for line_stats in output_stats:
    lines.append([stats[stat] for stat in line_stats])
  return values_to_tsv(lines)


def values_to_tsv(value_lines):
  str_lines = []
  for value_line in value_lines:
    strs = [VALUES_TO_STRS.get(value, str(value)) for value in value_line]
    str_lines.append('\t'.join(strs))
  return '\n'.join(str_lines)


def format_summary_stats_human(stats):
  output = []
  # Summary stats table.
  if stats['reads'] > 0 or stats['overlap_bp'] > 0:
    output.append('\t\tMinimum\tAverage\tMedian\tMaximum')
  if stats['reads'] > 0:
    line = 'Read lengths:'
    for summary in 'min', 'avg', 'med', 'max':
      line += '\t{}'.format(round(stats[summary+'_rlen'], 3))
    output.append(line)
  if stats['overlap_bp'] > 0:
    line = 'Overlap lens:'
    for summary in 'min', 'avg', 'med', 'max':
      line += '\t{}'.format(round(stats[summary+'_overlap'], 3))
    output.append(line)
    # Start of descriptive text.
    output.append(
      '{error_rate:0.4f} errors per base: {errors} errors in {overlap_bp}bp of overlap.'
      .format(**stats)
    )
  else:
    output.append('No overlaps were detected.')
  if stats['reads'] > 0:
    output.append(
      '{:0.3f}% of reads were in well-mapped pairs: {pairs} pairs out of {reads} total reads.'
      .format(100*stats['paired_read_frac'], **stats)
    )
  else:
    output.append('No reads were found.')
  if stats['pair_bases'] > 0:
    output.append(
      'These pairs contained {pair_bases} bases, {:0.4f}% of which were in overlaps.'
      .format(100*stats['overlap_rate'], **stats)
    )
  else:
    output.append('No paired reads were found.')
  return '\n'.join(output)


def add_computed_stats(stats, counters, precision=6):
  for stat in (
    'error_rate', 'paired_read_frac', 'overlap_rate', 'max_rlen', 'min_rlen', 'avg_rlen', 'med_rlen'
  ):
    stats[stat] = None
  for stat in 'max_overlap', 'min_overlap', 'avg_overlap', 'med_overlap':
    stats[stat] = 0
  if stats['overlap_bp'] > 0:
    stats['error_rate'] = round(stats['errors']/stats['overlap_bp'], precision)
    summarize_list(stats, 'overlap', counters['overlap_lens'])
  if stats['reads'] > 0:
    stats['paired_read_frac'] = round(2*stats['pairs']/stats['reads'], precision)
    summarize_list(stats, 'rlen', counters['read_lens'])
  if stats['pair_bases'] > 0:
    stats['overlap_rate'] = round(2*stats['overlap_bp']/stats['pair_bases'], precision)


def summarize_list(stats, datatype, counter, precision=6):
  summaries = {
    'max': lambda counter: max(counter.keys()),
    'min': lambda counter: min(counter.keys()),
    'avg': get_average,
    'med': get_median,
  }
  for stat_type, summary_fxn in summaries.items():
    value = round(summary_fxn(counter), precision)
    stat_name = f'{stat_type}_{datatype}'
    stats[stat_name] = value


def get_average(counter):
  """Get the average of a series of counts of how often each value occurred."""
  values_total = 0
  count_total = 0
  for value, count in counter.items():
    values_total += value*count
    count_total += count
  if count_total > 0:
    return values_total/count_total
  else:
    return None


def get_median(counter):
  """Get the median from a series of counts of how often each value occurred.
  This honors the strict definition of a median, returning the average of the middle values
  if there's an even number of values (and there are different values on either side on the
  halfway mark)."""
  total = sum(counter.values())
  progress = 0
  last = None
  avg_with_next = False
  for value in sorted(counter.keys()):
    count = counter[value]
    progress += count
    if avg_with_next:
      return (value+last)/2
    elif progress*2 == total:
      # There's an even number of values and we just saw the one on one side of the halfway mark.
      avg_with_next = True
    elif progress*2 >= total+1:
      # We just saw the middle value.
      return value
    last = value


def print_alignment(align_file, detailed_read=None):
  for read in samreader.read(align_file):
    name = mated_name(read)
    line = name+':'
    ref_pos = get_reference_positions(read)
    started = False
    i = 0
    for pos, base in zip(ref_pos, read.seq):
      if name == detailed_read:
        if pos is None:
          print(f'None: {base}')
        else:
          print(f'{pos:4d}: {base}')
      while pos is not None and i < pos:
        if started:
          line += '-'
        else:
          line += ' '
        i += 1
      line += base
      i += 1
      started = True
    print(line)


def format_pair(pair):
  if pair is None:
    return 'None'
  read_strs = []
  for read in pair:
    if read is None:
      read_strs.append('None')
    else:
      read_strs.append(mated_name(read))
  return read_strs


def mated_name(read):
  return f'{read.qname}/{read.mate}'


#TODO: Replace with more efficient implementation in cigarlib itself, if necessary.
#      Maybe take the blocks and compute chunks of coordinates at once.
def get_reference_positions(read):
  positions = []
  readlen = len(read.seq)
  cigar_list = cigarlib.split_cigar(read.cigar)
  blocks = cigarlib.get_contiguous_blocks(read.pos, cigar_list, read.reversed, readlen)
  # This loop takes 96% of the time spent in this function, and 83% of total script time.
  # Note: Performance was measured with time.perf_counter() and accumulating elapsed times in a
  # global dict. These measurement operations took plenty of time themselves: overall, the script
  # took about 1.9x the normal amount of time.
  for read_coord in range(1, readlen+1):
    # This takes 42% of the time spent in this function, and 36% of total script time.
    ref_coord = cigarlib.to_ref_coord(blocks, read_coord)
    # This takes 15% of the time spent in this function, and 13.5% of total script time.
    positions.append(ref_coord)
  if read.reversed:
    return list(reversed(positions))
  else:
    return positions


def human_time(sec):
  """Convert a number of seconds into a human-readable quantity of time.
  Example output: '26.5 minutes', '1.3 months', '5 weeks' (gives an integer if float ends in '.0').
  """
  if sec < 60:
    return format_time(sec, 'second')
  elif sec < 60*60:
    return format_time(sec/60, 'minute')
  elif sec < 24*60*60:
    return format_time(sec/60/60, 'hour')
  elif sec < 10*24*60*60:
    return format_time(sec/60/60/24, 'day')
  elif sec < 40*24*60*60:
    return format_time(sec/60/60/24/7, 'week')
  elif sec < 365*24*60*60:
    return format_time(sec/60/60/24/30.5, 'month')
  else:
    return format_time(sec/60/60/24/365, 'year')


def format_time(quantity, unit):
  rounded = round(quantity, 1)
  if rounded == int(quantity):
    rounded = int(quantity)
  output = f'{rounded} {unit}'
  if rounded != 1:
    output += 's'
  return output


def fail(message):
  logging.critical('Error: '+message)
  if __name__ == '__main__':
    sys.exit(1)
  else:
    raise Exception('Unrecoverable error')


class PairError(RuntimeError):
  pass


class InvalidState(RuntimeError):
  pass


class Pair:

  def __init__(self, first=None, second=None):
    self._pair = [None, None]
    self[0] = first
    self[1] = second
    self._is_well_mapped_cache = None

  @property
  def first(self):
    return self._pair[0]

  @property
  def second(self):
    return self._pair[1]

  @first.setter
  def first(self, value):
    self[0] = value

  @second.setter
  def second(self, value):
    self[1] = value

  @property
  def is_full(self):
    """Are both slots in this pair filled with reads?"""
    if self[0] is None or self[1] is None:
      return False
    return True

  @property
  def num_reads(self):
    num = 0
    if self[0] is not None:
      num += 1
    if self[1] is not None:
      num += 1
    return num

  @property
  def forward(self):
    """Returns True if mate 1 is forward and mate 2 is reverse.
    Returns False if mate1 is reverse and mate 2 is forward.
    Returns None otherwise (including if the pair isn't full yet)."""
    mate1 = mate2 = None
    for read in self:
      if read is None:
        return None
      if read.mate == 1:
        mate1 = read
      elif read.mate == 2:
        mate2 = read
    if mate1 is None or mate2 is None:
      return None
    if mate1.forward and mate2.reverse:
      return True
    elif mate1.reverse and mate2.forward:
      return False

  def __len__(self):
    reads = 0
    for read in self:
      if read is not None:
        reads += 1
    return reads

  def __getitem__(self, index):
    return self._pair[index]

  def __setitem__(self, index, value):
    if value is not None and not isinstance(value, samreader.Alignment):
      raise ValueError(
        f'{type(self).__name__} can only hold None or instances of {samreader.Alignment.__name__}.'
      )
    self._pair[index] = value
    self._is_well_mapped_cache = None

  def __repr__(self):
    return repr(self._pair)

  def __str__(self):
    return str(self._pair)

  def get(self, property, throw=True):
    """Get the value of a SAM column for both reads in the pair.
    Only valid for columns that could be expected to be identical for both reads:
    `qname` and `rname`.
    If the pair isn't full or the values disagree between the mates, this will either return `None`
    (if `throw` is `False`) or throw a `PairError` (if `throw` is `True`)."""
    if property not in ('qname', 'rname'):
      raise ValueError(f'Invalid `property` {property!r}')
    if not self.is_full:
      if throw:
        raise PairError(f'Pair only contains {len(self)} reads.')
      else:
        return None
    value1 = getattr(self[0], property)
    value2 = getattr(self[1], property)
    if value1 != value2:
      if throw:
        raise PairError(f'Values of property {property!r} differ: {value1!r} vs {value2!r}')
      else:
        return None
    return value1

  def add(self, value):
    if self[1] is None:
      if self[0] is None:
        self[0] = value
      else:
        self[1] = value
    else:
      raise InvalidState('Pair is already full.')

  def shift(self):
    removed = self._pair[0]
    self._pair = [self._pair[1], None]
    self._is_well_mapped_cache = None
    return removed

  def is_well_mapped(self, mapq_thres=None, cached=False):
    if self._is_well_mapped_cache is None:
      if cached:
        logging.warning('Warning: Cached value of is_well_mapped requested, but none is cached.')
        return None
      self._is_well_mapped_cache = self._is_well_mapped(mapq_thres=mapq_thres)
    return self._is_well_mapped_cache

  def _is_well_mapped(self, mapq_thres=None):
    for read in self:
      if not self._read_is_well_mapped(read, mapq_thres=mapq_thres):
        return False
    return True

  GOOD_FLAGS = 1 | 2
  BAD_FLAGS = 4 | 8 | 256 | 512 | 1024 | 2048
  @classmethod
  def _read_is_well_mapped(cls, read, mapq_thres=None):
    if read is None:
      return False
    if read.rnext is None:
      return False
    if read.flag & cls.GOOD_FLAGS != cls.GOOD_FLAGS:
      return False
    if read.flag & cls.BAD_FLAGS:
      return False
    if mapq_thres is not None and read.mapq < mapq_thres:
      return False
    return True


class Error(typing.NamedTuple):
  type: str
  rname: str
  ref: str
  ref_coord: int
  read_coord1: int
  read_coord2: int
  base1: str
  base2: str
  qual1: str
  qual2: str


if __name__ == '__main__':
  try:
    sys.exit(main(sys.argv))
  except BrokenPipeError:
    pass
