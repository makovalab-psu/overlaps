#!/usr/bin/env python3
import argparse
import collections
import logging
import pathlib
import subprocess
import sys
import time
from bfx import cigarlib
from bfx import samreader
assert sys.version_info.major >= 3, 'Python 3 required'

DESCRIPTION = """"""

class Error(object):
  __slots__ = ('type', 'rname', 'ref_coord', 'coord1', 'coord2', 'alt1', 'alt2')
  defaults = {'type':'snv'}
  def __init__(self, **kwargs):
    for name in self.__slots__:
      if name in kwargs:
        setattr(self, name, kwargs[name])
      else:
        default = self.defaults.get(name, None)
        setattr(self, name, default)


def make_argparser():
  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.add_argument('align', type=pathlib.Path, nargs='?', default=sys.stdin,
    help='Name-sorted BAM or SAM file. Omit to read from stdin.')
  parser.add_argument('-q', '--mapq', type=int,
    help='Mapping quality threshold. Alignments worse than this MAPQ will be ignored.')
  parser.add_argument('-H', '--human', dest='format', action='store_const', const='human',
    default='tsv',
    help='Print human-readable text instead of tab-delimited output.')
  parser.add_argument('-s', '--summary', action='store_true',
    help='Print summary stats for the entire file. This is the default, unless --details is given. '
      'Giving --summary AND --details forces both to be printed.')
  parser.add_argument('-d', '--details', action='store_true',
    help='Print detailed info for each read pair.')
  parser.add_argument('-p', '--progress', type=int, default=15,
    help='Print periodic updates to stderr while processing large files. Give a number X and this '
      'will print human-readable stats every X minutes (as well as an initial update X/10 minutes '
      'in). Give 0 to turn off these messages. Default: %(default)s minutes')
  parser.add_argument('-P', '--print', action='store_true',
    help='Debug print the alignment as cigarlib sees it.')
  parser.add_argument('-R', '--detailed-read',
    help='When using --print, single out this read and print every reference coordinate of every '
      'base. Give the full read name, ending in a /1 or /2 for the mate.')
  parser.add_argument('-l', '--log', type=argparse.FileType('w'), default=sys.stderr,
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  volume = parser.add_mutually_exclusive_group()
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
  align_file = open_input(args.align)
  summary = True
  details = False
  if args.details:
    details = True
    if not args.summary:
      summary = False

  if args.print:
    print_alignment(align_file, args.detailed_read)
    return

  overlapper = Overlapper(align_file)
  process_file(overlapper, args.mapq, args.format, details, args.progress)

  if summary:
    if args.progress and args.format == 'human':
      print('Finished!', file=sys.stderr)
    print(format_summary_stats(overlapper.stats, overlapper.counters, args.format))


def process_file(overlapper, mapq_thres, format, details, progress):
  progress_sec = progress*60
  start = int(time.time())
  last = None
  for errors, pair, overlap_len in overlapper.analyze_overlaps(mapq_thres):
    if details:
      print(format_read_stats(errors, pair, overlap_len, format=format))
    else:
      is_progress_time, last = is_it_progress_time(progress_sec, start, last)
      if is_progress_time:
        now = int(time.time())
        output = (
          '\tProcessed {} reads after {}:\n'
          .format(overlapper.stats['reads'], human_time(now-start))
        )
        output += format_summary_stats(overlapper.stats, overlapper.counters, format)
        print(output, file=sys.stderr)


def open_input(align_path):
  if align_path is sys.stdin:
    return align_path
  elif not align_path.is_file():
    fail(f'Error: Input file must be a regular file. {str(align_path)!r} is not.')
  if align_path is None or str(align_path) == '-':
    return sys.stdin
  else:
    if align_path.suffix == '.bam':
      return open_bam(align_path)
    else:
      return open(align_path)


def open_bam(bam_path):
  process = subprocess.Popen(('samtools', 'view', bam_path), stdout=subprocess.PIPE)
  for line in process.stdout:
    yield str(line, 'utf8')


class Overlapper:

  def __init__(self, align_file):
    self.align_file = align_file
    self.stats = {'reads':0, 'pairs':0, 'pair_bases':0, 'overlap_bp':0, 'errors':0}
    self.counters = {'overlap_lens':collections.Counter(), 'read_lens':collections.Counter()}

  def analyze_overlaps(self, mapq_thres=None):
    for pair in self.get_pairs(mapq_thres=mapq_thres):
      errors, overlap_len = get_mismatches(pair)
      self.stats['pair_bases'] += len(pair[0].seq) + len(pair[1].seq)
      self.stats['overlap_bp'] += overlap_len
      self.stats['errors'] += len(errors)
      self.stats['pairs'] += 1
      self.counters['overlap_lens'][overlap_len] += 1
      yield errors, pair, overlap_len


  def get_pairs(self, mapq_thres=None):
    pair = [None, None]
    last_name = None
    for read in samreader.read(self.align_file):
      self.stats['reads'] += 1
      self.counters['read_lens'][read.length] += 1
      # logging.debug(mated_name(read))
      name = read.qname
      assert not (name.endswith('/1') or name.endswith('/2')), name
      mate_num = which_mate(read)
      assert mate_num in (1, 2), read.flag
      # logging.debug(f'pair: {format_pair(pair)}')
      if pair[0] is None:
        pair[0] = read
        continue
      elif pair[1] is None:
        pair[1] = read
      else:
        raise fail(
          'Error: Invalid state. Encountered a full pair too early. Pair:\n  {}\n  {}'
          .format(pair[0].qname, pair[1].qname)
        )
      # We should only have full pairs at this point.
      # If the names don't match, the first read can't have a mate in the file, since this is name-
      # sorted. Discard it, shift the second read to the first position, and loop again.
      if pair[0].qname != pair[1].qname:
        pair = [pair[1], None]
        continue
      if pair_is_well_mapped(pair, mapq_thres):
        yield pair
      pair = [None, None]


def pair_is_well_mapped(pair, mapq_thres=None):
  read1, read2 = pair
  if not(read1.flag & 2 and read2.flag & 2):
    return False
  if read1.rnext is None or read2.rnext is None:
    return False
  if mapq_thres is not None and (read1.mapq < mapq_thres or read2.mapq < mapq_thres):
    return False
  return True


def get_mismatches(pair):
  errors = []
  overlap_len = 0
  read1, read2 = pair
  base_map = get_base_map(read2)
  positions = get_reference_positions(read1)
  started = False
  for pos, base1 in zip(positions, read1.seq):
    base2 = base_map.get(pos)
    if pos is None:
      pass #logging.debug(f'None: {base1} -> <INS>')
    elif base2 is None:
      if started:
        pass #logging.debug(f'{pos:4d}: {base1} -> <DEL>')
        #TODO: Oops, insertions don't have a reference position. How to detect them?
        #      Doesn't look like pysam gives an easy way. Looks like more CIGAR parsing.
    else:
      started = True
      overlap_len += 1
      if base1 != base2:
        pass #logging.debug(f'{pos:4d}: {base1} -> {base2}')
        error = Error(type='snv', rname=read1.qname, ref_coord=pos, alt1=base1, alt2=base2)
        #TODO: Get the read coordinates (coord1, coord2).
        errors.append(error)
  return errors, overlap_len


def get_base_map(read):
  """Return a mapping of reference coordinates to read bases."""
  positions = get_reference_positions(read)
  return dict(zip(positions, read.seq))


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


def format_read_stats(errors, pair, overlap_len, format='tsv'):
  if format == 'tsv':
    return format_read_stats_tsv(errors, pair, overlap_len)
  elif format == 'human':
    return format_read_stats_human(errors, pair, overlap_len)


def format_read_stats_tsv(errors, pair, overlap_len):
  fields = [pair[0].qname, len(errors), overlap_len]
  for error in errors:
    fields.append(error.ref_coord)
  return '\t'.join(map(str, fields))


def format_read_stats_human(errors, pair, overlap_len):
  output = f'{pair[0].qname}: {len(errors)} errors in {overlap_len}bp'
  if errors:
    output += ':'
  for error in errors:
    output += f'\n{error.ref_coord:4d} {error.type.upper()}: {error.alt1} -> {error.alt2}'
  return output


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
  output = []
  for line_stats in output_stats:
    line = '\t'.join([str(stats[stat]) for stat in line_stats])
    output.append(line)
  return '\n'.join(output)


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
    stats['error_rate'] = round(stats['errors']/stats['overlap_bp'], 6)
    summarize_list(stats, 'overlap', counters['overlap_lens'])
  if stats['reads'] > 0:
    stats['paired_read_frac'] = round(2*stats['pairs']/stats['reads'], 6)
    summarize_list(stats, 'rlen', counters['read_lens'])
  if stats['pair_bases'] > 0:
    stats['overlap_rate'] = round(2*stats['overlap_bp']/stats['pair_bases'], 6)


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
  return f'{read.qname}/{which_mate(read)}'


#TODO: Replace with more efficient implementation in cigarlib itself, if necessary.
#      Maybe take the blocks and compute chunks of coordinates at once.
def get_reference_positions(read):
  positions = []
  readlen = len(read.seq)
  cigar_list = cigarlib.split_cigar(read.cigar)
  reverse = read_is_reversed(read)
  blocks = cigarlib.get_contiguous_blocks(read.pos, cigar_list, reverse, readlen)
  # This loop takes 96% of the time spent in this function, and 83% of total script time.
  # Note: Performance was measured with time.perf_counter() and accumulating elapsed times in a
  # global dict. These measurement operations took plenty of time themselves: overall, the script
  # took about 1.9x the normal amount of time.
  for read_coord in range(1, readlen+1):
    # This takes 42% of the time spent in this function, and 36% of total script time.
    ref_coord = cigarlib.to_ref_coord(blocks, read_coord)
    # This takes 15% of the time spent in this function, and 13.5% of total script time.
    positions.append(ref_coord)
  if reverse:
    return list(reversed(positions))
  else:
    return positions


def which_mate(read):
  if read.flag & 64:
    return 1
  elif read.flag & 128:
    return 2


def read_is_reversed(read):
  if read.flag & 16:
    return True
  else:
    return False


def mate_is_mapped(read):
  return read.flag & 2


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
  logging.critical(message)
  if __name__ == '__main__':
    sys.exit(1)
  else:
    raise Exception('Unrecoverable error')


if __name__ == '__main__':
  try:
    sys.exit(main(sys.argv))
  except BrokenPipeError:
    pass
