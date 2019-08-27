#!/usr/bin/env python3
import argparse
import logging
import pathlib
import subprocess
import sys
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
  parser.add_argument('align', type=pathlib.Path,
    help='Name-sorted BAM or SAM file.')
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
  parser.add_argument('-p', '--print', action='store_true',
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

  stats = print_andor_compile_stats(align_file, args.format, details, args.mapq)

  if summary:
    print(format_summary_stats(stats, args.format))


def print_andor_compile_stats(align_file, format, details, mapq_thres):
  stats = {'reads':0, 'pairs':0, 'bases':0, 'overlap':0, 'errors':0}
  for pair in get_pairs(align_file, stats, mapq_thres):
    errors, overlap_len = get_mismatches(pair)
    stats['bases'] += len(pair[0].seq) + len(pair[1].seq)
    stats['overlap'] += overlap_len
    stats['errors'] += len(errors)
    stats['pairs'] += 1
    if details:
      print(format_read_stats(errors, pair, overlap_len, format=format))
  return stats


def open_input(align_path):
  if align_path:
    if align_path.suffix == '.bam':
      return open_bam(align_path)
    else:
      return open(align_path)
  else:
    return sys.stdin


def open_bam(bam_path):
  process = subprocess.Popen(('samtools', 'view', bam_path), stdout=subprocess.PIPE)
  for line in process.stdout:
    yield str(line, 'utf8')


def get_pairs(align_file, stats, mapq_thres=None):
  pair = [None, None]
  last_name = None
  for read in samreader.read(align_file):
    stats['reads'] += 1
    logging.debug(mated_name(read))
    name = read.qname
    assert not (name.endswith('/1') or name.endswith('/2')), name
    mate_num = which_mate(read)
    assert mate_num in (1, 2), read.flag
    logging.debug(f'pair: {format_pair(pair)}')
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
      logging.debug(f'None: {base1} -> <INS>')
    elif base2 is None:
      if started:
        logging.debug(f'{pos:4d}: {base1} -> <DEL>')
        #TODO: Oops, insertions don't have a reference position. How to detect them?
        #      Doesn't look like pysam gives an easy way. Looks like more CIGAR parsing.
    else:
      started = True
      overlap_len += 1
      if base1 != base2:
        logging.debug(f'{pos:4d}: {base1} -> {base2}')
        error = Error(type='snv', rname=read1.qname, ref_coord=pos, alt1=base1, alt2=base2)
        #TODO: Get the read coordinates (coord1, coord2).
        errors.append(error)
  return errors, overlap_len


def get_base_map(read):
  """Return a mapping of reference coordinates to read bases."""
  positions = get_reference_positions(read)
  return dict(zip(positions, read.seq))


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


def format_summary_stats(stats, format='tsv'):
  if format == 'tsv':
    return (
      '{}\t{errors}\t{overlap}\t{pairs}\t{reads}\t{bases}'
      .format(round(stats['errors']/stats['overlap'], 6), **stats)
    )
  elif format == 'human':
    return (
      '{} errors per base: {errors} errors in {overlap}bp of overlap.\n'
      '{reads} total reads in {pairs} well-mapped pairs totalling {bases} bases (in the pairs).'
      .format(round(stats['errors']/stats['overlap'], 4), **stats)
    )


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
#      Take the blocks and compute chunks of coordinates at once.
def get_reference_positions(read):
  positions = []
  readlen = len(read.seq)
  cigar_list = cigarlib.split_cigar(read.cigar)
  reverse = read_is_reversed(read)
  blocks = cigarlib.get_contiguous_blocks(read.pos, cigar_list, reverse, readlen)
  for read_coord in range(1, readlen+1):
    ref_coord = cigarlib.to_ref_coord(blocks, read_coord)
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
