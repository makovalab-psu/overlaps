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
  parser.add_argument('-p', '--print', action='store_true',
    help='Debug print the alignment as cigarlib sees it.')
  parser.add_argument('-l', '--log', type=argparse.FileType('w'), default=sys.stderr,
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  volume = parser.add_mutually_exclusive_group()
  volume.add_argument('-q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
    default=logging.WARNING)
  volume.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  volume.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)
  return parser


def main(argv):

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')

  if args.align:
    if args.align.suffix == '.bam':
      align_file = open_bam(args.align)
    else:
      align_file = open(args.align)
  else:
    align_file = sys.stdin

  if args.print:
    print_alignment(align_file)
    return

  for pair in get_pairs(align_file):
    print(pair[0].qname+':')
    get_mismatches(pair)


def open_bam(bam_path):
  process = subprocess.Popen(('samtools', 'view', bam_path), stdout=subprocess.PIPE)
  for line in process.stdout:
    yield str(line, 'utf8')


def get_pairs(align_file):
  pair = None
  last_name = None
  for read in samreader.read(align_file):
    logging.debug(mated_name(read))
    name = read.qname
    assert not (name.endswith('/1') or name.endswith('/2')), name
    mate_num = which_mate(read)
    assert mate_num in (1, 2), read.flag
    logging.debug(f'pair: {format_pair(pair)}')
    #TODO: Rely on something other than flag 2. Pairs can be both mapped, but not in their proper
    #      pair. Also, you might eventually want to filter by MAPQ.
    #      Best to gather the pair if they both have the same name, then judge later whether they're
    #      both mapped acceptably.
    if pair is None:
      if mate_is_mapped(read):
        if mate_num == 1:
          pair = [read, None]
        elif mate_num == 2:
          pair = [None, read]
    else:
      if pair[0] is None:
        logging.debug(f'  setting pair[0] = {mated_name(read)}')
        pair[0] = read
      elif pair[1] is None:
        logging.debug(f'  setting pair[1] = {mated_name(read)}')
        pair[1] = read
      else:
        fail('Error: Invalid pair.')
      validate_pair(pair)
      yield pair
      pair = None


def validate_pair(pair):
  read1, read2 = pair
  if read1.qname != read2.qname:
    mate_states = []
    for read in pair:
      if mate_is_mapped(read):
        mate_states.append('mate mapped')
      else:
        mate_states.append('mate unmapped')
    fail(
      f'Error: Sequential reads have different names: '
      f'{read1.qname!r} ({mate_states[0]}) != '
      f'{read2.qname!r} ({mate_states[1]})'
    )
  if which_mate(read1) != 1 or which_mate(read2) != 2:
    fail(
      f'Error: Mates were paired up incorrectly: mate 1 != {which_mate(read1)} or '
      f'mate 2 != {which_mate(read2)}'
    )


def get_mismatches(pair):
  read1, read2 = pair
  base_map = get_base_map(read2)
  positions = get_reference_positions(read1)
  started = False
  for pos, base1 in zip(positions, read1.seq):
    base2 = base_map.get(pos)
    if base2 is None:
      if started:
        print(f'  {pos:2d}: {base1} -> <DEL>')
        #TODO: Oops, insertions don't have a reference position. How to detect them?
        #      Doesn't look like pysam gives an easy way. Looks like more CIGAR parsing?
    else:
      started = True
      if base1 != base2:
        print(f'  {pos:2d}: {base1} -> {base2}')
  # 'type', 'rname', 'ref_coord', 'coord1', 'coord2', 'alt1', 'alt2'


def get_base_map(read):
  """Return a mapping of reference coordinates to read bases."""
  positions = get_reference_positions(read)
  return dict(zip(positions, read.seq))


def print_alignment(align_file):
  for read in samreader.read(align_file):
    name = mated_name(read)
    line = name+':'
    ref_pos = get_reference_positions(read)
    started = False
    i = 0
    for pos, base in zip(ref_pos, read.seq):
      if name == 'BB/2':
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
