#!/usr/bin/env python3
import argparse
import logging
import pathlib
import sys
import pysam
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
  parser.add_argument('bam', type=pathlib.Path,
    help='Name-sorted BAM file.')
  parser.add_argument('-p', '--print', action='store_true',
    help='Debug print the alignment as pysam sees it.')
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

  if args.print:
    print_alignment(args.bam)
    return

  for pair in get_pairs(args.bam):
    print(pair[0].query_name+':')
    get_mismatches(pair)


def get_pairs(bam_path):
  bam = pysam.AlignmentFile(bam_path, 'rb')
  pair = None
  last_name = None
  for read in bam.fetch(until_eof=True):
    logging.debug(mated_name(read))
    name = read.query_name
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
  if read1.query_name != read2.query_name:
    mate_states = []
    for read in pair:
      if mate_is_mapped(read):
        mate_states.append('mate mapped')
      else:
        mate_states.append('mate unmapped')
    fail(
      f'Error: Sequential reads have different names: '
      f'{read1.query_name!r} ({mate_states[0]}) != '
      f'{read2.query_name!r} ({mate_states[1]})'
    )
  if which_mate(read1) != 1 or which_mate(read2) != 2:
    fail(
      f'Error: Mates were paired up incorrectly: mate 1 != {which_mate(read1)} or '
      f'mate 2 != {which_mate(read2)}'
    )


def get_mismatches(pair):
  read1, read2 = pair
  base_map = get_base_map(read2)
  positions = read1.get_reference_positions()
  sequence = read1.query_sequence
  started = False
  for pos, base1 in zip(positions, sequence):
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


def get_base_map(read):
  """Return a mapping of reference coordinates to read bases."""
  positions = read.get_reference_positions()
  sequence = read.query_sequence
  return dict(zip(positions, sequence))


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
  return f'{read.query_name}/{which_mate(read)}'


def print_alignment(bam_path):
  bam = pysam.AlignmentFile(bam_path, 'rb')
  for read in bam.fetch(until_eof=True):
    line = f'{read.query_name}/{which_mate(read)}:'
    seq = read.query_sequence
    ref_pos = read.get_reference_positions()
    started = False
    i = 0
    for pos, base in zip(ref_pos, seq):
      while i < pos:
        if started:
          line += '-'
        else:
          line += ' '
        i += 1
      line += base
      i += 1
      started = True
    print(line)


def which_mate(read):
  if read.flag & 64:
    return 1
  elif read.flag & 128:
    return 2


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
