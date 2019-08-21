#!/usr/bin/env python3
import argparse
import logging
import pathlib
import sys
import pysam
assert sys.version_info.major >= 3, 'Python 3 required'

DESCRIPTION = """"""


def make_argparser():
  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.add_argument('bam', type=pathlib.Path,
    help='')
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

  print_alignment(args.bam)


def print_alignment(bam_path):
  bam = pysam.AlignmentFile(bam_path, 'rb')
  for read in bam.fetch():
    line = f'{read.query_name}/{which_mate(read.flag)}:'
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


def which_mate(flag):
  if flag & 64:
    return 1
  elif flag & 128:
    return 2


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
