#!/usr/bin/env python3
import argparse
import logging
import pathlib
import sys
# Path hack to load modules from the parent directory.
script_dir = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(script_dir.parent))
import read_formats

DESCRIPTION = """"""


def make_argparser():
  parser = argparse.ArgumentParser(add_help=False, description=DESCRIPTION)
  options = parser.add_argument_group('Options')
  options.add_argument('function', choices=('read_errors_file','read_errors'))
  options.add_argument('infile', type=argparse.FileType('r'))
  options.add_argument('-n', '--null', default='.',
    help='The string to be considered None. Default: %(default)s')
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

  if args.function in ('read_errors_file', 'read_errors'):
    read_errors(args.function, args.infile, args.null)


def read_errors(function, infile, null):
    if function == 'read_errors_file':
      reader = read_formats.read_errors_file
    elif function == 'read_errors':
      reader = read_formats.read_errors
    for pair in reader(infile, null=null):
      pair_simplified = omit_fields(pair, ('errors','overlap'))
      print(pair_simplified)
      print(' ', pair.overlap)
      for error in pair.errors:
        print(' ', error)


def omit_fields(value, omissions):
  kwargs = {}
  for attr in omissions:
    kwargs[attr] = None
  return value._replace(**kwargs)


def fail(message):
  logging.critical(f'Error: {message}')
  if __name__ == '__main__':
    sys.exit(1)
  else:
    raise Exception(message)


if __name__ == '__main__':
  try:
    sys.exit(main(sys.argv))
  except BrokenPipeError:
    pass
