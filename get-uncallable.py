#!/usr/bin/env python3
import argparse
import importlib
import logging
import pathlib
import sys
import time
assert sys.version_info.major >= 3, 'Python 3 required'
summarize_context = importlib.import_module('summarize-context')

DESCRIPTION = """Calculate how many substitutions whose direction we can't determine.
E.g. the reference doesn't match either of the read bases."""


def make_argparser():
  parser = argparse.ArgumentParser(add_help=False, description=DESCRIPTION)
  options = parser.add_argument_group('Options')
  options.add_argument('errors', metavar='errors.tsv', type=pathlib.Path)
  options.add_argument('contexts', metavar='seq-context.tsv', type=pathlib.Path)
  options.add_argument('outfile', type=argparse.FileType('w'), default=sys.stdout, nargs='?')
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

  start = time.perf_counter()
  called, uncalled = get_callable_counts(args.errors, args.contexts)
  elapsed = round(time.perf_counter() - start)
  if called == 0:
    pct = '.'
  else:
    pct = round(100*uncalled/called, 7)
  print(called, uncalled, pct, elapsed, sep='\t', file=args.outfile)


def get_callable_counts(errors_path, contexts_path):
  called = uncalled = 0
  all_read_bases = summarize_context.get_error_bases(errors_path)
  unique_contexts = set(summarize_context.read_context(contexts_path))
  for context in unique_contexts:
    read_bases = all_read_bases[context.chr][context.pos]
    sample_base = summarize_context.call_base(read_bases, context.base)
    for base1, base2 in read_bases:
      error_base = summarize_context.get_error_base(sample_base, base1, base2)
      if error_base is None:
        uncalled += 1
      else:
        called += 1
  return called, uncalled


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
