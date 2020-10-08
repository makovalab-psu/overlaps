#!/usr/bin/env python3
import argparse
import collections
import csv
import gzip
import logging
import math
import pathlib
import sys
from typing import Dict, List
assert sys.version_info.major >= 3, 'Python 3 required'

DESCRIPTION = """Produce statistics on the sequence context surrounding every error in a sample.
This will count how often each base appears at each position relative to the site of each error."""

# Note: This is based on code from Jupyter notebooks, especially
# 2019-12-17-ecoli-auto-seq-context-debug.ipynb.


def make_argparser():
  parser = argparse.ArgumentParser(add_help=False, description=DESCRIPTION)
  io = parser.add_argument_group('Options')
  io.add_argument('errors_path', metavar='errors.tsv', type=pathlib.Path,
    help='File containing data on intra-pair errors.')
  io.add_argument('context_path', metavar='seq-context.tsv', type=pathlib.Path,
    help='File containing data on reference sequence context near errors.')
  io.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
    help='Destination file for output. Warning: Any existing file will be overwritten. '
      'Default: stdout')
  options = parser.add_argument_group('Options')
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

  counts_by_substitution = get_base_counts_by_substitution(args.errors_path, args.context_path)

  for data in make_table(counts_by_substitution):
    print(*data, sep='\t', file=args.output)


def get_base_counts_by_substitution(
    errors_path: pathlib.Path, context_path: pathlib.Path, before_bases=None
  ) -> Dict[str,Dict[int,Dict[str,int]]]:
  """Returns a mapping of:
    substitution -> [context coord -> [context base -> count]]
  e.g.
    counts_by_error[('G','T')][-7]['A'] == 966400
  The substitution format is a tuple of (before_base, after_base).
  Filter by providing a list (or any iterable) of `before_bases`.
  This will skip any substitution where the starting base doesn't appear in `before_bases`.
  For errors where neither read agrees with the reference, the substitution is `None`."""
  counts_by_error = collections.defaultdict(lambda: collections.defaultdict(collections.Counter))
  all_read_bases = get_error_bases(errors_path)
  # Count the bases at every site affected by errors.
  unique_contexts = set(read_context(context_path))
  for context in unique_contexts:
    read_bases = all_read_bases[context.chr][context.pos]
    sample_base = call_base(read_bases, context.base)
    # Count votes for the most common error base that `sample_base` got mutated into.
    error_base_counts = collections.Counter()
    for base1, base2 in read_bases:
      error_base = get_error_base(sample_base, base1, base2)
      if before_bases is not None and (sample_base not in before_bases):
        continue
      error_base_counts[error_base] += 1
    if not error_base_counts:
      continue
    most_common_error_base = error_base_counts.most_common(1)[0][0]
    substitution = (sample_base, most_common_error_base)
    # Count up the reference bases around this site.
    i = context.index
    for j, base in enumerate(context.seq):
      rel_pos = j - i
      counts_by_error[substitution][rel_pos][base] += 1
  return counts_by_error


def tuple_eq(tup1, tup2):
  for attr in dir(tup1):
    if attr.startswith('_') or attr == 'count':
      continue
    if getattr(tup1, attr) != getattr(tup2, attr):
      return False
  return True


Context = collections.namedtuple('Context', ('chr', 'pos', 'index', 'base', 'seq', 'gc'))
Context.__eq__ = tuple_eq


def read_context(context_path):
  with context_path.open() as context_file:
    for line_raw in context_file:
      fields = line_raw.rstrip('\r\n').split('\t')
      if len(fields) != 6:
        continue
      if fields[5] == '.':
        gc = None
      else:
        gc = float(fields[5])
      yield Context(
        chr=fields[0],
        pos=int(fields[1]),
        index=int(fields[2]),
        base=fields[3],
        seq=fields[4],
        gc=gc,
      )


def read_errors(errors_path):
  pair_line = None
  with errors_path.open() as errors_file:
    for line_raw in errors_file:
      if line_raw.startswith('pair\t'):
        pair_line = line_raw
      elif line_raw.startswith('error\t'):
        pair_fields = pair_line.rstrip('\r\n').split('\t')
        pair = parse_fields(pair_fields[1:], Pair, PAIR_METADATA)
        error_fields = line_raw.rstrip('\r\n').split('\t')
        error = parse_fields(error_fields[1:], Error, ERROR_METADATA, pair=pair)
        yield error
      else:
        raise ValueError(f'Invalid first field in line {line_raw!r}')


def bool_literal(raw_value):
  if raw_value == 'True':
    return True
  elif raw_value == 'False':
    return False
  elif raw_value == 'None':
    return None
  else:
    raise ValueError(f'Invalid bool literal {raw_value}')


ERROR_METADATA = [
  {'key':'type', 'type':str},
  {'key':'chr', 'type':str},
  {'key':'pos', 'type':int},
  {'key':'read_coord1', 'type':int},
  {'key':'read_coord2', 'type':int},
  {'key':'base1', 'type':str},
  {'key':'base2', 'type':str},
]
PAIR_METADATA = [
  {'key':'name', 'type':str},
  {'key':'mapped', 'type':bool_literal},
  {'key':'rlen1', 'type':int},
  {'key':'rlen2', 'type':int},
  {'key':'overlap', 'type':int},
  {'key':'errors', 'type':int},
]

Error = collections.namedtuple('Error', [meta['key'] for meta in ERROR_METADATA]+['pair'])
Pair = collections.namedtuple('Pair', [meta['key'] for meta in PAIR_METADATA])
Error.__eq__ = Pair.__eq__ = tuple_eq

def parse_fields(raw_fields, tuple_type, metadata, null='.', **kwargs):
  data = {}
  for raw_value, meta in zip(raw_fields, metadata):
    key = meta['key']
    if raw_value == null:
      data[key] = None
    else:
      data[key] = meta['type'](raw_value)
  for key, value in kwargs.items():
    data[key] = value
  return tuple_type(**data)


def get_error_bases(errors_path):
  """Returns all the read alleles for all errors in the format:
  error_bases[ref][coord] == error_list
  error_list == [('A', 'C'), ('G', 'C'), ...]"""
  error_bases = {}
  for error in read_errors(errors_path):
    bases_by_coord = error_bases.get(error.chr, {})
    bases_list = bases_by_coord.get(error.pos, [])
    bases_list.append((error.base1, error.base2))
    bases_by_coord[error.pos] = bases_list
    error_bases[error.chr] = bases_by_coord
  return error_bases


def get_error_base(ref_base, read_base1, read_base2):
  """Decide which is the error out of two disagreeing bases.
  This takes two discordant bases (`read_base1` and `read_base2`) from a pair of reads, plus the
  reference base at that position `ref_base`, and determines which is the error base. It assumes
  the read base that matches the reference base is the original, "correct" base, and the other read
  base is the error. If neither matches the reference base, it returns `None`.
  `ref_base` can also be the sample base instead of the actual base from the reference sequence."""
  if read_base1 == ref_base:
    return read_base2
  elif read_base2 == ref_base:
    return read_base1
  else:
    return None


def call_base(error_list, ref_base=None):
  """Do (naive) variant calling at this site.
  Provide a series of bases and this will return the most common one.
  In the case of a tie, it will pick arbitrarily. Or, if a `ref_base` is given, and it's among the
  bases tied for first place, it will pick that.
  `error_list`: A sequence of 2-tuples. Each tuple is normally a pair of disagreeing bases from one
    location in a pair of reads.
  `ref_base`:   The genomic base at that position from the reference file. If given, only used to
    break a tie."""
  if not error_list:
    return ref_base
  votes = collections.Counter()
  for base1, base2 in error_list:
    votes[base1] += 1
    votes[base2] += 1
  # Find the winning base. If there's a tie, and the ref_base is among the winners, pick that.
  top_count = None
  top_bases = []
  for base, count in votes.most_common():
    if top_count is None or count == top_count:
      # This is the winner or tied for winner.
      top_count = count
      top_bases.append(base)
    elif count < top_count:
      # We've finished going through the winners.
      break
    else:
      raise RuntimeError(f'Invalid state: count > top_count ({count} > {top_count})')
  # Return the winner.
  if len(top_bases) == 1:
    return top_bases[0]
  elif len(top_bases) > 1:
    # It's a tie (multiple winners).
    if ref_base in top_bases:
      # Prefer the reference base, if it's among the winners.
      return ref_base
    else:
      # If the reference base isn't among the winners, just pick any of them.
      return top_bases[0]
  else:
    raise RuntimeError(f'Invalid state: len(top_bases) < 1 ({top_bases})')


def make_table(counts_by_substitution):
  for (before_base, after_base), counts_by_pos in counts_by_substitution.items():
    for pos, counts_by_base in counts_by_pos.items():
      for ref_base, count in counts_by_base.items():
        yield (before_base, after_base, pos, ref_base, count)


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
