#!/usr/bin/env python3
import argparse
import collections
import logging
import sys
assert sys.version_info.major >= 3, 'Python 3 required'

# Adapted from 2020-01-28-ecoli-auto-post-homopol2.ipynb

DESCRIPTION = """Sum up counts of post-homopolymer errors."""


def make_argparser():
  parser = argparse.ArgumentParser(add_help=False, description=DESCRIPTION)
  options = parser.add_argument_group('Options')
  options.add_argument('contexts', metavar='seq-context.tsv', type=argparse.FileType('r'))
  options.add_argument('errors', metavar='errors.tsv', type=argparse.FileType('r'))
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

  hpol_counts = count_hpol_errors(args.contexts, args.errors)

  for line in format_table(sorted(make_table(hpol_counts), key=get_comparable_list)):
    print(line)


def count_hpol_errors(contexts_file, errors_file):
  hpol_counts = collections.defaultdict(collections.Counter)
  context_index = compile_contexts(read_context(contexts_file))
  for error in read_errors(errors_file):
    context = context_index[error.chr][error.pos]
    run_base, run_len = get_prerun(context.seq, context.index)
    if run_base is None:
      continue
    error_base = get_error_base(context.base, error.base1, error.base2)
    if error_base == run_base:
      hpol_counts[error_base][run_len] += 1
    else:
      hpol_counts[error_base][None] += 1
  return hpol_counts


def compile_contexts(contexts):
  contexts_by_chr = {}
  duplicates = 0
  for context in contexts:
    contexts_by_pos = contexts_by_chr.get(context.chr, {})
    if context.pos in contexts_by_pos:
      other = contexts_by_pos[context.pos]
      if context == other:
        duplicates += 1
        continue
      else:
        raise ValueError(f'Different contexts for pos {context.pos}:\n{context}\n{other}')
    contexts_by_pos[context.pos] = context
    contexts_by_chr[context.chr] = contexts_by_pos
  return contexts_by_chr


def get_prerun(seq, index):
  """Look for homopolymers preceding the error base.
  Arguments should be the `seq` and `index` attributes of a `Context`.
  Returns a tuple of `run_base, run_len`:
  `run_base`: The base preceding the error position.
  `run_len`:  How many times that base occurred in a row right before the error.
  Returns `(None, 0)` if the index is 0 (the error base is the first base in the seq)."""
  run_len = 0
  run_base = None
  if not 0 <= index < len(seq):
    raise ValueError(f'index {index!r} out of bounts for seq {seq!r}.')
  for i, char in enumerate(seq):
    if index == i:
      return run_base, run_len
    if char == run_base:
      run_len += 1
    else:
      run_len = 1
      run_base = char


def get_error_base(ref_base, error_base1, error_base2):
  if error_base1 == ref_base:
    return error_base2
  elif error_base2 == ref_base:
    return error_base1
  else:
    return None


def get_comparable_list(iterable):
  return [get_comparable_value(value) for value in iterable]


def get_comparable_value(raw_value):
  if raw_value is None:
    is_none = 0
    value = -1
  else:
    is_none = 1
    value = raw_value
  return (is_none, value)


def make_table(hpol_counts):
  for error_base, run_len_counts in hpol_counts.items():
    for run_len, count in run_len_counts.items():
      yield (error_base, run_len, count)


def format_table(rows):
  for raw_row in rows:
    row = []
    for raw_value in raw_row:
      if raw_value is None:
        value = '.'
      else:
        value = str(raw_value)
      row.append(value)
    yield '\t'.join(row)


########## Contexts ##########

def read_context(context_file):
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

def format_context(context):
  output = []
  for i, base in enumerate(context.seq):
    if i == context.index:
      output.append('|'+base+'|')
    else:
      output.append(base)
  return ''.join(output)

def are_tuples_equal(tup1, tup2):
  for attr in dir(tup1):
    if attr.startswith('_') or attr == 'count':
      continue
    if getattr(tup1, attr) != getattr(tup2, attr):
      return False
  return True

Context = collections.namedtuple('Context', ('chr', 'pos', 'index', 'base', 'seq', 'gc'))
Context.__eq__ = are_tuples_equal
Context.__str__ = format_context


########## Errors ##########

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

def read_errors(errors_file):
  pair_line = None
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

Error = collections.namedtuple('Error', [meta['key'] for meta in ERROR_METADATA]+['pair'])
Pair = collections.namedtuple('Pair', [meta['key'] for meta in PAIR_METADATA])
Error.__eq__ = Pair.__eq__ = are_tuples_equal

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
