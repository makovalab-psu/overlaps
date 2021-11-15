"""Parse errors.tsv (the output of overlap.py --details)."""
import collections


def tuple_eq(tup1, tup2):
  for attr in dir(tup1):
    if attr.startswith('_') or attr == 'count':
      continue
    if getattr(tup1, attr) != getattr(tup2, attr):
      return False
  return True


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
  {'key':'qual1', 'type':str},
  {'key':'qual2', 'type':str},
]
PAIR_METADATA = [
  {'key':'name', 'type':str},
  {'key':'mapped', 'type':bool_literal},
  {'key':'rlen1', 'type':int},
  {'key':'rlen2', 'type':int},
  {'key':'overlap', 'type':int},
  {'key':'errors', 'type':int},
  {'key':'forward', 'type':bool_literal},
  {'key':'ostart', 'type':int},
  {'key':'oend', 'type':int},
]


Error = collections.namedtuple('Error', [meta['key'] for meta in ERROR_METADATA]+['pair'])
Pair = collections.namedtuple('Pair', [meta['key'] for meta in PAIR_METADATA])
Error.__eq__ = Pair.__eq__ = tuple_eq


def read_errors(errors_file):
  pair_line = pair = None
  for line_raw in errors_file:
    if line_raw.startswith('pair\t'):
      # Don't immediately parse this and create a new Pair.
      # There are often way more pair lines than errors, so that would usually be a waste of time.
      pair_line = line_raw
      # Invalidate any existing pair we've parsed, so that on the next error
      # line we'll make sure to parse and use this one.
      pair = None
    elif line_raw.startswith('error\t'):
      if pair is None:
        pair_fields = pair_line.rstrip('\r\n').split('\t')
        pair = parse_fields(pair_fields[1:], Pair, PAIR_METADATA)
      error_fields = line_raw.rstrip('\r\n').split('\t')
      error = parse_fields(error_fields[1:], Error, ERROR_METADATA, pair=pair)
      yield error
    else:
      raise ValueError(f'Invalid first field in line {line_raw!r}')


def parse_fields(raw_fields, tuple_class, metadata, null='.', **kwargs):
  data = {}
  for raw_value, meta in zip(raw_fields, metadata):
    key = meta['key']
    if raw_value == null:
      data[key] = None
    else:
      data[key] = meta['type'](raw_value)
  for key, value in kwargs.items():
    data[key] = value
  return tuple_class(**data)
