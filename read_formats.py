"""Parse formats produced by the scripts in this repo."""
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


PAIR_METADATA = [                         # Column
  None,                                   # 1
  {'key':'name', 'type':str},             # 2
  {'key':'mapped', 'type':bool_literal},  # 3
  {'key':'rlen1', 'type':int},            # 4
  {'key':'rlen2', 'type':int},            # 5
  None,                                   # 6
  None,                                   # 7
  {'key':'forward', 'type':bool_literal}, # 8
  None,                                   # 9
  None,                                   # 10
]
OVERLAP_METADATA = [           # Column
  None,                        # 1
  None,                        # 2
  None,                        # 3
  None,                        # 4
  None,                        # 5
  {'key':'len', 'type':int},   # 6
  None,                        # 7
  None,                        # 8
  {'key':'start', 'type':int}, # 9
  {'key':'end', 'type':int},   # 10
]
ERROR_METADATA = [                   # Column
  None,                              # 1
  {'key':'type', 'type':str},        # 2
  {'key':'chr', 'type':str},         # 3
  {'key':'pos', 'type':int},         # 4
  {'key':'read_coord1', 'type':int}, # 5
  {'key':'read_coord2', 'type':int}, # 6
  {'key':'base1', 'type':str},       # 7
  {'key':'base2', 'type':str},       # 8
  {'key':'qual1', 'type':str},       # 9
  {'key':'qual2', 'type':str},       # 10
]

def get_fields(metadata):
  fields = []
  for meta in metadata:
    if meta is not None:
      fields.append(meta['key'])
  return fields

Pair = collections.namedtuple('Pair', get_fields(PAIR_METADATA)+['overlap','errors'])
Overlap = collections.namedtuple('Overlap', get_fields(OVERLAP_METADATA))
Error = collections.namedtuple('Error', get_fields(ERROR_METADATA))
Error.__eq__ = Pair.__eq__ = Overlap.__eq__ = tuple_eq


def read_errors_file(errors_file):
  """Parse the errors file (the output of overlap.py --details) and return the data.
  This will yield each Pair in the errors file. This includes the data in the "error" lines, since
  each Pair includes every associated Error as a list in its `errors` field."""
  pair = None
  for line_num, line_raw in enumerate(errors_file,1):
    value_strs = line_raw.rstrip('\r\n').split('\t')
    if value_strs[0] == 'pair':
      if pair is not None:
        yield pair
      if len(value_strs) != len(PAIR_METADATA):
        raise ValueError(
          f'Error on line {line_num}: Wrong number of columns (expected {len(PAIR_METADATA)}, '
          f'saw {len(value_strs)})'
        )
      overlap_values = parse_fields(value_strs, OVERLAP_METADATA)
      overlap = Overlap(**overlap_values)
      pair_values = parse_fields(value_strs, PAIR_METADATA)
      pair = Pair(overlap=overlap, errors=[], **pair_values)
    elif value_strs[0] == 'error':
      if len(value_strs) != len(ERROR_METADATA):
        raise ValueError(
          f'Error on line {line_num}: Wrong number of columns (expected {len(ERROR_METADATA)}, '
          f'saw {len(value_strs)})'
        )
      error_values = parse_fields(value_strs, ERROR_METADATA)
      error = Error(**error_values)
      pair.errors.append(error)
  if pair is not None:
    yield pair


def read_errors(errors_file):
  """Only read the errors from the errors file.
  For performance, this skips parsing "pair" lines with no associated errors.
  This still yields a series of Pairs, but only the Pairs with at least one Error."""
  pair_line = pair = None
  for line_raw in errors_file:
    if line_raw.startswith('pair\t'):
      if pair is not None:
        yield pair
      # Don't immediately parse this and create a new Pair.
      # There are often way more pair lines than errors, so that would usually be a waste of time.
      pair_line = line_raw
      # Invalidate any existing pair we've parsed, so that on the next error
      # line we'll make sure to parse and use this one.
      pair = None
    elif line_raw.startswith('error\t'):
      if pair is None:
        value_strs = pair_line.rstrip('\r\n').split('\t')
        overlap_values = parse_fields(value_strs, OVERLAP_METADATA)
        overlap = Overlap(**overlap_values)
        pair_values = parse_fields(value_strs, PAIR_METADATA)
        pair = Pair(overlap=overlap, errors=[], **pair_values)
      value_strs = line_raw.rstrip('\r\n').split('\t')
      error_values = parse_fields(value_strs, ERROR_METADATA)
      error = Error(**error_values)
      pair.errors.append(error)
    else:
      raise ValueError(f'Invalid first field in line {line_raw!r}')
  if pair is not None:
    yield pair


def parse_fields(raw_fields, metadata, null='.'):
  data = {}
  for raw_value, meta in zip(raw_fields, metadata):
    if meta is None:
      continue
    key = meta['key']
    if raw_value == null:
      data[key] = None
    else:
      data[key] = meta['type'](raw_value)
  return data
