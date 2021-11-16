#!/usr/bin/env python3
import argparse
import difflib
import enum
import logging
import pathlib
import subprocess
import sys
import types
assert sys.version_info.major >= 3, 'Python 3 required'

TESTS_DIR = pathlib.Path(__file__).resolve().parent
ROOT_DIR = TESTS_DIR.parent
DESCRIPTION = """"""


def make_argparser():
  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.add_argument('tests', metavar='test_name', nargs='*',
    help='The tests to run.')
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

  # Create the dicts holding all defined tests.

  meta_tests = get_objects_diff(GlobalsInitial, GlobalsAfterMeta)
  meta_tests['all'] = lambda name: run_test_group(simple_tests)
  meta_tests['active'] = lambda name: run_test_group(active_tests)
  meta_tests['inactive'] = lambda name: run_test_group(inactive_tests)

  active_tests = get_objects_diff(GlobalsAfterMeta, GlobalsAfterActive)
  inactive_tests = get_objects_diff(GlobalsAfterActive, GlobalsAfterInactive)
  simple_tests = add_dicts(active_tests, inactive_tests)

  all_tests = add_dicts(meta_tests, simple_tests)

  if not args.tests:
    print('Meta tests:')
    for test_name in meta_tests.keys():
      print('  '+test_name)
    print('Active tests:')
    for test_name in active_tests.keys():
      print('  '+test_name)
    print('Inactive tests:')
    for test_name in inactive_tests.keys():
      print('  '+test_name)
    return 1

  unknown_tests = []
  for test_name in args.tests:
    if test_name not in all_tests:
      unknown_tests.append(test_name)
  if unknown_tests:
    fail('Error: Test(s) "{}" unrecognized.'.format('", "'.join(unknown_tests)))

  for test_name in args.tests:
    test_fxn = all_tests[test_name]
    test_fxn(test_name)

GlobalsInitial = globals().copy()


##### Meta tests #####

GlobalsAfterMeta = globals().copy()


##### Active tests #####

def overlapper(test_name):
  do_simple_test(
    test_name, 'overlaps.py', 'overlap.align.bam', 'overlap.details.out.tsv',
    ('--details', Placeholders.INPUT)
  )


def read_errors(test_name):
  do_simple_test(
    test_name, 'read_formats.py', 'overlap.details.out.tsv', 'read-errors.out.tsv',
    ('read_errors', Placeholders.INPUT), script=ROOT_DIR/'tests/read-formats-helper.py'
  )


def read_errors_file(test_name):
  do_simple_test(
    test_name, 'read_formats.py', 'overlap.details.out.tsv', 'read-errors-file.out.tsv',
    ('read_errors_file', Placeholders.INPUT), script=ROOT_DIR/'tests/read-formats-helper.py'
  )


GlobalsAfterActive = globals().copy()


##### Inactive tests #####

GlobalsAfterInactive = globals().copy()


##### Helper functions #####

class Placeholders(enum.Enum):
  INPUT = 1
  OUTPUT = 2


def do_simple_test(test_name, script_name, input_name, output_name, cmd_args, script=None):
  cmd_args = list(cmd_args)
  print(f'{test_name} ::: {script_name} ::: {input_name}\t', end='')
  if script is None:
    script = ROOT_DIR / script_name
  for i in range(len(cmd_args)):
    if cmd_args[i] == Placeholders.INPUT:
      cmd_args[i] = TESTS_DIR/input_name
    elif cmd_args[i] == Placeholders.OUTPUT:
      cmd_args[i] = TESTS_DIR/output_name
  result, exit_code = run_command_and_capture([script]+cmd_args, onerror='stderr')
  if exit_code != 0:
    print('FAILED')
  else:
    expected = read_file(TESTS_DIR/output_name)
    if result != expected:
      print('FAILED')
      for line in trimmed_diff(expected.splitlines(), result.splitlines()):
        print(line)
    else:
      print('success')


def add_dicts(*dicts):
  combined = {}
  for d in dicts:
    for key, value in d.items():
      combined[key] = value
  return combined


def run_test_group(test_group):
  for test_name, test_fxn in test_group.items():
    test_fxn(test_name)


def get_objects_diff(objects_before, objects_after, object_type=types.FunctionType):
  diff = {}
  for name, obj in objects_after.items():
    if name not in objects_before and isinstance(obj, object_type):
      diff[name] = obj
  return diff


def head(string, lines=10):
  string_lines = string.splitlines()
  output = '\n'.join(string_lines[:lines])
  if len(string_lines) > lines:
    output += '\n\t...'
  return output


def run_command(command, onerror='warn'):
  if onerror == 'stderr':
    result = run_command_and_catch(command, onerror=onerror, stderr=subprocess.PIPE)
    if result.returncode != 0:
      logging.error(str(result.stderr, 'utf8'))
  else:
    result = run_command_and_catch(command, onerror=onerror, stderr=subprocess.DEVNULL)
  return result.returncode


def run_command_and_capture(command, onerror='warn'):
  result = run_command_and_catch(
    command, onerror=onerror, stdout=subprocess.PIPE, stderr=subprocess.PIPE
  )
  if result.returncode != 0 and onerror == 'stderr':
    logging.error(str(result.stderr, 'utf8'))
  return str(result.stdout, 'utf8'), result.returncode


def run_command_and_catch(command, onerror='warn', **kwargs):
  try:
    result = subprocess.run(command, **kwargs)
  except (OSError, subprocess.CalledProcessError) as error:
    if onerror == 'stderr':
      pass
    elif onerror == 'warn':
      logging.error(f'Error: ({type(error).__name__}) {error}')
    elif onerror == 'raise':
      raise
  return result


def read_file(path):
  try:
    with path.open('r') as file:
      return file.read()
  except OSError as error:
    logging.error('Error: {}'.format(error))
    return None


def trimmed_diff(lines1, lines2, lineterm=''):
  """Get a trimmed diff.
  Input lines should be newline-free."""
  diff_lines = difflib.unified_diff(
    lines1, lines2, n=1, fromfile='a', tofile='b', fromfiledate='c', tofiledate='d',
    lineterm=lineterm
  )
  header_line = 0
  for line in diff_lines:
    if header_line == 0 and line == '--- a\tc'+lineterm:
      header_line = 1
    elif header_line == 1 and line == '+++ b\td'+lineterm:
      header_line = 2
    elif header_line == 2:
      header_line = None
    if header_line is None:
      yield line


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
