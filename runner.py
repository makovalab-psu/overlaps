#!/usr/bin/env python3
import argparse
import logging
import os
from pathlib import Path
import subprocess
import sys
import time
from typing import Union, Optional, cast, List, Tuple, Set
assert sys.version_info.major >= 3, 'Python 3 required'

SCRIPT_DIR = Path(__file__).resolve().parent
JOB_ARGS = ['--verbose', '--min-ref-size', 2000000, '--mem-ratio', 400, '--slurm']
DESCRIPTION = """"""


def make_argparser() -> argparse.ArgumentParser:
  parser = argparse.ArgumentParser(add_help=False, description=DESCRIPTION)
  io = parser.add_argument_group('I/O')
  io.add_argument('todo', metavar='todo.txt', type=Path,
    help='A file full of accession numbers we want to process.')
  io.add_argument('done', metavar='done.txt', type=Path,
    help='A file to track the accession numbers we\'ve finished.')
  io.add_argument('launched', metavar='launched.txt', type=Path,
    help='A file to track the accession numbers of jobs we\'ve launched.')
  io.add_argument('parent_dir', metavar='path/to/parent/dir', type=Path,
    help='The top of the project directory structure. This should be the direct parent of "runs" '
      'and "refs".')
  io.add_argument('this_config', metavar='slurm-wait.ini', type=Path,
    help='Control the pace of this script with slurm-wait.py.')
  io.add_argument('job_config', metavar='slurm-wait-job.ini', type=Path,
    help='Control the pace of the subprocesses with slurm-wait.py.')
  options = parser.add_argument_group('Options')
  options.add_argument('-b', '--begin', type=int,
    help='Start jobs at this step instead of the beginning.')
  options.add_argument('-t', '--threads', type=int, default=32,
    help='Default: %(default)s')
  options.add_argument('-T', '--dl-threads', type=int, default=4,
    help='Number of threads to use when downloading from SRA. Default: %(default)s')
  options.add_argument('-n', '--pick-node', action='store_true',
    help='Let slurm-wait.py specify which node to run each job on. Passes this option to '
      'dl-and-process.py.')
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


def main(argv: List[str]) -> Optional[int]:

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')

  last_acc = None
  while True:

    # Read in the accessions from their files to understand what's left to do.
    todo_in_file = read_file_as_list(args.todo)
    done = read_file_as_list(args.done)
    launched = read_file_as_list(args.launched)
    todo = subtract_lists(todo_in_file, done, launched)
    if not todo:
      print('Done! No more runs in the todo list.')
      return 0
    next_acc = todo[0]

    node = wait_for_node(args.this_config, args.threads, last_acc=last_acc)
    if node is False:
      # slurm-wait.py failed. Wait and try again.
      time.sleep(5)
      continue
    # Just as a safety to mitigate the effects of a runaway loop.
    time.sleep(1)

    run_dir = args.parent_dir/'runs'/next_acc
    try:
      os.mkdir(run_dir)
    except FileExistsError:
      pass
    print(f'Launching {next_acc}')
    launch_job(
      next_acc, args.parent_dir, args.job_config, args.threads, args.dl_threads, args.pick_node,
      args.begin,
    )
    launched.append(next_acc)
    write_list_to_file(args.launched, launched)
    last_acc = next_acc
    #TODO: Add launched jobs to "done" if their progress.ini says they're finished?
    #TODO: Relaunch failed jobs?
    #      Read progress.ini to find out how far they got.

  return 0


def read_file_as_list(list_path: Path) -> List[str]:
  data = []
  with list_path.open() as list_file:
    for line in list_file:
      data.append(line.strip())
  return data


def write_list_to_file(list_path: Path, data: List[str]) -> None:
  with list_path.open('w') as list_file:
    for item in data:
      print(item, file=list_file)


def subtract_lists(list1: list, *lists) -> list:
  """Return a copy of `list1` but remove any items present in the other lists."""
  new_list = []
  remove_these: Set[list] = set()
  for list2 in lists:
    remove_these |= set(list2)
  for item in list1:
    if item not in remove_these:
      new_list.append(item)
  return new_list


def wait_for_node(config_path: Path, threads: int, last_acc: str=None) -> Union[str,bool,None]:
  cmd_raw: list = ['slurm-wait.py', '--config', config_path, '--cpus', threads]
  if last_acc:
    #TODO: This could be caught in an infinite loop if the download for the last job finishes too
    #      quickly or this script gets held up (perhaps because of one failed slurm-wait.py command).
    #      Make sure to check whether the job already finished.
    #TODO: Also, do I want to allow for re-running incomplete jobs? Then the download step will be
    #      skipped entirely.
    cmd_raw += ['--wait-for-job', last_acc+':download']
  cmd = list(map(str, cmd_raw))
  result = subprocess.run(cmd, stdout=subprocess.PIPE, encoding='utf8')
  if result.returncode != 0:
    logging.warning(f'Warning: slurm-wait.py failed with exit code {result.returncode}.')
    return False
  node = result.stdout.strip()
  if node:
    return node
  else:
    return None


def launch_job(
    acc: str, parent_dir: Path, config: Path, threads: int, dl_threads: int, pick_node: bool,
    begin: Optional[int],
  ) -> None:
  run_dir = parent_dir/'runs'/acc
  cmd_raw = cast(list, [SCRIPT_DIR/'dl-and-process.py']) + JOB_ARGS + [
    '--threads', threads, '--dl-threads', dl_threads, '--wait-config', config, '--progress-file',
    run_dir/'progress.ini', '--refs-dir', parent_dir/'refs/individual',
    parent_dir/'refs/all.complete.fa', parent_dir/'refs/seqs_to_refs.tsv', acc, run_dir
  ]
  if pick_node:
    cmd_raw.insert(1, '--pick-node')
  if begin:
    cmd_raw[1:1] = ['--begin', begin]
  cmd = list(map(str, cmd_raw))
  print('$ '+' '.join(cmd))
  out_log = run_dir/'dl-and-process.out.log'
  err_log = run_dir/'dl-and-process.err.log'
  with out_log.open('w') as out_file, err_log.open('w') as err_file:
    subprocess.Popen(cmd, stdout=out_file, stderr=err_file)


def fail(message: str):
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
