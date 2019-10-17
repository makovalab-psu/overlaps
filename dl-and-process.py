#!/usr/bin/env python3
import argparse
import configparser
import logging
import os
from pathlib import Path
from typing import Dict, List, Tuple, Any, Union, Optional, cast
import subprocess
import sys
import time
assert sys.version_info.major >= 3, 'Python 3 required'

SCRIPT_DIR = Path(__file__).resolve().parent
TMP_DIR = Path('~/tmp').expanduser()
DESCRIPTION = """Automatically download and analyze an SRA run."""


def make_argparser() -> argparse.ArgumentParser:
  parser = argparse.ArgumentParser(add_help=False, description=DESCRIPTION)
  io = parser.add_argument_group('I/O')
  io.add_argument('meta_ref', type=Path)
  io.add_argument('seqs_to_refs', type=Path)
  io.add_argument('acc', metavar='SRA accession',
    help='Accession number of the run to process.')
  io.add_argument('outdir', type=Path)
  io.add_argument('-r', '--refs-dir', type=Path, required=True)
  io.add_argument('-p', '--progress-file', type=Path,
    help='File to track how much of the pipeline we\'ve completed. If it exists, this will read '
      'the file and assume it records the last step completed by a previous run. It will also '
      'record how far we got in this run.')
  params = parser.add_argument_group('Parameters')
  params.add_argument('-q', '--mapq', type=int,
    help='Minimum MAPQ required when examining aligned reads.')
  params.add_argument('--min-ref-size', type=int,
    help='Minimum size of reference sequence to consider when choosing between them.')
  options = parser.add_argument_group('Options')
  options.add_argument('-b', '--begin', type=int,
    help='Start at this step instead of the beginning.')
  options.add_argument('-t', '--threads', type=int, default=1,
    help='Number of threads to use when aligning to the reference. Default: %(default)s')
  options.add_argument('-S', '--slurm', action='store_true',
    help='Run subcommands on Slurm cluster.')
  options.add_argument('--mem-ratio', type=int, default=500,
    help='Default: %(default)s bytes per base')
  options.add_argument('--min-mem', type=int, default=16*1024*1024*1024,
    help='Default: %(default)s bytes')
  options.add_argument('--max-mem', type=int, default=503*1024*1024*1024,
    help='Default: %(default)s bytes')
  options.add_argument('-h', '--help', action='help',
    help='Print this argument help text and exit.')
  logs = parser.add_argument_group('Logging')
  logs.add_argument('-l', '--log', type=argparse.FileType('w'), default=sys.stderr,
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  volume = logs.add_mutually_exclusive_group()
  volume.add_argument('-Q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
    default=logging.WARNING)
  volume.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  volume.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)
  return parser


def main(argv: List[str]) -> int:

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')

  # Process arguments.
  fq1_path, fq2_path = get_fq_paths(args.outdir)
  mem_params = {'min_mem':args.min_mem, 'max_mem':args.max_mem, 'ratio':args.mem_ratio}
  begin = 1
  if args.progress_file:
    if args.begin:
      begin = args.begin
    elif args.progress_file.is_file():
      progress = read_config_section(args.progress_file, 'end', {'step':int})
      begin = progress.get('step', 0)+1
    update_config(args.progress_file, 'start', step=begin, when=int(time.time()))
    del_config_section(args.progress_file, 'end')

  # Step 1: Download.
  if begin <= 1:
    download(args.acc, args.outdir, slurm=args.slurm)
    record_progress(args.progress_file, 1)

  # Step 2: Align.
  if begin <= 2:
    # Check the size of the data.
    fq_bases, fq_bytes = get_fq_size(fq1_path, fq2_path)
    if args.slurm:
      # Determine how much memory to request from Slurm.
      mem_req = format_mem_req(get_mem_req(fq_bases, fq_bytes, **mem_params))
    align(
      args.outdir, args.refs_dir, args.seqs_to_refs, args.meta_ref, args.mapq, args.min_ref_size,
      args.threads, args.slurm, mem_req, args.acc
    )
    record_progress(args.progress_file, 2)

  # Step 3: Find errors in overlaps.
  if begin <= 3:
    overlap(args.outdir/'align.auto.bam', args.outdir, args.mapq, args.slurm, args.acc)
    record_progress(args.progress_file, 3)

  # Step 4: Calculate statistics on errors.
  if begin <= 4:
    analyze(args.outdir/'errors.tsv', args.outdir, args.acc, args.slurm)
    record_progress(args.progress_file, 4)

  return 0


def get_fq_paths(outdir: Path) -> Tuple[Path,Path]:
  return outdir/'reads_1.fastq', outdir/'reads_2.fastq'


def record_progress(progress_path: Optional[Path], step: int) -> None:
  if progress_path is None:
    return
  update_config(progress_path, 'end', step=step, when=int(time.time()))


def update_config(config_path: Path, section: str, **kwargs) -> None:
  config = configparser.ConfigParser(interpolation=None)
  if config_path.is_file():
    config.read(config_path)
  data: Union[dict,configparser.SectionProxy]
  try:
    data = config[section]
  except KeyError:
    config.add_section(section)
    data = config[section]
  for key, value in kwargs.items():
    data[key] = str(value)
  with config_path.open('w') as config_file:
    config.write(config_file)


def del_config_section(config_path: Path, section: str) -> None:
  if not config_path.is_file():
    return
  config = configparser.ConfigParser(interpolation=None)
  config.read(config_path)
  try:
    del config[section]
  except KeyError:
    return
  with config_path.open('w') as config_file:
    config.write(config_file)


def read_config_section(
    config_path: Path, section: str, types: Optional[Dict[str,type]]=None
) -> Dict[str,Any]:
  data = {}
  config = configparser.ConfigParser(interpolation=None)
  try:
    config.read(config_path)
    for key, raw_value in config.items(section):
      if types and key in types:
        value = types[key](raw_value)
      else:
        value = raw_value
      data[key] = value
  except configparser.Error:
    fail(f'Invalid config file format in {config_path!r}.')
  return data


def download(acc: str, outdir: Path, slurm: bool=False) -> None:
  cmd: List = ['fasterq-dump', '--threads', '8', '--temp', TMP_DIR, acc, '-o', outdir/'reads']
  if slurm:
    cmd = ['srun', '-C', 'new', '-J', acc+':download', '--cpus-per-task', '8', '--mem', '10G'] + cmd
  run_command(cmd, onerror='fail', exe='fasterq-dump')


def align(
    outdir: Path, refs_dir: Path, seqs_to_refs: Path, meta_ref: Path, mapq: Optional[int],
    min_ref_size: Optional[int], threads: int, slurm: bool=False, mem_req: Optional[str]=None,
    acc: Optional[str]=None,
) -> None:
  fq1_path, fq2_path = get_fq_paths(outdir)
  cmd: List = [
    SCRIPT_DIR/'align-multi.py', '--clobber', '--threads', threads, '--name-sort', '--keep-tmp',
    '--refs-dir', refs_dir, seqs_to_refs, meta_ref, fq1_path, fq2_path,
    '--ref-counts', outdir/'ref-counts.tsv', '--output', outdir/'align.auto.bam'
  ]
  if mapq is not None:
    cmd[1:1] = ['--mapq', mapq]
  if min_ref_size is not None:
    cmd[1:1] = ['--min-size', min_ref_size]
  if slurm:
    acc = cast(str, acc)
    cmd = [
      'srun', '-C', 'new', '-J', acc+':align', '--cpus-per-task', threads, '--mem', mem_req
    ] + cmd
  run_command(cmd, onerror='fail', exe='align-multi.py')


def overlap(
    align_path: Path, outdir: Path, mapq: Optional[int], slurm: bool=False, acc: Optional[str]=None
) -> None:
  cmd: List = [
    SCRIPT_DIR/'overlaps.py', '--details', align_path, '--progress', '0',
    '--output2', 'summary', outdir/'errors.summary.tsv', '--output', outdir/'errors.tsv'
  ]
  if mapq is not None:
    cmd[1:1] = ['--mapq', mapq]
  if slurm:
    acc = cast(str, acc)
    cmd = ['srun', '-C', 'new', '-J', acc+':overlaps', '--mem', '24G'] + cmd
  run_command(cmd, onerror='fail', exe='overlaps.py')


def analyze(errors_path: Path, outdir: Path, acc: str, slurm: bool=False) -> None:
  cmd: List = [
    SCRIPT_DIR/'analyze.py', '--tsv', '--errors', acc, errors_path, '--output', outdir/'analysis.tsv'
  ]
  if slurm:
    cmd = ['srun', '-C', 'new', '-J', acc+':analyze', '--mem', '24G'] + cmd
  run_command(cmd, onerror='fail', exe='analyze.py')


def get_fq_size(fq1_path: Path, fq2_path: Path) -> Tuple[Optional[int],Optional[int]]:
  cmd = ('bioawk', '-c', 'fastx', '{tot+=length($seq)} END {print tot}', fq1_path, fq2_path)
  result = subprocess.run(cmd, stdout=subprocess.PIPE, encoding='utf8')
  if result.returncode != 0:
    logging.warning(f'Warning: Could not determine number of bases in sample.')
    bases = None
  else:
    bases = int(result.stdout.strip())
  bytes_ = os.path.getsize(fq1_path) + os.path.getsize(fq2_path)
  logging.info(f'Info: Sample is {bases} bases and {bytes_} bytes.')
  return bases, bytes_


def get_mem_req(
    bases: Optional[int], bytes_: Optional[int], ratio: int=400, min_mem: Optional[int]=None,
    max_mem: Optional[int]=None
) -> int:
  if bases is not None: 
    mem = ratio*bases
  elif bytes_ is not None:
    mem = int(ratio*bytes_/2.6)
  elif max_mem is not None:
    mem = max_mem
  else:
    mem = 128*1024*1024*1024
  if min_mem is not None and mem < min_mem:
    return cast(int, min_mem)
  elif max_mem is not None and mem > max_mem:
    return cast(int, max_mem)
  else:
    return mem


def format_mem_req(mem_bytes: int) -> str:
  mem_kb = round(mem_bytes/1024)
  return f'{mem_kb}K'


def run_command(cmd_raw: List, onerror: str='warn', exe: Optional[str]=None) -> int:
  cmd: List[str] = list(map(str, cmd_raw))
  logging.warning('+ $ '+' '.join(cmd))
  result = subprocess.run(cmd)
  if result.returncode != 0:
    if exe is None:
      exe = os.path.basename(cmd[0])
    message = f'{exe} failed with exit code {result.returncode}.'
    if onerror == 'warn':
      logging.warning('Warning: '+message)
    elif onerror == 'fail':
      fail(message)
  return result.returncode


def fail(message: str) -> None:
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
