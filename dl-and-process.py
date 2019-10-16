#!/usr/bin/env python3
import argparse
import logging
import os
from pathlib import Path
from typing import List, Tuple, Any, Optional, cast
import subprocess
import sys
assert sys.version_info.major >= 3, 'Python 3 required'

SCRIPT_DIR = Path(__file__).resolve().parent
TMP_DIR = Path('~/tmp').expanduser()
DESCRIPTION = """"""


def make_argparser() -> argparse.ArgumentParser:
  parser = argparse.ArgumentParser(add_help=False, description=DESCRIPTION)
  io = parser.add_argument_group('I/O')
  io.add_argument('acc', metavar='SRA accession',
    help='Accession number of the run to process.')
  io.add_argument('outdir', type=Path)
  io.add_argument('meta_ref', type=Path)
  io.add_argument('seqs_to_refs', type=Path)
  io.add_argument('-r', '--refs-dir', type=Path, required=True)
  params = parser.add_argument_group('Parameters')
  params.add_argument('-q', '--mapq', type=int,
    help='Minimum MAPQ required when examining aligned reads.')
  params.add_argument('--min-ref-size', type=int,
    help='Minimum size of reference sequence to consider when choosing between them.')
  options = parser.add_argument_group('Options')
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

  # Download.
  download(args.acc, args.outdir, slurm=args.slurm)

  # Check the size of the data.
  if args.slurm:
    fq_bases, fq_bytes = get_fq_size(fq1_path, fq2_path)
    mem_req = format_mem_req(get_mem_req(fq_bases, fq_bytes, **mem_params))

  # Align.
  align(
    args.outdir, args.refs_dir, args.seqs_to_refs, args.meta_ref, args.mapq, args.min_ref_size,
    args.threads, args.slurm, mem_req, args.acc
  )

  # Find errors in overlaps.
  overlap(args.outdir/'align.auto.bam', args.outdir, args.mapq, args.slurm, args.acc)

  # Calculate statistics on errors.
  analyze(args.outdir/'errors.tsv', args.outdir, args.acc, args.slurm)

  return 0


def get_fq_paths(outdir: Path) -> Tuple[Path,Path]:
  return outdir/'reads_1.fastq', outdir/'reads_2.fastq'


def download(acc: str, outdir: Path, slurm: bool=False) -> None:
  cmd: List = ['fasterq-dump', '--threads', '8', '--temp', TMP_DIR, acc, '-o', outdir/'reads']
  if slurm:
    cmd = ['srun', '-C', 'new', '-J', acc+':down', '--cpus-per-task', '8', '--mem', '10G'] + cmd
  exit_code = run_command(cmd)
  if exit_code != 0:
    fail(f'fasterq-dump failed with exit code {exit_code}.')


def align(
    outdir: Path, refs_dir: Path, seqs_to_refs: Path, meta_ref: Path, mapq: Optional[int],
    min_ref_size: Optional[int], threads: int, slurm: bool=False, mem_req: Optional[str]=None,
    acc: Optional[str]=None,
) -> None:
  fq1_path, fq2_path = get_fq_paths(outdir)
  cmd: List = [
    SCRIPT_DIR/'align-multi.py', '--threads', threads, '--name-sort', '--keep-tmp',
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
  exit_code = run_command(cmd)
  if exit_code != 0:
    fail(f'align-multi.py failed with exit code {exit_code}.')


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
    cmd = ['srun', '-C', 'new', '-J', acc+':over', '--mem', '24G'] + cmd
  exit_code = run_command(cmd)
  if exit_code != 0:
    fail(f'overlaps.py failed with exit code {exit_code}.')


def analyze(errors_path: Path, outdir: Path, acc: str, slurm: bool=False) -> None:
  cmd: List = [
    SCRIPT_DIR/'analyze.py', '--tsv', '--errors', acc, errors_path, '--output', outdir/'analysis.tsv'
  ]
  if slurm:
    cmd = ['srun', '-C', 'new', '-J', acc+':analyze', '--mem', '24G'] + cmd
  exit_code = run_command(cmd)
  if exit_code != 0:
    fail(f'analyze.py failed with exit code {exit_code}.')


def get_fq_size(fq1_path: Path, fq2_path: Path) -> Tuple[Optional[int],Optional[int]]:
  cmd = ('bioawk', '-c', 'fastx', '{tot+=length($seq)} END {print tot}', fq1_path, fq2_path)
  result = subprocess.run(cmd, stdout=subprocess.PIPE, encoding='utf8')
  if result.returncode != 0:
    return None, None
  bases = int(result.stdout.strip())
  bytes_ = os.path.getsize(fq1_path) + os.path.getsize(fq2_path)
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


def run_command(cmd_raw: List) -> None:
  cmd: List[str] = list(map(str, cmd_raw))
  logging.warning('+ $ '+' '.join(cmd))
  result = subprocess.run(cmd)
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
