#!/usr/bin/env python3
import argparse
import collections
import configparser
import logging
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple, Mapping, Sequence, Generator, Any, Union, Optional, cast
assert sys.version_info.major >= 3, 'Python 3 required'

Scalar = Union[str,int,float]
Progress = Mapping[str,Mapping[str,Scalar]]
PROGRESS_TYPES: Dict[str,type] = {
  'start_step':int, 'start_time':int, 'end_step':int, 'end_time':int, 'commit_time':int,
}
SCRIPT_DIR = Path(__file__).resolve().parent
TMP_DIR = Path('~/tmp').expanduser()
FQ_NAMES = [Path('reads_1.fastq'), Path('reads_2.fastq')]
DESCRIPTION = """Automatically download and analyze an SRA run."""

class FormatError(Exception):
  pass


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
  params.add_argument('-c', '--context', type=int, default=20,
    help='How much sequence context to extract around each error.')
  options = parser.add_argument_group('Options')
  options.add_argument('-b', '--begin', type=int,
    help='Start at this step instead of the beginning.')
  options.add_argument('-t', '--threads', type=int, default=1,
    help='Number of threads to use when aligning to the reference. Default: %(default)s')
  options.add_argument('-T', '--dl-threads', type=int, default=4,
    help='Number of threads to use when downloading from SRA. Default: %(default)s')
  options.add_argument('-S', '--slurm', action='store_true',
    help='Run subcommands on Slurm cluster.')
  options.add_argument('-w', '--wait-config', type=Path,
    help='The config file for slurm-wait.py, if you want to invoke it before launching each '
      'slurm job.')
  options.add_argument('-n', '--pick-node',
    help='Use slurm-wait.py to choose the node to run the job on, instead of letting slurm '
      'decide.')
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
  mem_params = {'min_mem':args.min_mem, 'max_mem':args.max_mem, 'ratio':args.mem_ratio}
  slurm_params: Optional[Dict[str,Any]] = None
  if args.slurm:
    slurm_params = {'pick_node':args.pick_node}
    if args.wait_config:
      slurm_params['config'] = args.wait_config
  # Figure out what step to begin on.
  begin = None
  if args.begin:
    begin = args.begin
    logging.info(f'Info: Starting at step {begin}, as directed by --begin.')
  elif args.progress_file and args.progress_file.is_file():
    logging.info('About to get begin from progress.')
    progress = read_progress(args.progress_file)
    last_end = get_last_step(progress)
    if isinstance(last_end, int):
      begin = last_end+1
      logging.info(
        f'Info: Starting at step {begin}, the step after the last recorded completed step.'
      )
  if not begin:
    begin = 1
    logging.info(f'Starting at step {begin} by default.')
  # Write starting data to progress file
  if args.progress_file:
    init_progress_file(args.progress_file, begin, SCRIPT_DIR)

  # Describe each step and how to perform it.
  steps: Sequence[Dict[str,Any]] = (
    {  # Step 1: Download.
      'fxn':download,
      'inputs':[],
      'args':[args.acc, args.outdir, args.dl_threads, slurm_params],
      'outputs':FQ_NAMES,
    },
    {  # Step 2: Align.
      'fxn':size_and_align,
      'inputs':[{'step':1, 'output':1, 'arg':1}, {'step':1, 'output':2, 'arg':2}],
      'args':[
        'dummy', 'dummy', mem_params, args.outdir, args.refs_dir,
        args.seqs_to_refs, args.meta_ref, args.mapq, args.min_ref_size, args.threads, args.acc,
        slurm_params
      ],
      'outputs':['align.auto.bam'],
    },
    {  # Step 3: Find errors in overlaps.
      'fxn':overlap,
      'inputs':[{'step':2, 'output':1, 'arg':1}],
      'args':['dummy', args.outdir, args.mapq, slurm_params, args.acc],
      'outputs':['errors.tsv'],
    },
    {  # Step 4: Extract the sequence context from the reference.
      'fxn':get_context,
      'inputs':[{'step':3, 'output':1, 'arg':1}, {'step':4, 'output':1, 'arg':4}],
      'args':['dummy', args.outdir/'metrics.tsv', args.refs_dir, 'dummy', args.context],
      'outputs':['seq-context.tsv']
    },
    {  # Step 5: Calculate statistics on errors.
      'fxn':analyze,
      'inputs':[{'step':3, 'output':1, 'arg':1}, {'step':5, 'output':1, 'arg':2}],
      'args':['dummy', 'dummy', args.acc, slurm_params],
      'outputs':['analysis.tsv'],
    },
  )

  run_steps(steps, begin, args.outdir, args.progress_file)

  return 0


########## MAIN HELPERS ##########

def run_steps(steps: Sequence[Dict[str,Any]], begin: int, outdir: Path, progress_file: Path):
  # Execute each step.
  for step_num, step in enumerate(steps, 1):
    if begin > step_num:
      continue
    if step_num > 1:
      fail_if_bad_output(outdir, *steps[step_num-2]['outputs'])
    args = insert_input_paths(steps, step['args'], step['inputs'], outdir)
    step['fxn'](*args)
    fail_if_bad_output(outdir, *step['outputs'])
    record_progress(progress_file, step_num)


def insert_input_paths(
    steps: Sequence[Dict[str,Any]], args_raw: List[Any], inputs: Sequence[Dict[str,int]],
    outdir: Path,
  ) -> List[Any]:
  args_list = args_raw.copy()
  for input_ in inputs:
    step_i = input_['step']-1
    output_i = input_['output']-1
    arg_i = input_['arg']-1
    output_name = steps[step_i]['outputs'][output_i]
    if args_list[arg_i] != 'dummy':
      fail(
        f'Problem with input file specification: argument {arg_i} of step {step_i} is '
        f"{args_list[arg_i]!r}, not 'dummy'. (Attempted to replace it with {output_name!r}.)"
      )
    args_list[arg_i] = outdir/output_name
  return args_list


def init_progress_file(progress_file: Path, start: int, script_dir: Path) -> None:
  progress = read_progress(progress_file)
  last_run = get_last_run(progress)
  if last_run is None:
    run = 1
  else:
    run = last_run+1
  section = {
    'start_step': start,
    'start_time': int(time.time())
  }
  git_info = get_git_info(script_dir)
  if git_info is None:
    logging.warning(f'Warning: Failed to find git commit info in {script_dir}')
  else:
    section['commit'] = git_info['commit']
    section['commit_time'] = git_info['timestamp']
  progress[f'run{run}'] = section
  write_config(progress_file, progress)


def record_progress(progress_path: Optional[Path], step: int) -> None:
  if progress_path is None:
    return
  progress = read_progress(progress_path)
  last_run = get_last_run(progress)
  new_data = {
    f'run{last_run}': {
      'step': step,
      'when': int(time.time())
    }
  }
  update_config(progress_path, new_data)


def get_last_step(progress):
  last_run = get_last_run(progress)
  if last_run is None:
    return None
  last_section = progress.get(f'run{last_run}')
  if last_section:
    return last_section.get('end_step')
  else:
    return None


def get_last_run(progress):
  last_run = None
  for section in progress.keys():
    if section.startswith('run'):
      run = int(section[3:])
      if last_run is None or run > last_run:
        last_run = run
  return last_run


def read_progress(progress_path: Path) -> Progress:
  raw_progress = read_config(progress_path, PROGRESS_TYPES)
  return convert_progress(raw_progress)


def convert_progress(progress: Progress) -> Progress:
  """Convert old progress structure to the new one, if necessary."""
  if any([section.startswith('run') for section in progress.keys()]):
    # It's the new format.
    return progress
  mapping = {
    ('start', 'step'): 'start_step',
    ('start', 'when'): 'start_time',
    ('end', 'step'): 'end_step',
    ('end', 'when'): 'end_time',
    ('version', 'timestamp'): 'commit_time',
    ('version', 'commit'): 'commit',
  }
  run0 = {}
  for section_name, section in progress.items():
    for key, value in section.items():
      new_key = mapping.get((section_name, key), f'{section_name}_{key}')
      run0[new_key] = value
  if run0:
    return {'run0':run0}
  else:
    return {}


def read_config(config_path: Path, types: Dict[str,type]=None) -> Progress:
  data: Progress = {}
  config = configparser.ConfigParser(interpolation=None)
  try:
    config.read(config_path)
    for section in config.sections():
      for key, raw_value in config.items(section):
        if types and key in types:
          value = types[key](raw_value)
        else:
          value = raw_value
        try:
          data[section][key] = value
        except KeyError:
          data[section] = {key:value}
  except configparser.Error:
    logging.critical(f'Error: Invalid config file format in {config_path!r}.')
    raise
  return data


def update_config(config_path: Path, new_data: Progress) -> None:
  config = configparser.ConfigParser(interpolation=None)
  if config_path.is_file():
    config.read(config_path)
  for section, section_data in new_data.items():
    data: Union[dict,configparser.SectionProxy]
    try:
      data = config[section]
    except KeyError:
      config.add_section(section)
      data = config[section]
    for key, value in section_data.items():
      data[key] = str(value)
  with config_path.open('w') as config_file:
    config.write(config_file)


def write_config(config_path: Path, data: Progress) -> None:
  config = configparser.ConfigParser(interpolation=None)
  for section_name, section_data in data.items():
    config.add_section(section_name)
    section = config[section_name]
    for key, value in section_data.items():
      section[key] = str(value)
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


def dict_get(d: Mapping, *keys: Any) -> Optional[Any]:
  current_value = d
  for key in keys:
    try:
      if key in current_value:
        current_value = current_value[key]
      else:
        return None
    except TypeError:
      return None
  return current_value


def get_git_info(git_dir: Path) -> Optional[Dict[str,Any]]:
  cmd = (
    'git', f'--work-tree={git_dir}', f'--git-dir={git_dir}/.git', 'log', '-n', '1', '--pretty=%ct %h'
  )
  try:
    result = subprocess.run(cmd, stdout=subprocess.PIPE, encoding='utf8')
  except OSError:
    return None
  if result.returncode != 0:
    return None
  fields = result.stdout.split()
  try:
    return {'timestamp':int(fields[0]), 'commit':fields[1]}
  except ValueError:
    return None


########## STEPS ##########

def download(acc: str, outdir: Path, threads: int=4, slurm_params: Dict[str,Any]=None) -> None:
  cmd: List = ['fasterq-dump', '--threads', threads, '--temp', TMP_DIR, acc, '-o', outdir/'reads']
  if slurm_params is not None:
    if 'config' in slurm_params:
      node = slurm_wait(config=slurm_params['config'], cpus=threads, mem='10G')
      specifier = get_slurm_specifier(node, slurm_params['pick_node'])
    cmd = [
      'srun', '-J', acc+':download', '--cpus-per-task', threads, '--mem', '10G'
    ] + specifier + cmd
  run_command(cmd, onerror='fail', exe='fasterq-dump')


def size_and_align(fq1_path: Path, fq2_path: Path, mem_params: Dict[str,int], outdir: Path,
    refs_dir: Path, seqs_to_refs: Path, meta_ref: Path, mapq: Optional[int],
    min_ref_size: Optional[int], threads: int, acc: str=None, slurm_params: Dict[str,Any]=None
  ) -> None:
  fq_bases, fq_bytes = get_fq_size(fq1_path, fq2_path)
  write_fq_size(outdir/'metrics.tsv', fq_bases, fq_bytes)
  mem_req = format_mem_req(get_mem_req(fq_bases, fq_bytes, **mem_params))
  align(
    outdir, refs_dir, seqs_to_refs, meta_ref, mapq, min_ref_size, threads, slurm_params, mem_req, acc
  )


def align(
    outdir: Path, refs_dir: Path, seqs_to_refs: Path, meta_ref: Path, mapq: Optional[int],
    min_ref_size: Optional[int], threads: int, slurm_params: Dict[str,Any]=None, mem_req: str=None,
    acc: str=None,
) -> None:
  fq1_path, fq2_path = [outdir/name for name in FQ_NAMES]
  cmd: List = [
    SCRIPT_DIR/'align-multi.py', '--clobber', '--threads', threads, '--name-sort', '--keep-tmp',
    '--refs-dir', refs_dir, seqs_to_refs, meta_ref, fq1_path, fq2_path,
    '--ref-counts', outdir/'ref-counts.tsv', '--run-info', outdir/'metrics.tsv',
    '--output', outdir/'align.auto.bam'
  ]
  if mapq is not None:
    cmd[1:1] = ['--mapq', mapq]
  if min_ref_size is not None:
    cmd[1:1] = ['--min-size', min_ref_size]
  if slurm_params is not None:
    if 'config' in slurm_params:
      node = slurm_wait(config=slurm_params['config'], cpus=threads, mem=mem_req)
      specifier = get_slurm_specifier(node, slurm_params['pick_node'])
    acc = cast(str, acc)
    cmd = [
      'srun', '-J', acc+':align', '--cpus-per-task', threads, '--mem', mem_req
    ] + specifier + cmd
  run_command(cmd, onerror='fail', exe='align-multi.py')


def overlap(
    align_path: Path, outdir: Path, mapq: Optional[int], slurm_params: Dict[str,Any]=None,
    acc: str=None
) -> None:
  cmd: List = [
    SCRIPT_DIR/'overlaps.py', '--details', align_path, '--progress', '0',
    '--output2', 'summary', outdir/'errors.summary.tsv', '--output', outdir/'errors.tsv'
  ]
  if mapq is not None:
    cmd[1:1] = ['--mapq', mapq]
  if slurm_params is not None:
    if 'config' in slurm_params:
      node = slurm_wait(config=slurm_params['config'], cpus=1, mem='24G')
      specifier = get_slurm_specifier(node, slurm_params['pick_node'])
    acc = cast(str, acc)
    cmd = ['srun', '-J', acc+':overlaps', '--mem', '24G'] + specifier + cmd
  run_command(cmd, onerror='fail', exe='overlaps.py')


def get_context(
    errors_path: Path, metrics_path: Path, refs_dir: Path, out_path: Path, con_size: int
  ):
  metrics = read_metrics(metrics_path)
  ref_path = refs_dir/metrics['ref']
  cmd_raw: List = [
    'getcontext.py', ref_path, '--chrom-field', 1, '--coord-field', 2, '--window', con_size,
    '--output', out_path
  ]
  cmd = list(map(str, cmd_raw))
  logging.warning('+ $ '+' '.join(cmd))
  process = subprocess.Popen(cmd, stdin=subprocess.PIPE, encoding='utf8')
  for ref, coord in extract_sites(errors_path):
    process.stdin.write(f'{ref}\t{coord}\n')
  process.stdin.close()
  exit_code = process.wait()
  if exit_code!= 0:
    fail(f'getcontext.py failed with exit code {exit_code}.')


def analyze(
    errors_path: Path, analysis_path: Path, acc: str, slurm_params: Dict[str,Any]=None
  ) -> None:
  cmd: List = [
    SCRIPT_DIR/'analyze.py', '--tsv', '--errors', acc, errors_path, '--output', analysis_path
  ]
  if slurm_params is not None:
    if 'config' in slurm_params:
      node = slurm_wait(config=slurm_params['config'], cpus=1, mem='24G')
      specifier = get_slurm_specifier(node, slurm_params['pick_node'])
    cmd = ['srun', '-J', acc+':analyze', '--mem', '24G'] + specifier + cmd
  run_command(cmd, onerror='fail', exe='analyze.py')


########## STEP HELPERS ##########

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


def write_fq_size(out_path: Path, bases: Optional[int], bytes_: Optional[int]) -> None:
  if bases is None:
    bases_str = '.'
  else:
    bases_str = str(bases)
  if bytes_ is None:
    bytes_str = '.'
  else:
    bytes_str = str(bytes_)
  with out_path.open('w') as out_file:
    print('bases', bases_str, sep='\t', file=out_file)
    print('bytes', bytes_str, sep='\t', file=out_file)


def get_mem_req(
    bases: Optional[int], bytes_: Optional[int], ratio: int=400, min_mem: int=None,
    max_mem: int=None
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


def extract_sites(errors_path: Path) -> Generator[Tuple[str,int],None,None]:
  with errors_path.open() as errors_file:
    for line_raw in errors_file:
      line = line_raw.rstrip('\r\n')
      if not line:
        continue
      fields = line.split('\t')
      if fields[0] != 'error':
        continue
      if len(fields) != 8:
        raise FormatError(
          f'Errors file {str(errors_path)!r} has wrong number of fields on line: {line_raw!r}.'
        )
      ref = fields[2]
      try:
        coord = int(fields[3])
      except ValueError:
        raise FormatError(
          f'Errors file {str(errors_path)!r} has invalid coordinate on line: {line_raw!r}'
        )
      yield ref, coord


def read_metrics(metrics_path: Path) -> Dict[str,Any]:
  metrics = {}
  with metrics_path.open() as metrics_file:
    for line_raw in metrics_file:
      line = line_raw.rstrip('\r\n')
      if not line:
        continue
      fields = line.split('\t')
      value: Union[Scalar,List[Scalar]]
      if len(fields) < 2:
        raise FormatError(
          f'Invalid format in {str(metrics_path)!r}: Too few fields in line {line_raw!r}'
        )
      elif len(fields) == 2:
        key = fields[0]
        value = parse_value(fields[1])
      elif len(fields) > 2:
        key = fields[0]
        value = [parse_value(raw) for raw in fields[1:]]
      metrics[key] = value
  return metrics


def parse_value(raw_value: str) -> Scalar:
  value: Union[str,int,float]
  try:
    value = int(raw_value)
  except ValueError:
    try:
      value = float(raw_value)
    except ValueError:
      value = raw_value
  return value


def slurm_wait(config: Path=None, cpus: int=None, mem: str=None) -> Optional[str]:
  cmd_raw: List[Any] = ['slurm-wait.py']
  if config:
    cmd_raw.extend(['--config', config])
  if cpus:
    cmd_raw.extend(['--cpus', cpus])
  if mem:
    cmd_raw.extend(['--mem', mem])
  cmd = list(map(str, cmd_raw))
  result = subprocess.run(cmd, stdout=subprocess.PIPE, encoding='utf8')
  if result.returncode != 0:
    fail(f'slurm-wait.py failed with exit code {result.returncode}.')
  node = result.stdout.strip()
  if node:
    return node
  else:
    return None


def get_slurm_specifier(node: Optional[str], pick_node: bool=False):
  if pick_node and node is not None:
    return ['-w', node]
  else:
    return ['-C', 'new']


def run_command(cmd_raw: List, onerror: str='warn', exe: str=None) -> int:
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


def fail_if_bad_output(outdir: Path, *filenames: str) -> None:
  for filename in filenames:
    path = outdir/filename
    if not path.is_file():
      fail(f'Error: File {str(path)!r} not found.')
    if os.path.getsize(path) == 0:
      fail(f'Error: File {str(path)!r} is empty.')


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
