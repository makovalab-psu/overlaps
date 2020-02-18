#!/usr/bin/env python3
import argparse
import logging
import os
import requests
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Optional, Sequence, List
assert sys.version_info.major >= 3, 'Python 3 required'

EBI_SRA_URL = 'https://www.ebi.ac.uk/ena/data/warehouse/filereport?result=read_run&fields=fastq_ftp&accession='
DESCRIPTION = """Download a sample from the SRA."""


def make_argparser():
  parser = argparse.ArgumentParser(add_help=False, description=DESCRIPTION)
  options = parser.add_argument_group('Options')
  options.add_argument('accession', metavar='dispname',
    help='Sample SRA accession number.')
  options.add_argument('outdir', default=Path('.'), type=Path, nargs='?',
    help='Output directory.')
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

  download(args.accession, args.outdir)


def download(accession: str, outdir: Path) -> None:
  for i, url in enumerate(get_ftp_urls(accession), 1):
    cmd: List = ['curl', '--silent', '--show-error', url]
    proc1 = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    with (outdir/f'reads_{i}.fastq').open('w') as fq_file:
      proc2 = subprocess.Popen(['gunzip', '-c', '-'], stdin=proc1.stdout, stdout=fq_file)
      run_pipeline((proc1, proc2), onerror='fail', exes=('curl', 'gunzip'))


def get_ftp_urls(accession: str) -> List[str]:
  ftp_urls: List[str] = []
  url = EBI_SRA_URL+accession
  response = make_request(url)
  for line_num, line in enumerate(response.text.splitlines(), 1):
    if line_num == 1:
      if line != 'fastq_ftp':
        raise FormatError(f'Invalid response from {url!r}: wrong first line ({line!r})')
    elif line_num == 2:
      fields = line.split(';')
      if len(fields) == 1 and fields[0].endswith('.fastq.gz'):
        raise FormatError(f'Invalid response from {url!r}: Looks single-ended? ({fields[0]!r})')
      elif len(fields) == 2 and all([u.startswith('ftp.sra.ebi.ac.uk') for u in fields]):
        ftp_urls = ['ftp://'+u for u in fields]
      else:
        raise FormatError(f'Invalid response from {url!r}: bad second line ({line!r})')
    else:
      if line_num == 3:
        logging.warning(f'Extra lines in response from {url!r}:')
      if 3 <= line_num < 10:
        logging.warning(f'  {line}')
      elif line_num == 10:
        logging.warning(f'  ...')
  if line_num < 2:
    raise FormatError(
      f'Invalid response from {url!r}: Not enough lines! Saw {line_num} instead of 2.'
    )
  return ftp_urls


def make_request(url: str, retry_wait=5, max_tries=5) -> requests.models.Response:
  response = None
  tries = 0
  while response is None:
    tries += 1
    try:
      response = requests.get(url)
    except requests.exceptions.RequestException as error:
      logging.error(f'Error: Issue with HTTP request to {url!r}: {error}')
    if response is not None and response.status_code != 200:
      logging.error(f'Error: HTTP {response.status_code} returned from {url!r}')
      response = None
    if response is None:
      if tries >= max_tries:
        raise RuntimeError(f'Failed to request {url!r}')
      time.sleep(retry_wait)
      retry_wait *= 4
  return response


def run_pipeline(procs: Sequence, onerror='warn', exes: Sequence[Optional[str]]=None):
  if exes is None:
    exes = [None]*len(procs)
  elif len(exes) != len(procs):
    raise RuntimeError(f'exes list must be as long as procs list. Received {exes!r}')
  procs[0].stdout.close()
  for proc, exe in zip(procs, exes):
    return_code = proc.wait()
    warn_or_fail(return_code, onerror, proc.args[0], exe)
  return return_code


def warn_or_fail(return_code: int, onerror: str, arg0: str, exe: str=None):
  if return_code == 0:
    return
  if exe is None:
    exe = os.path.basename(arg0)
  message = f'{exe} failed with exit code {return_code}.'
  if onerror == 'warn':
    logging.warning(f'Warning: {message}')
  elif onerror == 'fail':
    fail(message)


def fail(message):
  logging.critical('Error: '+str(message))
  if __name__ == '__main__':
    sys.exit(1)
  else:
    raise Exception(message)


class FormatError(Exception):
  pass


if __name__ == '__main__':
  try:
    sys.exit(main(sys.argv))
  except BrokenPipeError:
    pass
