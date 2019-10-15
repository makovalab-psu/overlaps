#!/usr/bin/env python3
import argparse
import collections
import functools
import logging
import operator
import os
from pathlib import Path
import random
import subprocess
import sys
from typing import List, Dict, Tuple
from bfx import samreader
assert sys.version_info.major >= 3, 'Python 3 required'

BAD_FLAGS = (4, 256, 512, 1024, 2048)
BAD_FLAG_INT = functools.reduce(operator.or_, BAD_FLAGS)
DESCRIPTION = """Align a sample to the correct reference sequence.
This is for when it's unclear which reference sequence is appropriate for a sample.
It will determine this by aligning the reads to a "meta-genome" of all candidate reference sequences
concatenated, then choose the reference where the most reads align.
Then it will do another alignment against that reference only."""


def make_argparser() -> argparse.ArgumentParser:
  parser = argparse.ArgumentParser(add_help=False, description=DESCRIPTION)
  options = parser.add_argument_group('Options')
  options.add_argument('meta_ref', type=Path,
    help='Meta-reference file containing all possible sequences.')
  options.add_argument('fastq1', type=Path,
    help='Input reads, mate 1')
  options.add_argument('fastq2', type=Path,
    help='Input reads, mate 2')
  options.add_argument('seqs_to_refs', type=Path,
    help='File mapping sequence names to reference files. Should contain one line per sequence, '
      'with two tab-delimited columns: the sequence name, and the path to the reference file. '
      'If --refs is given, the path will be interpreted as a relative path under the --refs '
      'directory.')
  options.add_argument('-o', '--output', type=Path,
    help='Write the output BAM to this file. Default: choose a filename based on the first fastq '
      'file.')
  options.add_argument('--refs-dir', type=Path,
    help='Directory containing all reference files.')
  options.add_argument('-m', '--mapq', type=int,
    help='MAPQ threshold when evaluating meta-alignment.')
  options.add_argument('-s', '--min-size', type=int,
    help="Don't consider sequences smaller than this many bases when finding the sequence with the "
      'most alignments.')
  options.add_argument('-k', '--keep-tmp', action='store_true',
    help='Keep intermediate temporary files.')
  options.add_argument('-t', '--threads', type=int, default=1)
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


def main(argv: List[str]) -> int:

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')

  # Process arguments and create paths.
  base = get_reads_base(args.fastq1)
  align_tmp_path = get_tmp_path(base+'.align.meta', 'bam')
  if args.output:
    out_path = args.output
  else:
    out_path = get_tmp_path(base+'.align', 'bam')

  # Read sequence-to-reference mapping.
  seqs_to_refs, seq_sizes = read_map(args.seqs_to_refs)

  # Align to meta-reference.
  align(args.meta_ref, args.fastq1, args.fastq2, align_tmp_path, threads=args.threads)

  # Count alignments to each sequence.
  logging.warning('Counting alignments to each reference sequence..')
  aln_counts = count_alignments(align_tmp_path, args.mapq)
  best_seq = get_best_seq(aln_counts, seq_sizes, args.min_size)
  if best_seq is None:
    fail(f'No quality alignments found. Could not determine best reference.')

  # Find the right reference file and align to it.
  ref_path = get_ref_path(best_seq, seqs_to_refs, args.refs_dir)
  logging.warning(f'Found {ref_path.name} to be the best reference choice.')
  align(ref_path, args.fastq1, args.fastq2, out_path, threads=args.threads)

  if not args.keep_tmp:
    os.remove(align_tmp_path)


def get_reads_base(reads_path: Path) -> str:
  """Get the basename of a file, omitting any `_1` or `_2` before the extension.
  Takes a path, possibly including directories, and returns a string of the whole
  path, but with the file extension and any `_1` or `_2` removed."""
  basename = reads_path.stem
  if basename.endswith('_1') or basename.endswith('_2'):
    basename = basename[:-2]
  return str(reads_path.parent / basename)


def get_tmp_path(base: str, ext: str) -> Path:
  tries = 1
  name = base+'.'+ext
  while os.path.exists(name):
    tries += 1
    nonce = random.randint(100000, 999999)
    name = f'{base}.{nonce}.{ext}'
    if tries % 1000 == 0:
      logging.warning(
        f'Warning: Problem finding tmp filename. Tried {tries} times. Example: {name!r} exists'
      )
  return Path(name)


def align(ref_path: Path, fq1_path: Path, fq2_path: Path, out_path: Path, threads: int=1) -> None:
  cmd_raw = ('align.py', 'bwa', '--threads', threads, ref_path, fq1_path, fq2_path, '-o', out_path)
  cmd = list(map(str, cmd_raw))
  logging.warning('+ $ '+' '.join(cmd))
  subprocess.run(cmd)


def read_map(map_path: Path) -> Tuple[Dict[str,str], Dict[str,int]]:
  seqs_to_refs = {}
  seq_sizes = {}
  with map_path.open() as map_file:
    for line in map_file:
      fields = line.rstrip('\r\n').split('\t')
      seq_name, ref_path, size_str = fields
      size = int(size_str)
      seqs_to_refs[seq_name] = ref_path
      seq_sizes[seq_name] = size
  return seqs_to_refs, seq_sizes


def count_alignments(align_path: Path, mapq_thres: int=None) -> collections.Counter:
  aln_counts = collections.Counter()
  proc = subprocess.Popen(('samtools', 'view', align_path), text=True, stdout=subprocess.PIPE)
  for aln in samreader.read(proc.stdout):
    if is_good_alignment(aln, mapq_thres):
      aln_counts[aln.rname] += 1
  return aln_counts


def get_best_seq(
    aln_counts: collections.Counter, seq_sizes: Dict[str,int], min_size: int=None
) -> str:
  for seq_name, count in sorted(aln_counts.items(), reverse=True, key=lambda item: item[1]):
    if min_size is None or seq_sizes[seq_name] >= min_size:
      return seq_name
  return None


def get_ref_path(best_seq: str, seqs_to_refs: Dict[str,str], refs_dir: Path=None) -> Path:
  try:
    ref_path_rel = seqs_to_refs[best_seq]
  except KeyError:
    fail(f'Did not find sequence {best_seq} in mapping file.')
  if refs_dir:
    ref_path = refs_dir / ref_path_rel
  else:
    ref_path = Path(ref_path_rel)
  if not ref_path.is_file():
    fail(f'Reference file not found or not a regular file: {ref_path}')
  return ref_path


def is_good_alignment(
    aln: samreader.Alignment, mapq_thres: int=None, bad_flags: int=BAD_FLAG_INT
) -> bool:
  if aln is None:
    return False
  if aln.rname is None:
    return False
  if aln.flag & bad_flags:
    return False
  if mapq_thres is not None and aln.mapq < mapq_thres:
    return False
  return True


def fail(message: str, exit_code: int=1) -> None:
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
