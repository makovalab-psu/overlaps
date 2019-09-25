#!/usr/bin/env python3
import argparse
import logging
import os
import pathlib
import sys
import unittest
# Path hack to load modules from the parent directory.
script_dir = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(script_dir.parent))
import analyze

DESCRIPTION = """Run unit(ish) tests."""


def make_argparser():
  parser = argparse.ArgumentParser(description=DESCRIPTION)
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
  unittest.main()


class TestFactory(unittest.TestCase):

  @classmethod
  def make_tests(cls):
    for i, data in enumerate(cls.test_data):
      test_function = cls.make_test(**data)
      name = data.get('name', str(i+1))
      setattr(cls, 'test_'+name, test_function)


########## analyze.py ##########

class GetBinTest(TestFactory):

  @classmethod
  def make_test(cls, num=None, denom=None, total_bins=None, bin=None, left=None, right=None):
    def test(self):
      bin_result, left_result, right_result = analyze.get_bin(num, denom, total_bins=total_bins)
      self.assertEqual(bin, bin_result)
      self.assertEqual(left, left_result)
      self.assertEqual(right, right_result)
    return test

  test_data = (
    {'denom':91,  'total_bins':10, 'num':10, 'bin':0, 'left':0,  'right':10},
    {'denom':91,  'total_bins':10, 'num':11, 'bin':1, 'left':10, 'right':19},
    {'denom':91,  'total_bins':11, 'num':1,  'bin':0, 'left':0,  'right':9},
    {'denom':91,  'total_bins':11, 'num':9,  'bin':0, 'left':0,  'right':9},
    {'denom':91,  'total_bins':11, 'num':10, 'bin':1, 'left':9,  'right':17},
    {'denom':91,  'total_bins':11, 'num':83, 'bin':9, 'left':75, 'right':83},
    {'denom':91,  'total_bins':11, 'num':84, 'bin':10,'left':83, 'right':91},
    {'denom':91,  'total_bins':11, 'num':91, 'bin':10,'left':83, 'right':91},
    {'denom':63,  'total_bins':9,  'num':1,  'bin':0, 'left':0,  'right':7},
    {'denom':63,  'total_bins':9,  'num':7,  'bin':0, 'left':0,  'right':7},
    {'denom':63,  'total_bins':9,  'num':8,  'bin':1, 'left':7,  'right':14},
    {'denom':63,  'total_bins':9,  'num':56, 'bin':7, 'left':49, 'right':56},
    {'denom':63,  'total_bins':9,  'num':57, 'bin':8, 'left':56, 'right':63},
    {'denom':63,  'total_bins':9,  'num':63, 'bin':8, 'left':56, 'right':63},
    {'denom':48,  'total_bins':6,  'num':39, 'bin':4, 'left':32, 'right':40},
    {'denom':48,  'total_bins':6,  'num':40, 'bin':4, 'left':32, 'right':40},
    {'denom':48,  'total_bins':6,  'num':41, 'bin':5, 'left':40, 'right':48},
    {'denom':251, 'total_bins':10, 'num':51, 'bin':1, 'left':26, 'right':51},
    {'denom':251, 'total_bins':10, 'num':52, 'bin':2, 'left':51, 'right':76},
  )

GetBinTest.make_tests()


class BinOverlapTest(TestFactory):

  @classmethod
  def make_test(cls, overlap=None, readlen=None, total_bins=None, expected=None):
    def test(self):
      observed = analyze.bin_overlap(overlap, readlen, total_bins=total_bins)
      self.assertEqual(expected, observed)
    return test

  test_data = (
    {'readlen':91,  'total_bins':11, 'overlap':35, 'expected':[0, 0, 0, 0, 0, 0, 2, 9, 8, 8, 8]},
  )

BinOverlapTest.make_tests()



if __name__ == '__main__':
  main(sys.argv)
