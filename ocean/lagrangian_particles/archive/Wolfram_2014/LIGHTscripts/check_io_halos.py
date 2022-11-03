#!/usr/bin/env python

import numpy as np
import glob
import re
import subprocess

def numberarray_between_strings(file, start, end, N=10000):
  #grep1 = subprocess.Popen(('grep', '-A%d'%(N), start, file),stdout=subprocess.PIPE)
  #grep2 = subprocess.check_output(('grep', '-B%d'%(N), end), stdin=grep1.stdout)
  #grep1.wait()
  #numbers = np.asarray(re.sub(end,'',re.sub(start,'',grep2)).split(),dtype='int')
  f = open(file,'r')
  data = f.read()
  numbers = np.asarray(re.findall(start + r'(.*?)' + end, data, re.DOTALL)[0].split(),dtype='int')
  f.close()
  return numbers

def main():
  files = glob.glob('log*err')
  files.sort()

  # get total number of processors
  procs = []
  for afile in files:
    number = int(re.findall('[0-9].*[0-9]',afile)[0])
    procs.append(number)
  procs = np.asarray(procs)
  procs.sort()

  assert procs[-1]+1 == procs.shape[0], 'Output file for log*err is missing!'
  nprocs = procs.shape[0]
  halos = np.zeros((nprocs,nprocs), dtype='int')
  for afile in files:
    print 'processing ', afile, 
    number = int(re.findall('[0-9].*[0-9]',afile)[0])
    initalhalo= numberarray_between_strings(afile, 'ioProcNeighs = ', 'ioProcNeighs = Done')
    computedhalo = numberarray_between_strings(afile, 'ioProcNeighs after = ', 'ioProcNeighs after = Done')
    print 'difference = ',np.max(initalhalo-computedhalo)


if __name__ == "__main__":
  main()
