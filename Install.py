#!/usr/bin/env python

"""
Install.py tool to build the Tadah.MD library.
The library is required by the ML-TADAH package.

Issue below command from the lammps/src directory:

    make lib-tadah.md

"""

from __future__ import print_function
import sys, os, subprocess
from argparse import ArgumentParser

sys.path.append('..')
from install_helpers import get_cpus, fullpath

parser = ArgumentParser(prog='Install.py',
                        description="LAMMPS Tadah.MD library build wrapper script")

args = parser.parse_args()

build_dir='build'
cwd = fullpath('.')
lib = os.path.basename(cwd)

makefile=os.path.join(cwd,'Makefile')

if not os.path.exists(makefile):
  sys.exit("lib/%s/Makefile.%s does not exist" % (lib, machine))

# make the library with parallel make
n_cpus = get_cpus()

print("Configuring lib%s..." % lib)
cmd = "./configure"
try:
  txt = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
  print(txt.decode('UTF-8'))
except subprocess.CalledProcessError as e:
  print("Configure failed with:\n %s" % e.output.decode('UTF-8'))
  sys.exit(1)

print("Building lib%s..." % lib)
cmd = "make -f %s clean; make -f %s -j%d" % (makefile, makefile, n_cpus)
try:
  txt = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
  print(txt.decode('UTF-8'))
except subprocess.CalledProcessError as e:
  print("Make failed with:\n %s" % e.output.decode('UTF-8'))
  sys.exit(1)

if os.path.exists("%s/lib%s.a" % (os.path.join(cwd,build_dir),lib)):
  print("Build was successful.")
else:
  sys.exit("Build of lib/%s was NOT successful" % (lib, lib))
