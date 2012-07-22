#!/usr/bin/env python

from read_config import getConfigArguments
from os import system

cmd = "gcc -Wall -o covariances covariances.c ../../theory/theory_all.c ../../theory/cosmology.c ../../theory/tomo.c ../../theory/EB_functions.c ../../theory/maths.c ../../emu_11//emu.c ../../emu_11//hubble.c -lfftw3 -lgsl -lgslcblas -lm"  

system(cmd)

format = "nokeys"
config = "covariances.conf"
prog = "./covariances"

system("rm covariances" )
params = getConfigArguments(config,format)
print "# Combined Probes prediction code:"
print "# running " + prog + " " +params
system(prog + " " + params)
