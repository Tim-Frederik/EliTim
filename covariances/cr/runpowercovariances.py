#!/usr/bin/env python

from read_config import getConfigArguments
from os import system

system("export LD_LIBRARY_PATH=/usr/local/lib")
system("rm power_covariances")


cmd = "gcc -Wall -o power_covariances power_covariances.c ../../theory/theory_all.c ../../theory/cosmology.c ../../theory/tomo.c ../../theory/survey.c ../../theory/halo.c ../../theory/EB_functions.c ../../theory/maths.c ../../emu_11/emu.c ../../emu_11/hubble.c -lfftw3 -lgsl -lgslcblas -lm"  

system(cmd)

format = "nokeys"
config = "power_covariances.conf"
prog = "./power_covariances"

params = getConfigArguments(config,format)
print "# Power Spectrum Covariance code:"
print "# running " + prog + " " +params
system(prog + " " + params)
#system("nohup "+ prog + " " + params+" /dev/null &")
