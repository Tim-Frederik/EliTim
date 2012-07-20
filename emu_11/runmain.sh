#!/bin/bash -xv

rm main
gcc -Wall -o main main.c emu.c hubble.c -lfftw3 -lgsl -lgslcblas -lm
./main

