!/bin/bash -xv
export LD_LIBRARY_PATH=/usr/local/lib

rm pred_shear_power

gcc -Wall -o pred_shear_power predictions_shear.c ../../theory/survey.c ../../theory/halo.c ../../theory/theory_all.c ../../theory/tomo.c ../../theory/cosmology.c ../../theory/EB_functions.c ../../theory/maths.c ../../emu_11/emu.c ../../emu_11/hubble.c -lfftw3 -lgsl -lgslcblas -lm   

./pred_shear_power ../paralist_WMAP7_1 ../pkappa_Halofit_ 1 0 1 100 20 5000 1 0 2 3 0.11 1.99 95 ../z_histo_bin95_larger_z01 0 5 0.09 0.3325 0.5025 0.6725 0.8725 1.99

