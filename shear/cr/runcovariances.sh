!/bin/bash -xv
export LD_LIBRARY_PATH=/usr/local/lib

rm cov

gcc -Wall -o cov covariances.c ../../theory/survey.c ../../theory/halo.c ../../theory/theory_all.c ../../theory/tomo.c ../../theory/cosmology.c ../../theory/EB_functions.c ../../theory/maths.c ../../emu_11/emu.c ../../emu_11/hubble.c -lfftw3 -lgsl -lgslcblas -lm   

./cov 1 ../combined/paralist_om_s8_zoom_Ngrid100 predictions/COV_COYOTE_ 40 0 40 100 20 5000 1 1 5 3 0.09 1.99 96 z_histo_bin100_DES_BCC.txt 0 5000.0 10.0 0.42 5 0.09 0.3325 0.5025 0.6725 0.8725 1.99

