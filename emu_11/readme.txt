Cosmic Emu v1.1

* NOW WITH REDSHIFT! *

This code gives predictions for the nonlinear matter power spectrum
for a given set of cosmological parameters (flat wCDM) .  The emulation 
process is covered in "The Coyote Universe III: Simulation Suite and 
Precision Emulator for the Nonlinear Matter Power Spectrum" by Lawrence et al.
The paper describes how power spectra from 37 simulated cosmologies
are used to construct a Gaussian process emulator for a spectrum given
inputs.  This code implements the Kriging equations for Gaussian
process predictions using parameter estimates obtained through the
process described in the paper.

The cosmological parameters and their ranges are:

0.120  < omega_m < 0.155 (including CDM and baryons)
0.0215 < omega_b < 0.0235
0.85   < n_s     < 1.05
0.61   < sigma_8 < 0.9
-1.30  < w       < -0.70
0      < z       < 1

The code will not produce an answer outside these ranges. The hubble
parameter for a given cosmology is determined from the best-fit value
for the distance to the last scattering surface and can be determined
with a separate function.

The code works in both C and Fortran (emu.c provides the base
function written in C and a wrapper available to Fortran).  The Gnu
Scientific Library (GSL) is required for compilation. In
addition to the GSL, make sure you have all of

corrlengths.h
design.h
emu.c
hubble.c
kemu.h
ksim.h
meansd.h
pcbasis.h
pcweights.h
precisions.h

A working program that demonstrates how to use everything is given
in main.c.   A makefile is provided that assumes the GSL is 
available and stored in a usual place.  If you have stored it 
somewhere else, then you're probably more qualified than I am to fix
the compilation so that it works.  To compile the C example, type

$ make

and run "emu.exe"  which will prompt you for an output file, 
parameter values for each input, and the desired type of output.
The header information contains the parameters and other
information calculated from the parameters.  The main output 
will be 2 columns with 1995 rows plus a header row.  The first column 
is k (in Mpc^{1} ) and the second column is the spectrum for the requested cosmology
at the requested redshift.

The main emulation function, emu,  takes 3 parameters.  The first is a vector
of length 6 with the cosmological parameters and redshift.  The second is a
placeholder vector of length 2*1995 which contains the k values and
the spectrume.  Thethird parameter is the output type where 0 indicates the 
log(Delta^2 / k^(1.5)), 1 indicates Delta^2, and 2 indicates P(k) =
2*(pi^2)*(Delta^2) / (k^3).

The function getH0fromCMB(double *x, double *stuff)
determines the H0 from the CMB constraint l_A=d_lss*pi/r_s=302.4 given
the cosmology specified in x.  This is the value used by the emulator.

Regarding the logo photo (CosmicEmu.jpg), the background comes from
NASA, ESA, and R. Massey at CalTech (
http://hubblesite.org/newscenter/archive/releases/2007/01/image/j/ ).
The emu is all over the Web and I don't know where it originated.  It
looks happy in space, I think.
