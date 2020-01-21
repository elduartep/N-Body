# N-Body
Particle-mesh code for cosmological n-body simulations

# First You need to install

FFTW3

openmp

# INPUTS

## A CAMB power spectrum

## or a set of initial conditions in GADGET units
An ASCII file with 6 columns: x, y, z: comoving coordinates in Mpc/h, u, v, w: hte peculiar velocity in km/s times sqrt(a)

# OUTPUTS

xyz.dat
ASCII file with comoving x, y z in Mpc/h

uvw.dat
ASCII file with comoving vx, vy, vz in km/s

# parameters
parameters.h

There you need to specify the number of threads for paralelization, the initial and final z, the cosmological parameters, the CAMB file for computing the inicial conditions, a file with the initial conditions in GADGET units if you used 2LPT for generating the IC, the box lenght, number of particles and number of cells for computing the FFT



