# Testing isotropization from UrQMD transport simulation of heavy ion collisions

This is the code used to make calculations for this paper: https://arxiv.org/pdf/1508.04378.
It mainly does the following:
1. Reading in particles from the UrQMD simulation of heavy ion collisions (see https://urqmd.org/)
2. Coarse-graining these particles, calculating energy-momentum tensor on a space-time grid
3. Boosting the energy-momentum tensor to Landau frame in every cell
4. Checking, how isotropic the energy momentum tensor is

The idea of this procedure is to check, if particles produces by UrQMD simulation can be
justifiably transformed into a coarse-grained ideal hydrodynamic picture. In simpler words,
"how far is UrQMD from ideal hydrodynamics". Of course, the answer depends on space, time, and many technical
factors, so see the paper for details --- https://arxiv.org/pdf/1508.04378.

Codewise, this is an old fortran 90 code. To compile just type "make". I don't really expect this code to be useful
to anybody, but for a sake of reference, conservation, and reproducibility I save it here.
