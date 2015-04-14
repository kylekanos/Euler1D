Euler1D
=======

Solves 1D Euler equations using RK2 integration with slope-limited piecewise linear reconstruction and either a Central Upwind flux or Kurganov-Tadmor flux (users choice, requires commenting rather than file IO flag). Output is to a *.curve file which is easily viewed in LLNL's VisIt software.

Initial and boundary conditions are hard-coded by user, who, presumably, knows what he/she is doing.
  
This code was written some time ago by someone else (to which I've done some modifications) and may be awfully hideous compared to some of my more recent work. Note that this comment reflects the actual appearance and not function. Obviously a PPM + HLLC method would be superior for resolution compared to the current PLM + CU/KT method employed.

Compilation
===========

This is Fortran code, so you need a Fortran compiler. I have only ever used gnu's Fortran compiler (gfortran) and Intel's Fortran compiler. These would be compiled as

    gfortran Euler1D.f90 -o Euler1D
    
and

    ifort Euler1D.f90 -o Euler1D
    
Use of other compilers should be fine, but consult your man page for compiler flag options that you may want to use (e.g., `-O3`)

Disclaimer
----------

I offer no guarantees that this works for you or is anything worth using for research. This was designed for simple 1D cases (e.g., the Sod shock tube) to learn about computational hydrodynamics. 
