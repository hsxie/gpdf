Supplementary material for the MATLAB numerical routines used in GPDF paper.

Author: H. S. Xie, huashengxie@gmail.com, IFTS-ZJU
Title: Generalized Plasma Dispersion Function: One-Solve-All Treatment, Visualizations,
and Application to Landau Damping
Journal: Physics of Plasmas (or see http://arxiv.org/abs/1305.6476)

Files:
./gpdf/zetaph.m       -- main routine for calculate GPDF
./test/test_1d.m      -- 1d plot of GPDF, using zetaph.m
      /test_2d.m      -- 2d plot of GPDF, using zetaph.m
	  /faddeeva.m     -- for calculate original PDF, i.e., Maxwellian distribution
./root/rootfinding.m  -- for root finding, using zetaph.m
./ivp/vp_landau_sim.m -- initial value simulation code for benchmark GPDF


2013-05-28 17:56