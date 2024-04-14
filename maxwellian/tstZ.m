% 24-04-13 14:05 test Zfun
close all; clear; clc;

z=0.1-1i;

Z1=zfun(z)
Z2=zfunJpole(z)
Z2-Z1