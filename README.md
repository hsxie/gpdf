# gpdf

Fast calculation of plasma dispersion functions for Maxwellian distribution and arbitrary distributions.

# For arbitrary distributions

A fast and accurate approach is provided, where FFT is used.

Ref and cite: Hua-Sheng Xie, Generalized plasma dispersion function: One-solve-all treatment, visualizations, and application to Landau damping,  Phys. Plasmas 20, 092125 (2013), doi: 10.1063/1.4822332.

# Maxwellian distribution (i.e., the most widely used)

The plasma dispersion function $Z(z)$ can also be calculated by J-pole Pade approximation (probably the quickest method), i.e.,

$Z(z)={\sum}_{j=1}^{J} \frac{b_j}{z-c_j}$.

For J=8, one can use 

b1= -0.017340112270401 - 0.046306439626294i; 

b2= -0.739917811220052 + 0.839518284620274i; 

b3= 5.840632105105495 + 0.953602751322040i; 

b4= -5.583374181615043 -11.208550459628098i; 

c1= 2.237687725134293 - 1.625941024120362i; 

c2= 1.465234091939142 - 1.789620299603315i; 

c3= 0.839253966367922 - 1.891995211531426i; 

c4= 0.273936218055381 - 1.941787037576095i; 

b(5:8)=(b(1:4))*; 

c(5:8)=-(c(1:4))*,

where * denotes complex conjugate. The above approximation is valid for the upper plane. For $Im(z)$ weakly close to the real axis, i.e., not far from the real axis of the lower plane, the above approximation is also valid. For accurate calculation of $Im(z)<0$, one can use

$Z(z)=[Z(z*)]*+2i\sqrt{\pi}\exp(-z^2)$.

Other J, such as J=4,8,12,16,24 can also be found.

Ref and cite: Huasheng Xie, BO: A unified tool for plasma waves and instabilities analysis, Computer Physics Communications 244 (2019) 343–371.
