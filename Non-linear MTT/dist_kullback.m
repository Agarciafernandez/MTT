function d=dist_kullback(mu0,var0,mu1,var1)
%Calculates the KLD between two Gaussians. State space of N dimensions.
%Author: Angel F. Garcia-Fernandez

N=length(mu0);
inv1=inv(var1);
sub=mu1-mu0;
d=1/2*(log(det(var1)/det(var0))+trace(inv1*var0)+sub'*inv1*sub-N);