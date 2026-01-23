function kappa=KappaPowerGaussian(cov,omega)
%Calculates the kappa parameter of the Power of a Gaussian, see Eq. (41) in
%[1] Battistelli, G., Chisci, L., Fantacci, C., Farina, A., and Graziano, A., "Consensus CPHD Filter for Distributed Multitarget Tracking", 
% IEEE Journal of Selected Topics in Signal Processing 7, 3 (2013), pp. 508-520.

%cov is the covariance matrix of the input Gaussian

cov_mult=2*pi*cov;
den=det(cov_mult)^(omega/2);
num=sqrt(det(cov_mult/omega));

kappa=num/den;


end