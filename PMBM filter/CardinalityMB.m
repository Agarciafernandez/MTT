function pcard=CardinalityMB(r)

%Two options to compute the cardinality distribution of a multi-Bernoulli
%RFS, one using FFT and the other one direct convolution

%We calculate the cardinality of a MB distribution using FFT
N = length(r);
exp_omega = exp(-1i*(0:N)/(N+1)*2*pi);
F = ones(1,N+1);
for i = 1:N
    F = F .* ((1-r(i)) + r(i)*exp_omega);
end
pcard = real(ifft(F));

%We calculate the cardinality of a MB distribution using direct convolution
% N = length(r);
% pcard = zeros(1,N+1);
% pcard(1) = 1;
% for i = 1:N
%   pcard(2:i+1) = (1-r(i))*pcard(2:i+1) + r(i)*pcard(1:i);
%   pcard(1) = (1-r(i))*pcard(1);
% end
