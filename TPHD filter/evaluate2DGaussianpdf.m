function output=evaluate2DGaussianpdf(value,mean,cov)

resta=value-mean;
quad=-0.5*resta'/cov*resta;
d=length(mean);

output=exp(quad-1/2*log(det(cov))-d*log(2*pi)/2);