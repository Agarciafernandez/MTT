function esfs=Compute_esf(z)
%This function computes all the elementary symmetric functions evaluated
%at z using Eq. (2.3) in 
%Halil Oruç, Hakan K. Akmaz,Symmetric functions and the Vandermonde matrix,Journal of Computational and Applied Mathematics,Volume 172, Issue 1,Pages 49-64,2004,

%Author: Ángel F. García-Fernández

esfs_prev=zeros(length(z)+1,1);
esfs_prev(1)=1;
esfs=zeros(length(z)+1,1);
esfs(1)=1;
for i=1:length(z)   
    for j=1:i-1
       esfs(j+1)=esfs_prev(j+1)+z(i)*esfs_prev(j);
    end
    %The case j=i can be computed directly as the product
    esfs(i+1)=prod(z(1:i));    
    %We continue the recursion
    esfs_prev=esfs;   
end
