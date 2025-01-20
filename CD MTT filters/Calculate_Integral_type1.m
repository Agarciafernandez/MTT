function Integral=Calculate_Integral_type1(A,B,delta_t)

%This function calculates an integral of type 1 as described in 
%A. F. García-Fernández, S. Sarkka, "Gaussian multi-target filtering with
%target dynamics driven by a stochastic differential equation", in IEEE
%Transactions on Signal Processing, 2025.

[n,p]=size(B);
Matrix=[A,B;...
    zeros(p,n+p)];

Matrix_exp=expm(Matrix*delta_t);

Integral=Matrix_exp(1:n,n+1:n+p);


end