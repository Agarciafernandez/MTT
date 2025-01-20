function Integral=Calculate_Integral_type2(A,Q_c,delta_t)

%This function calculates an integral of type 2 as described in 
%A. F. García-Fernández, S. Sarkka, "Gaussian multi-target filtering with
%target dynamics driven by a stochastic differential equation", in IEEE
%Transactions on Signal Processing, 2025.

n=size(A,1);

Matrix=[-A, Q_c;...
    zeros(n), A'];

Matrix_exp=expm(Matrix*delta_t);
F2=Matrix_exp(n+1:2*n,n+1:2*n);
G1=Matrix_exp(1:n,n+1:2*n);

Integral=F2'*G1;

end