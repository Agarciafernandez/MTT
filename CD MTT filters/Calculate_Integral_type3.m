function Integral=Calculate_Integral_type3(B,Q_c,A,delta_t)

%This function calculates an integral of type 3 as described in 
%A. F. García-Fernández, S. Sarkka, "Gaussian multi-target filtering with
%target dynamics driven by a stochastic differential equation", in IEEE
%Transactions on Signal Processing, 2025.

n=size(A,1);

Matrix=[-B,Q_c,zeros(n);...
    zeros(n),zeros(n),eye(n);...
    zeros(n),zeros(n),A];

Matrix_exp=expm(Matrix*delta_t);

exp_B=expm(B*delta_t);
H1=Matrix_exp(1:n,2*n+1:3*n);
Integral=exp_B*H1;

end