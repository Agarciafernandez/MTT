function  [F_k,b_k,Q_k]=CalculateTransition_linear_SDE(A,u,time_lag,Q_c)

%This function calculates the transition model (transition matrix, bias and
%covariance matrix) resulting from a linear SDE, see Section III.B in 
%A. F. García-Fernández, S. Sarkka, "Gaussian multi-target filtering with
%target dynamics driven by a stochastic differential equation", in IEEE
%Transactions on Signal Processing, 2025.

%See also C. Van Loan, “Computing integrals involving the matrix exponential,” IEEE Transactions on Automatic Control, vol. 23, no. 3, pp. 395–404, 1978


%Offset for the mean
n_x=size(A,1);
n_p=size(u,2);

H1=[A,u;
    zeros(n_p,n_x),zeros(n_p,n_p)]*time_lag;

expm_H1=expm(H1);

b_k=expm_H1(1:n_x,n_x+1);

%Offset for the covariance


H2=[-A, Q_c;...
    zeros(n_x), A']*time_lag;

expm_H2=expm(H2);

F2=expm_H2(n_x+1:2*n_x,n_x+1:2*n_x);
G1=expm_H2(1:n_x,n_x+1:2*n_x);

Q_k=F2'*G1;

F_k=expm(A*time_lag);









