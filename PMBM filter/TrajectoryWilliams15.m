function [X_truth,t_birth,t_death]=TrajectoryWilliams15(Scenario_number,Nsteps,F,numtruth,Q,Area)

%Code modified from the source files in https://arxiv.org/abs/1203.2995

if (Scenario_number == 1), % in this case, targets are present to begin with
    birthtime = zeros(1,numtruth);
    Pmid = 0.1*eye(4);
    t_birth=ones(1,numtruth);
    
else, % in other cases, a new target appears every 10 time steps
    birthtime = 10*(0:numtruth-1);
    Pmid = 0.25*eye(4);
    t_birth=1:10:10*numtruth;
end;

t_death=[(Nsteps+1)/2, (Nsteps+1)*ones(1,numtruth-1)];

%If one wants to change birth/death times, do it here
%t_birth=[1,5,1,7];
%t_death=[20,60,70,81];



simlen = Nsteps; % must be odd
midpoint = (simlen+1)/2;
numfb = midpoint-1;
chol_Q=chol(Q)';


% Initialise at time midpoint and propagate forward and backwards (I sum
% 100, 100 to center it in my region of interest)
x = chol(Pmid)'*randn(size(F,1),numtruth)+repmat([Area(1)/2;0;Area(2)/2;0],1,numtruth);
xf = x; xb = x;
X_truth=zeros(4*numtruth,Nsteps);

xlog{midpoint} = x;
X_truth(:,midpoint)=x(:);


for t = 1:numfb
    % Run forward and backward simulation process
    xf=F*xf + chol_Q*randn(size(F,1),size(x,2));
    xb = F\(xb + chol_Q*randn(size(F,1),size(x,2)));  
    xlog{midpoint-t} = xb(:,midpoint-t>birthtime);
    xlog{midpoint+t} = xf; % note that all targets exist after midpoint   
    X_truth(:,midpoint-t)=xb(:);
    X_truth(:,midpoint+t)=xf(:);
end

%We set zero the entries of X_truth that are not used 

for i=1:size(X_truth,1)/4
    X_truth_i=X_truth(4*i-3:4*i,:);
    X_truth_i(:,1:t_birth(i)-1)=0;
    X_truth_i(:,t_death(i):end)=0;   
    X_truth(4*i-3:4*i,:)=X_truth_i;
end
    
    



