function [X_truth,t_birth,t_death]=TrajectoryWilliams15_qmetric(Scenario_number,Nsteps,F,numtruth,Q,Area)

%Code modified from the source files in https://arxiv.org/abs/1203.2995

if (Scenario_number == 1), % in this case, targets are present to begin with
    birthtime = zeros(1,numtruth);
    Pmid = diag([4,0.01,4,0.01]);
    
else, % in other cases, a new target appears every 10 time steps
    birthtime = 10*(0:numtruth-1);
    Pmid = 0.25*eye(4);
end;


%If one wants to change birth/death times, do it here
t_birth=[1,1,1,6];
t_death=[30,75,80,100];



simlen = Nsteps; % must be odd
midpoint = (simlen+1)/2;
numfb = midpoint-1;
chol_Q=chol(Q)';


% Initialise at time midpoint and propagate forward and backwards (I sum
% [Area(1)/2;0;Area(2)/2;0] to center it in the region of interest)
x = chol(Pmid)'*randn(size(F,1),numtruth)+repmat([Area(1)/2;0;Area(2)/2;0],1,numtruth);

%We add some different velocities
x(:,1)=x(:,1)+[0;-1;0;-1];
x(:,2)=x(:,2)+[0;1;0;-1];
x(:,3)=x(:,3)+[0;-1;0;1];
x(:,4)=x(:,4)+[0;1;0;1];

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
    
    



