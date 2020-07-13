%This file provides a simple demo with one-dimensional target states to use the metric on sets of trajectories
% Á. F. García-Fernández, A. S. Rahmathullah and L. Svensson,"A Metric on the Space of Finite Sets of Trajectories for Evaluation of
%  Multi-Target Tracking Algorithms," in IEEE Transactions on Signal Processing, vol. 68, pp. 3917-3928, 2020

%Authors: Abu Sajana Rahmathullah and Ángel García-Fernández
clear

% Choose example number below
example_id=1; %Value from 1 to 7 (the first 6 are discussed in the paper)

% parameters for the metric
c=5; 
p=1; 
gamma=1;


% structures for X and Y. X will define the grounth truth set of
% trajectories. Y will define the estimated sets of trajectories.
% Meaning of the struct fields is explained in LPTrajMetric_cluster.m
X = struct('xState', [], 'tVec', [], 'iVec', []);
Y = struct('xState', [], 'tVec', [], 'iVec', []);

%Parameters to define sets of trajectories
delVal = 1;
DelVal = 0.1;
switch example_id
    case 1
        % Example 1
        %T is the total number of time steps
        %stDim Target state dimension
        %nx: number of trajectories in X
        %ny: number of trajectories in Y
        T = 5; stDim = 1; nx = 1; 
        X.xState = ones(stDim, T, nx); 
        X.tVec = 1; 
        X.iVec = T;
        ny = 1;
        Y.xState = X.xState + delVal; 
        Y.tVec = 1; 
        Y.iVec = T;
        
    case 2
        % Example 2
        T = 5; stDim = 1; nx = 1;
        X.xState = ones(stDim, T, nx); 
        X.tVec = [1]; 
        X.iVec = [T];
        ny = 1;
        Y.xState = ones(stDim, T, nx)+ delVal; Y.tVec = [1]; 
        Y.iVec = [T-1];
        
    case 3
        % Example 3
        T = 5; stDim = 1; nx = 1;
        X.xState = ones(stDim, T, nx); 
        X.tVec = [1]; 
        X.iVec = [T];
        ny  = 2;
        Y.xState = repmat(X.xState + DelVal, [1, 1, ny]); 
        Y.tVec = [1, 4]';
        Y.iVec = [3, 2]';
        
    case 4
        % Example 4
        T = 5; stDim = 1; nx = 3;
        X.xState = ones(stDim, T, nx); 
        X.tVec = [1, 1, 4]'; 
        X.iVec = [T, 3, 2]';
        X.xState(:, :, 2:3) = X.xState(:, :, 2:3)+ 2*DelVal+delVal;
        ny  = 3;
        Y.xState = repmat(X.xState(:,:,1)+DelVal, [1,1,ny]); Y.tVec = [1,4,1]';
        Y.xState(:, :, 3) = Y.xState(:, :, 3)+ delVal;
        Y.iVec = [3, 2, T]';
        
    case 5
        % Example 5
        T = 4; stDim = 1; nx = 2;
        X.xState = ones(stDim, T, nx); X.tVec = [1, 1]'; 
        X.iVec = [T, T]';
        X.xState(:, :, 2) = X.xState(:, :, 2)+ 2*DelVal+delVal;
        ny  = 2; 
        Tswi = 3;
        Y.xState = ones(stDim, T, nx); Y.tVec = [1, 1]'; Y.iVec = [T, T]';
        Y.xState(:, 1:Tswi-1, 1) = X.xState(:, 1:Tswi-1, 1)+ delVal;
        Y.xState(:, Tswi:T, 1) = X.xState(:, Tswi:T, 2) - delVal;
        Y.xState(:, 1:Tswi-1, 2) = X.xState(:, 1:Tswi-1, 2) - delVal;
        Y.xState(:, Tswi:T, 2) = X.xState(:, Tswi:T, 1)+ delVal;
        
    case 6
        % Example 6
        T = 4; stDim = 1; nx = 2;
        X.xState = ones(stDim, T, nx); X.tVec = [1, 1]'; X.iVec = [T, T]';
        X.xState(:, :, 2) = X.xState(:, :, 2)+ 2*DelVal+delVal;
        ny  = 3;
        Y.xState = ones(stDim, T, nx); Y.tVec = [1, 1, 3]'; Y.iVec = [T, 2, 2]';
        Y.xState(:, :, 1) = X.xState(:, :, 1)+ delVal;
        Y.xState(:, :, 2:3) = repmat(X.xState(:, :, 2)-delVal, [1 1 2]);
               
        
    case 7
        % Example where there can be holes in the trajectories, which are represented with nan symbol.
        T = 5; stDim = 1; nx = 2;
        X.xState = 0.5*ones(stDim, T, nx); X.tVec = [1, 1]'; X.iVec = [T, T]';
        X.xState(:, :, 2) = X.xState(:, :, 2)+ 5*DelVal+delVal;
        X.xState(:, 3, 2) = nan; %The second trajectory in X has a hole at time step 2. This is why it does not create track switch
        
        ny  = 2;
        Y.xState = ones(stDim, T, nx); Y.tVec = [1, 1]'; Y.iVec = [T, T]';
        Y.xState(:, :, 1) = X.xState(:, :, 1)+ delVal;
        Y.xState(:, :, 2) = zeros(stDim, T, 1);

        %Y.xState(:, :, 2:3) = repmat(X.xState(:, :, 2)-delVal, [1 1 2]);
                
end

%Plot scenario
figure(1)
clf
hold on
for i=1:nx
   start_time= X.tVec(i);
   end_time=X.iVec(i)+X.tVec(i)-1;
   plot(start_time:end_time, X.xState(1,start_time:end_time,i),'-ob')
end

for i=1:ny
   start_time= Y.tVec(i);
   end_time=Y.iVec(i)+Y.tVec(i)-1;
   plot(start_time:end_time, Y.xState(1,start_time:end_time,i),'-xr')
end

hold off
grid on
xlabel('Time step')
ylabel('Target state')
title('Sets of trajectories (X in blue, Y in red)')




% LP metric and its decomposition into costs for localisation, missed
% targets, false targets and track switches (See Sec. IV.D)
[dxy, loc_cost, miss_cost, fa_cost, switch_cost]=LPTrajMetric_cluster(X, Y, c, p, gamma);
%Note loc_cost,miss_cost, fa_cost and switch_cost are already given to the p-th
%power
disp('Metric value')
disp(dxy)
disp('Metric decomposition: localisation error, missed target cost, false target cost, switching cost')
disp([nthroot(sum(loc_cost),p), nthroot(sum(miss_cost),p), nthroot(sum(fa_cost),p), nthroot(sum(switch_cost),p)])

%Overall metric decomposed against time steps
dxy_time=nthroot(loc_cost+miss_cost+fa_cost+[0;switch_cost],p);

figure(2)
plot(1:T,dxy_time,'black',...
1:T,nthroot(loc_cost,p),'--red',...
1:T,nthroot(miss_cost,p),'-.green',...
1:T,nthroot(fa_cost,p),'-*blue',...
1:T,nthroot([0;switch_cost],p),'-+m','Linewidth',1.3)

grid on
legend('Total','Localisation','Missed','False','Switch')
ylabel('LP metric decomposition')
xlabel('Time step')



