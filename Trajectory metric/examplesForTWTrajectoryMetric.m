%This file provides a demo with one-dimensional target states to use
%the time weighted metric on sets of trajectories

% [1]   Á. F. García-Fernández, A. S. Rahmathullah and L. Svensson, "A time-weighted metric for sets of trajectories to assess multi-object
% tracking algorithms" in Proceedings of the 24th International Conference on Information Fusion, 2021.

% The metric on sets of trajectories without time weighting is proposed in

% [2] Á. F. García-Fernández, A. S. Rahmathullah and L. Svensson,"A Metric on the Space of Finite Sets of Trajectories for Evaluation of
%  Multi-Target Tracking Algorithms," in IEEE Transactions on Signal Processing, vol. 68, pp. 3917-3928, 2020

%Author: Ángel García-Fernández

clear


% parameters for the metric
c=5;
p=1;
gamma=10;
%gamma=10^8
%We use a rho forgetting factor for the metric
rho=0.995;

%User can chooose the ID of the estimate (from 1 to 4), see paper
ID_est=2;



% structures for X and Y. X will define the grounth truth set of
% trajectories. Y will define the estimated sets of trajectories.
% Meaning of the struct fields is explained in LPTrajMetric_cluster.m
X = struct('xState', [], 'tVec', [], 'iVec', []);

%Parameters to define sets of trajectories
delVal = 3;


%T is the total number of time steps
%stDim Target state dimension
%nx: number of trajectories in X
%ny: number of trajectories in Y
T = 800; stDim = 1; nx = 2;
x_far=100;
x_close=10;

X.xState = zeros(stDim, T, nx);
X.tVec = [1,1]';
X.iVec = [T,T]';

X.xState(stDim,1:100,2)=x_far;
X.xState(stDim,101:200,2)=x_close+(x_far-x_close)*(100:-1:1)/100;
X.xState(stDim,201:300,2)=x_close;
X.xState(stDim,301:400,2)=x_close+(x_far-x_close)*(1:1:100)/100;
X.xState(stDim,401:500,2)=x_far;
X.xState(stDim,501:600,2)=x_close+(x_far-x_close)*(100:-1:1)/100;
X.xState(stDim,601:700,2)=x_close;
X.xState(stDim,701:800,2)=x_close+(x_far-x_close)*(1:1:100)/100;


%Baseline estimate
Y1 = struct('xState', [], 'tVec', [], 'iVec', []);
ny = 2;
Y1.xState(1,:,1) = X.xState(1,:,1) - delVal;
Y1.xState(1,:,2) = X.xState(1,:,2) + delVal;

Y1.tVec = [1,1]';
Y1.iVec = [T,T]';


switch ID_est
    case 1
        %Estimate 1
        Y2=Y1;
    case 2
        %Estimate 2
        Y2=Y1;
        Y2.xState(1,250:800,2)=Y1.xState(1,250:800,1); %Track switch
        Y2.xState(1,250:800,1)=Y1.xState(1,250:800,2);
        
    case 3
        %Estimate 3
        Y2=Y1;
        Y2.xState(1,650:800,2)=Y1.xState(1,650:800,1); %Track switch
        Y2.xState(1,650:800,1)=Y1.xState(1,650:800,2);
    case 4
        %Estimate 4
        Y2 = Y1;
        Y2.xState(1,550:800,2)=x_far+delVal+(3*x_far-x_far)*(1:251)/151;
end




%Time weights for localisation, missed, false targets. These time weights
%are normalised to sum to one
time_weights1=(1-rho)/(1-rho^T)*rho.^(T-(1:T))';
%Time weights for switching
time_weights2=time_weights1(2:end);

%Plot scenario
figure(1)
clf
hold on
for i=1:nx
    start_time= X.tVec(i);
    end_time=X.iVec(i)+X.tVec(i)-1;
    plot(start_time:end_time, X.xState(1,start_time:end_time,i),'-b','Linewidth',1.3)
end

for i=1:ny
    start_time= Y2.tVec(i);
    end_time=Y2.iVec(i)+Y2.tVec(i)-1;
    plot(start_time:end_time, Y2.xState(1,start_time:end_time,i),'-r','Linewidth',1.3)
end

% for i=1:ny
%     start_time= Y3.tVec(i);
%     end_time=Y3.iVec(i)+Y2.tVec(i)-1;
%     plot(start_time:end_time, Y3.xState(1,start_time:end_time,i),'--g','Linewidth',1.3)
% end

hold off
grid on
xlabel('Time step')
ylabel('Target state (m)')
title(['Estimate ',num2str(ID_est)])

set(gca,'FontSize',15)


%title('Sets of trajectories (X in blue, Y in red)')


% LP metric and its decomposition into costs for localisation, missed
% targets, false targets and track switches (See Sec. IV.D in [2])
[dxy, loc_cost, miss_cost, fa_cost, switch_cost]=TimeWeightedLPTrajMetric_cluster(X, Y2, c, p, gamma,time_weights1,time_weights2);

%Trajectory metric without weights
[dxy_nw, loc_cost_nw, miss_cost_nw, fa_cost_nw, switch_cost_nw]=LPTrajMetric_cluster(X, Y2, c, p, gamma);



%Note loc_cost,miss_cost, fa_cost and switch_cost are already given to the p-th
%power
disp('Time weighted trajectory metric value')
disp(dxy)
disp('Time weighted metric decomposition: localisation error, missed target cost, false target cost, switching cost')
disp([nthroot(sum(loc_cost),p), nthroot(sum(miss_cost),p), nthroot(sum(fa_cost),p), nthroot(sum(switch_cost),p)])

%Trajectory metric without time weights

disp('Trajectory metric value (no time weights)')
disp(dxy_nw/T)
disp('Trajectory metric decomposition: localisation error, missed target cost, false target cost, switching cost')
disp([nthroot(sum(loc_cost_nw),p)/T, nthroot(sum(miss_cost_nw),p)/T, nthroot(sum(fa_cost_nw),p)/T, nthroot(sum(switch_cost_nw),p)/T])



%Overall metric decomposed against time steps
dxy_time=nthroot(loc_cost+miss_cost+fa_cost+[switch_cost;0],p);

dxy_nw_time=nthroot(loc_cost_nw+miss_cost_nw+fa_cost_nw+[switch_cost_nw;0],p);


figure(2)
plot(1:T,dxy_time,'black',...
    1:T,nthroot(loc_cost,p),'--red',...
    1:T,nthroot(miss_cost,p),'-.green',...
    1:T,nthroot(fa_cost,p),'-*blue',...
    1:T,nthroot([0;switch_cost],p),'-+m','Linewidth',1.3)

grid on
legend('Time-weighted total','Time-weighted localisation','Time-weighted missed','Time-weighted false','Time-weighted switch')
ylabel('TW-TM decomposition')
xlabel('Time step')
set(gca,'FontSize',15)


figure(3)
plot(1:T,dxy_nw_time/T,'black',...
    1:T,nthroot(loc_cost_nw,p)/T,'--red',...
    1:T,nthroot(miss_cost_nw,p)/T,'-.green',...
    1:T,nthroot(fa_cost_nw,p)/T,'-*blue',...
    1:T,nthroot([0;switch_cost_nw],p)/T,'-+m','Linewidth',1.3)

grid on
legend('Total','Localisation','Missed','False','Switch')
ylabel('TM decomposition')
xlabel('Time step')
set(gca,'FontSize',15)



