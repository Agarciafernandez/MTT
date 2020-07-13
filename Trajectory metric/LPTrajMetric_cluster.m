function [dxy, loc_cost, miss_cost, fa_cost, switch_cost]=LPTrajMetric_cluster(X, Y, c, p, gamma)
% This function computes the LP metric between sets of trajectories defined
% in
% Á. F. García-Fernández, A. S. Rahmathullah and L. Svensson,"A Metric on the Space of Finite Sets of Trajectories for Evaluation of
%  Multi-Target Tracking Algorithms," in IEEE Transactions on Signal Processing, vol. 68, pp. 3917-3928, 2020


% -------------------------------------------------------------------------
% Input:
% X, Y: sets of trajctories which are structs as follows:
%   X.tVec: 'nx x 1' dimensional vector that has start times of the 'nx'
%       trajectories in 'X'.
%   X.iVec: 'nx x 1' dimensional vector that has the duration of the 'nx'
%       trajectories in 'X'.
%   X.xState: 'stDim x T x nx' dimensional matrix, where 'stDim' is the
%       state dimension, 'T' is the max length of the trajectories. The
%       states of trajectory 'ind', 'X.xState(:, :, ind)' has '0' values
%       outisde '[X.tVec(ind), X.tVec(ind)+X.iVec(ind)-1]'. Note that within the
%       window where X.xState is valid can have 'holes', with 'nan' values.
% c: >0, cut-off parameter
% p: >= 1, exponent parameter
% gamma: >0, track switch penalty
% -------------------------------------------------------------------------
% Output:
% dxy: Metric value
% loc_cost: localisation cost (to the p-th power) for properly detected targets over time of dimension 'T x 1'
% miss_cost: cost (to the p-th power) for missed targets over time, dimension 'Tx1'
% fa_cost: cost (to the p-th power) for false targets over time, dimension 'Tx1'
% switch_cost: cost (to the p-th power) for switches over time, dimension '(T-1)x1'
% -------------------------------------------------------------------------
% Authors: Abu S. Rahmatullah and Ángel F. García-Fernández
% This code uses clustering, see Sec. IV.D, to compute the solution.


%%%%%%%%%% Input Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx = size(X.xState, 3);
ny = size(Y.xState, 3);
T = size(X.xState, 2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(nx==0 && ny==0)
    dxy=0;
    loc_cost=zeros(T,1);
    miss_cost=zeros(T,1);
    fa_cost=zeros(T,1);
    switch_cost=zeros(T-1,1);
    return;
end



%%%%%%%%%% localisation cost computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DAB = locCostComp_v2(X, Y, c, p);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Clustering
G=zeros(nx,ny);
T_min_a=zeros(nx,ny); %Matrices with minimum and maximum times in which trajectories i and j can be associated
T_max_a=zeros(nx,ny);
for i=1:nx
    t_i_x=X.tVec(i);
    t_f_x=X.tVec(i)+X.iVec(i)-1;
    for j=1:ny
        t_i_y=Y.tVec(j);
        t_f_y=Y.tVec(j)+Y.iVec(j)-1;
        t_i_max=max(t_i_x,t_i_y);
        t_f_min=min(t_f_x,t_f_y);
        %Option 1
        G_ij=squeeze(DAB(i,j,t_i_max:t_f_min))>=c^p;
        isnan_ij=or(isnan(X.xState(1,t_i_max:t_f_min,i)),isnan(Y.xState(1,t_i_max:t_f_min,j)))'; %We need to check when any of the trajectories isnan in this interval
        sum_isnan_ij=sum(isnan_ij);
        G(i,j)=sum(G_ij)==max(t_f_min-t_i_max+1-sum_isnan_ij,0);
        t_min_ij=find(and(G_ij==0,isnan_ij==0)); %They can be associated if this value is zero
        if(~isempty(t_min_ij))
            T_min_a(i,j)= t_min_ij(1)+t_i_max-1;
            T_max_a(i,j)= t_min_ij(end)+t_i_max-1;
        else
            T_min_a(i,j)= T+1; % We put infeasible values
            T_max_a(i,j)= 0;
            
        end
    end
end


%Graph connectivity
G_con=not(G(1:nx,1:ny));
Adj_matrix=[eye(nx),G_con;G_con',eye(ny)];
%Clustering
r = fliplr(symrcm(Adj_matrix));
Clusters = {r(1)};

max_length_c=1;
for i = 2:length(r)
    if any(Adj_matrix(Clusters{end},r(i)))
        Clusters{end}(end+1) = r(i);
    else
        Clusters{end+1} = r(i);
    end
    if(length(Clusters{end})>max_length_c)
        max_length_c=length(Clusters{end});
    end
end



loc_cost=zeros(T,1);
miss_cost=zeros(T,1);
fa_cost=zeros(T,1);
switch_cost=zeros(T-1,1);
dxy=0;
for i=1:length(Clusters)
    Cluster_i=Clusters{i};
    if(length(Cluster_i)==1)
        %The calculations simplify if there is only one trajectory in the
        %cluster
        i_x=find(Cluster_i<=nx);
        if(isempty(i_x))
            %Then it is a trajectory in Y
            i_y=find(Cluster_i>nx);
            list_y=Cluster_i(i_y)-nx;
            t_axis=Y.tVec(list_y):Y.tVec(list_y)+ Y.iVec(list_y)-1;
            isrealY_i=~isnan(Y.xState(1, t_axis,list_y))';
            dxy_i=c^p/2*sum(isrealY_i);
            fa_cost(t_axis)=fa_cost(t_axis)+c^p/2*isrealY_i;
        else
            list_x=Cluster_i(i_x);
            t_axis=X.tVec(list_x):X.tVec(list_x)+ X.iVec(list_x)-1;
            isrealX_i=~isnan(X.xState(1,t_axis,list_x))';
            dxy_i=c^p/2*sum(isrealX_i);
            miss_cost(t_axis)=miss_cost(t_axis)+c^p/2*isrealX_i;
        end
    elseif(length(Cluster_i)==2)
        i_x=find(Cluster_i<=nx);
        list_x=Cluster_i(i_x);
        i_y=find(Cluster_i>nx);
        list_y=Cluster_i(i_y)-nx;
        dxy_i=sum(DAB(list_x,list_y,:));
        
        t_axis=X.tVec(list_x):X.tVec(list_x)+ X.iVec(list_x)-1;
        isnan_X=isnan(X.xState(1,:,list_x))';
        no_exist_X=ones(T,1);
        no_exist_X(t_axis)=0;
        
        t_axis=Y.tVec(list_y):Y.tVec(list_y)+ Y.iVec(list_y)-1;
        isnan_Y=isnan(Y.xState(1,:,list_y))';
        no_exist_Y=ones(T,1);
        no_exist_Y(t_axis)=0;
        
        DAB_i=squeeze(DAB(list_x,list_y,:));
        
        %Errors equal to c^p correspond to false and missed target costs
        index1=DAB_i==c^p;
        fa_cost(index1)=fa_cost(index1)+c^p/2;
        miss_cost(index1)=miss_cost(index1)+c^p/2;
        
        %Missed targets
        index2=and(DAB_i==c^p/2,or(isnan_X,no_exist_X));
        fa_cost(index2)=fa_cost(index2)+c^p/2;
        
        %False targets
        index3=and(DAB_i==c^p/2,or(isnan_Y,no_exist_Y));
        miss_cost(index3)=miss_cost(index3)+c^p/2;
        
        %Localisation cost
        index4=and(DAB_i>0,not(or(index1, or(index2,index3))));
        
        loc_cost(index4)=loc_cost(index4)+DAB_i(index4);
        
    else
        
        
        i_x=find(Cluster_i<=nx);
        list_x=Cluster_i(i_x);
        i_y=find(Cluster_i>nx);
        list_y=Cluster_i(i_y)-nx;
        
        
        X_i.tVec=X.tVec(list_x);
        X_i.iVec=X.iVec(list_x);
        
        Y_i.tVec=Y.tVec(list_y);
        Y_i.iVec=Y.iVec(list_y);
        
        %Calculate minimum and maximum times when we need to consider the
        %assignments
        t_min=min(min(T_min_a(list_x,list_y)));
        t_max=max(max(T_max_a(list_x,list_y)));
        
        
        tf_X=X_i.tVec+X_i.iVec-1;
        tf_Y=Y_i.tVec+Y_i.iVec-1;
        
        
        %We sum the costs outside the considered window
        miss_cost_i=zeros(size(miss_cost));
        fa_cost_i=zeros(size(fa_cost));
        isrealX_i=~isnan(X.xState(1,:,list_x));
        isrealY_i=~isnan(Y.xState(1,:,list_y));
        
        for j=1:length(X_i.tVec)
            t_axis=[X_i.tVec(j):t_min-1, t_max+1:tf_X(j)];
            miss_cost_i(t_axis)=miss_cost_i(t_axis)+c^p/2*squeeze(isrealX_i(1,t_axis,j))';
        end
        
        for j=1:length(Y_i.tVec)
            t_axis=[Y_i.tVec(j):t_min-1, t_max+1:tf_Y(j)];
            fa_cost_i(t_axis)=fa_cost_i(t_axis)+c^p/2*squeeze(isrealY_i(1,t_axis,j))';
        end
        
        X_i.xState=X.xState(:,t_min:t_max,list_x);
        Y_i.xState=Y.xState(:,t_min:t_max,list_y);
        
        ti_X=X_i.tVec;
        ti_Y=Y_i.tVec;
        
        
        X_i.tVec=max(ti_X-t_min+1,1);
        Y_i.tVec=max(ti_Y-t_min+1,1);
        
        X_i.iVec=X_i.iVec-max(t_min-ti_X,0)-max(tf_X-t_max,0);
        Y_i.iVec=Y_i.iVec-max(t_min-ti_Y,0)-max(tf_Y-t_max,0);
        
        T_i=t_max-t_min+1;
        DAB_i=DAB([list_x,nx+1],[list_y,ny+1],t_min:t_max);
        
        nx_i=length(list_x);
        ny_i=length(list_y);
        nxny_i=nx_i*ny_i;
        nxny2_i=(nx_i+1)*(ny_i+1);
        
        
        
        [dxy_i,loc_cost_i2, miss_cost_i2, fa_cost_i2, switch_cost_i2]=LP_metric_cluster(X_i,Y_i,DAB_i,nx_i,ny_i,nxny_i,nxny2_i,T_i,c,p,gamma);
        dxy_i=dxy_i+sum(miss_cost_i)+sum(fa_cost_i);
        
        
        %We add the false and missed target costs of the tails
        
        miss_cost_i(t_min:t_max)=miss_cost_i(t_min:t_max)+miss_cost_i2;
        fa_cost_i(t_min:t_max)=fa_cost_i(t_min:t_max)+fa_cost_i2;
        
        loc_cost(t_min:t_max)=loc_cost(t_min:t_max)+loc_cost_i2;
        switch_cost(t_min:t_max-1)=switch_cost(t_min:t_max-1)+switch_cost_i2;
        miss_cost=miss_cost+miss_cost_i;
        fa_cost=fa_cost+fa_cost_i;
        
    end
    dxy=dxy+dxy_i;
    
end
dxy = dxy.^(1/p);

end




function [dxy,loc_cost, miss_cost, fa_cost, switch_cost]=LP_metric_cluster(X,Y,DAB,nx,ny,nxny,nxny2,T,c,p,gamma)
%%%%%%%%%%  variables to be calculated in LP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x = [W_1,1(1) W_2,1(1) .. W_nx+1,1(1), .. W_1,ny+1(1) W_2,ny+1(1) ...
% W_nx+1,ny+1(1), W_1,1(T) W_2,1(T) .. W_nx+1,1(T), .. W_1,ny+1(T)
% W_2,ny+1(T) ... W_nx+1,ny+1(T) e(1) .. e(T-1) h_11(1) .. h_nx,ny(1) ...
% h_1,1(T-1) ... h_nx,ny(T-1)]'

%%% Length of the variable components in x
WtLen = nxny2*T; etLen = T-1;  htLen = nxny * (T-1);
nParam = WtLen + etLen + htLen; % total number of variables

%%% Position of the variable components in x
WtPos = (1:WtLen);  etPos = WtLen+(1:etLen);
htPos = WtLen + etLen + (1:htLen);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%  objective function f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = zeros(nParam, 1);
f(WtPos) = reshape(DAB, [WtLen,1]); % for vec(W(1)) to vec(W(T)), loc cost
f(etPos) = 0.5 * gamma^p * ones(T-1,1); %for e(1) to e(T-1), switch cost
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%  equality constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Constraint 1
index_x=repmat(1:T*ny,nx+1,1);
index_x=index_x(:);
index_y=zeros(T*ny*(nx+1),1);
index_rep=1:ny*(nx+1);
for i=1:T
    index_y(index_rep+(i-1)*length(index_rep))=(ny+1)*(nx+1)*(i-1)+index_rep;
end
Aeq1=sparse(index_x,index_y,1,ny*T,nParam);
beq1 = ones(ny*T, 1);





%%%% Constraint 2 %%%%
index_x=repmat(reshape(1:T*nx,nx,T),ny+1,1);
index_x=index_x(:);
index_y=repmat(1:nx+1:(nx+1)/nx*length(index_x),nx,1);
index_y2=repmat((0:nx-1)',1,size(index_y,2));
index_y=index_y2+index_y;
index_y=index_y(:);
Aeq2=sparse(index_x,index_y,1,nx*T,nParam);
beq2 = ones(nx*T, 1);
Aeq = [Aeq1; Aeq2]; beq = [beq1; beq2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%  upper and lower bound constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0 <= W < inf, 0 <= e < inf, 0 <= h < inf
lb = zeros(nParam, 1);
ub = inf(nParam, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%  inequality constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Constraint 1 %%%%
index_minus_x=1:T-1;
index_minus_y=WtLen+index_minus_x;
value_minus=-1*ones(1,T-1);
index_one_x=repmat(1:T-1,nxny,1);
index_one_x=index_one_x(:);
index_one_y=WtLen + etLen +1:WtLen + etLen + (T-2)*nxny +nxny;
value_one=ones(1,length(index_one_y));
A1=sparse([index_minus_x';index_one_x],[index_minus_y';index_one_y'],[value_minus';value_one']);


%%%% Constraint 2 %%%%
index_m1_x=1: nxny*(T-1);
index_m1_y=htPos;
index_1_x=index_m1_x;
index_y=repmat(1:nx+1:(nx+1)*(ny+1)*(T-1),nx,1);
index_y(:,ny+1:ny+1:end)=[];
index_1_y=index_y+repmat((0:nx-1)',1,size(index_y,2));
index_1_y=index_1_y(:);
index_2_x=index_1_x;
index_2_y=index_1_y+(nx+1)*(ny+1);
A3=sparse([index_1_x';index_m1_x';index_2_x'],[index_1_y;index_m1_y';index_2_y],[ones(length(index_1_y),1);-ones(length(index_m1_y),1);-ones(length(index_2_y),1)]);



%%%% Constraint 3 %%%%
A4=sparse([index_1_x';index_m1_x';index_2_x'],[index_1_y;index_m1_y';index_2_y],[-ones(length(index_1_y),1);-ones(length(index_m1_y),1);ones(length(index_2_y),1)]);


A = [A1; A3; A4];
b = sparse((T-1)+nxny*(T-1)+nxny*(T-1),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%  optimisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
linProgOptions = optimoptions('linprog', 'Display','off'); %If this line returns an error, it may be required to install Matlab optimization toolbox
[x, dxy] = linprog(f, A, b, Aeq, beq, lb, ub, [], linProgOptions);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%  Metric and the assignment values to be returned %%%%%%%%%%%%%%%
wMat = reshape(x(1:nxny2*T), [nx+1, ny+1, T]);
[loc_cost, miss_cost, fa_cost, switch_cost] ...
    = computeLocFalseMissedSwitchCosts(wMat, DAB, X, Y, c, p, gamma);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function locCostMat = locCostComp_v2(X, Y, c, p)
% function locCostMat = locCostComp(stMat, X, Y, c, p)
% computing the localisation cost at each time 't' for every (i,j)

tmpCost = c^p / 2; % cost for being unassigned
T = size(X.xState, 2); nx = size(X.xState, 3); ny = size(Y.xState, 3);
locCostMat  = zeros(nx+1, ny+1, T);
for t = 1:T
    for xind = 1:nx+1
        if (xind <= nx) % xind not dummy
            if (t>=X.tVec(xind)) && (t<=(X.tVec(xind)+X.iVec(xind)-1))
                % if X_xind exists at t
                for yind = 1:ny+1
                    if (yind <= ny) && (t >= Y.tVec(yind)) && ...
                            (t <= (Y.iVec(yind) + Y.tVec(yind) - 1))
                        % if Y_yind exists at t
                        locCostMat(xind,yind,t) = computeLocCostPerTime( ...
                            X.xState(:,t,xind), Y.xState(:,t,yind), c, p);
                        
                    else % yind does not exist or yind is dummy
                        locCostMat(xind,yind,t) = tmpCost;
                    end
                end
            else        % if X_xind does not exist at t
                for yind = 1:ny
                    if (t >= Y.tVec(yind)) && ...
                            (t <= (Y.iVec(yind) + Y.tVec(yind) - 1))
                        % if Y_yind exists at t
                        locCostMat(xind,yind,t) = tmpCost;
                    end
                end
            end
        else    % xind is dummy
            for yind = 1:ny
                if (t >= Y.tVec(yind)) && ...
                        (t <= (Y.iVec(yind) + Y.tVec(yind) - 1))
                    % if Y_yind exists at t
                    locCostMat(xind,yind,t) = tmpCost;
                end
            end
        end
    end
end
end

function d = computeLocCostPerTime(x, y, c, p)
if all(~isnan(x)) && all(~isnan(y))
    % neither x nor y has hole
    d = min(norm(x-y, p)^p,c^p);
elseif any(isnan(x) & ~isnan(y)) || any(~isnan(x) & isnan(y))
    % exactly one of x and y has hole
    d = c^p/2;
else
    d = 0;
end
end

function [loc_cost, miss_cost, fa_cost, switch_cost] ...
    = computeLocFalseMissedSwitchCosts(w_mat, locCostMat, X, Y, c, p, gamma)
% computing the localisation cost, swtiching cost and cost for missed and false
% targets.

tmp_cost = c^p / 2; % cost for being unassigned
T = size(X.xState, 2); nx = size(X.xState, 3); ny = size(Y.xState, 3);


if(nx>1 && ny>1)
    switch_cost = 0.5 * gamma^p * ...
        squeeze(sum(sum(abs(diff(w_mat(1:nx, 1:ny, :), 1, 3)))));
elseif(nx==0 || ny==0)
    switch_cost=zeros(T-1,1);
elseif(nx==1 && ny==1)
    switch_cost = 0.5 * gamma^p * ...
        squeeze(abs(diff(w_mat(1:nx, 1:ny, :), 1, 3)));
    
else
    switch_cost = 0.5 * gamma^p * ...
        squeeze(sum(abs(diff(w_mat(1:nx, 1:ny, :), 1, 3))));
end

loc_mask = zeros(size(w_mat));
miss_mask = zeros(size(w_mat));
fa_mask = zeros(size(w_mat));
fa_miss_mask = zeros(size(w_mat)); %Accounts for false and missed target costs that arise for a localisation cost of c^p

for t = 1:T
    % localisation and miss cost calculations
    for xind = 1:nx
        if (t>=X.tVec(xind)) && (t<=(X.tVec(xind)+X.iVec(xind)-1))
            % if X_xind exists at t
            for yind = 1:ny
                if (yind <= ny) && (t >= Y.tVec(yind)) && ...
                        (t <= (Y.iVec(yind) + Y.tVec(yind) - 1))
                    % if Y_yind exists at t
                    if all(~isnan(X.xState(:,t,xind))) && ...
                            all(~isnan(Y.xState(:,t,yind)))
                        % there is no hole in x or y at this time
                        %%% add to localisation cost at the time based on
                        %%% weight (unless the weight is c^p
                        
                        if(locCostMat(xind, yind, t)<2*tmp_cost)
                            loc_mask(xind, yind, t) = 1;
                        else
                            fa_miss_mask(xind, yind, t) = 1;
                        end
                    elseif any(isnan(X.xState(:,t,xind))) && ...
                            all(~isnan(Y.xState(:,t,yind)))
                        % there is a hole in x but no hole in y
                        fa_mask(xind, yind, t) = 1;
                    elseif all(~isnan(X.xState(:,t,xind))) && ...
                            any(isnan(Y.xState(:,t,yind)))
                        % there is no hole in x but a hole in y
                        miss_mask(xind, yind, t) = 1;
                    end
                else % yind does not exist
                    
                    %%% add to miss cost at the time
                    miss_mask(xind, yind, t) = 1;
                end
            end
            %%% add to miss cost at the time for yind = ny+1
            yind = ny+1;
            miss_mask(xind, yind, t) = 1;
        end
    end
    
    for yind = 1:ny
        if (t>=Y.tVec(yind)) && (t<=(Y.tVec(yind)+Y.iVec(yind)-1))
            % if Y_yind exists at t
            if all(~isnan(Y.xState(:, t, yind))) % no hole in y
                for xind = 1:nx
                    if ~((xind <= nx) && (t >= X.tVec(xind)) && ...
                            (t <= (X.iVec(xind) + X.tVec(xind) - 1)))
                        % if X_xind does not exist at t
                        
                        %%% add to fa cost at the time
                        fa_mask(xind, yind, t) = 1;
                    end
                end
            end
            
            %%% add to fa cost at the time for xind = nx+1
            xind = nx+1;
            fa_mask(xind, yind, t) = 1;
        end
    end
end

loc_cost = squeeze(sum(sum(locCostMat .* w_mat .* loc_mask, 1), 2));
miss_cost = tmp_cost * squeeze(sum(sum(w_mat .* miss_mask, 1), 2))+ tmp_cost * squeeze(sum(sum(w_mat .* fa_miss_mask, 1), 2));
fa_cost = tmp_cost * squeeze(sum(sum(w_mat .* fa_mask, 1), 2))+ tmp_cost * squeeze(sum(sum(w_mat .* fa_miss_mask, 1), 2));
end
