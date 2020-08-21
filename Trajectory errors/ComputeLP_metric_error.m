function [squared_LP_metric, LP_metric_loc, LP_metric_mis, LP_metric_fal, LP_metric_switch]=ComputeLP_metric_error(X_estimate,t_b_estimate, length_estimate,X_truth,t_birth,t_death,c_gospa,gamma_track_metric,k,Nx)

%This function computes the squared error of the alive trajectories and its decomposition using the LP
%Trajectory metric (normalised by the current time window, which is from 1
%to k), to be able to compute Eq. (36) in
%Á. F. García-Fernández and L. Svensson, "Trajectory PHD and CPHD Filters," in IEEE Transactions on Signal Processing, vol. 67, no. 22, pp. 5702-5714, 15 Nov.15, 2019
%Author: Angel F. Garcia-Fernandez

X_multi_format_prov=reshape(X_truth(:,1:k),Nx,size(X_truth,1)/Nx,k);
X_multi_format= permute(X_multi_format_prov,[1,3,2]); %We use the correct format to call the function to compute the metric


%We construct the structures of X and Y

%First X
index_alive=and(k>=t_birth,k<t_death);
n_trajectories_alive=sum(index_alive);
X_multi_alive=X_multi_format(:,1:k,index_alive);

X.xState = X_multi_alive;
X.tVec = t_birth(index_alive)';
X.iVec = (repmat(k,1,n_trajectories_alive)-t_birth(index_alive)+1)';

%Second Y
Y.tVec=t_b_estimate';
Y.iVec=length_estimate';
Y.xState=zeros(Nx,k,length(X_estimate));

for i=1:length(X_estimate)
    Y.xState(:,t_b_estimate(i):t_b_estimate(i)+length_estimate(i)-1,i)=reshape(X_estimate{i},Nx,length_estimate(i));
end

%We only consider the position elements to compute the error
X.xState=X.xState([1,3],:,:);
Y.xState=Y.xState([1,3],:,:);

%We compute the LP metric and its decomposition
[dxy, loc_cost, miss_cost, fa_cost, switch_cost] = LPTrajMetric_cluster(X, Y, c_gospa, 2, gamma_track_metric);


squared_LP_metric=dxy^2/k; %We normalise by k
LP_metric_loc=sum(loc_cost)/k;
LP_metric_mis=sum(miss_cost)/k;
LP_metric_fal=sum(fa_cost)/k;
LP_metric_switch=sum(switch_cost)/k;