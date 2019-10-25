function [squared_LP_metric, LP_metric_loc, LP_metric_miss, LP_metric_fal, LP_metric_switch]=ComputeLP_metric_all_error(X_estimate,t_b_estimate, length_estimate,X_truth,t_birth,t_death,c_gospa,gamma_track_metric,k,Nx)



%%We use the LP Trajectory metric for all trajectories up to the current
%%time step
X_multi_abu_format_prov=reshape(X_truth(:,1:k),Nx,size(X_truth,1)/Nx,k);
X_multi_abu_format= permute(X_multi_abu_format_prov,[1,3,2]);


%We construct the structures of X and Y
% X = struct('xState', [], 'tVec', [], 'iVec', []);
% Y = struct('xState', [], 'tVec', [], 'iVec', []);

%First X

index_alive=logical(sum(X_truth(:,1:k)~=0,2));
index_alive=index_alive(1:4:end);
X_multi_alive=X_multi_abu_format(:,1:k,index_alive);

X.xState = X_multi_alive;
X.tVec = t_birth(index_alive)';
X.iVec = min((t_death(index_alive)-1)',k)-t_birth(index_alive)'+1;

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


[dxy, loc_cost, miss_cost, fa_cost, switch_cost] = LPTrajMetric_sparse_cluster(X, Y, c_gospa, 2, gamma_track_metric);

 
% [dxy2, wMat2, loc_cost2, miss_cost2, fa_cost2, switch_cost2] = LPTrajMetric_sparse(X, Y, c_gospa, 2, gamma_track_metric);
% 
%     
% 
% if(length(switch_cost)<k-1)
%     display('error')
% end
% 
% if(length(switch_cost2)<k-1)
%     display('error')
% end
% 
% if(abs(dxy-dxy2)>0.001)
%     display('error metric 1')
% end
% 
% if(sum(abs(loc_cost-loc_cost2)>0.001)>0)
%     display('error metric 2')
% end
% 
% if(sum(abs(miss_cost-miss_cost2)>0.001)>0)
%     display('error metric 3')
% end
% 
% if(sum(abs(fa_cost-fa_cost2)>0.001)>0)
%     display('error metric 4')
% end
% 
% if(dxy^2<sum(loc_cost)+sum(miss_cost)+sum(fa_cost)+sum(switch_cost)-10^(-5) || dxy^2>sum(loc_cost)+sum(miss_cost)+sum(fa_cost)+sum(switch_cost)+10^(-5))
%         display('error metric decomposition')
%   end
% 
% if(dxy2^2<sum(loc_cost2)+sum(miss_cost2)+sum(fa_cost2)+sum(switch_cost2)-10^(-5) || dxy2^2>sum(loc_cost2)+sum(miss_cost2)+sum(fa_cost2)+sum(switch_cost2)+10^(-5))
%         display('error metric decomposition')
%   end

squared_LP_metric=dxy^2/k; %normalised by k
LP_metric_loc=sum(loc_cost)/k;
LP_metric_miss=sum(miss_cost)/k;
LP_metric_fal=sum(fa_cost)/k;
LP_metric_switch=sum(switch_cost)/k;