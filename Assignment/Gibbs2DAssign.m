function     [assignments, costs]=Gibbs2DAssign(C,k_good)

% This function uses Gibbs sampling for data associations in a PMBM and PMB
% filter context to select k_good good data association hypotheses. 
% For an example of use in the update see the file "PoissonMBMtarget_update_Gibbs.m"

% It follows
% B. N. Vo, B. T. Vo, and H. G. Hoang, “An efficient implementation of the generalized labeled multi-Bernoulli filter,” IEEE Transactions on Signal
% Processing, vol. 65, no. 8, pp. 1975–1987, April 2017.
% and the code gibbswrap_jointpredupdt_custom.m in Prof. B. T. Vo webpage

% However, Gibbs sampling is used differently in PMBM/PMB filtering as we write the data association
% vector with Nm components, not with the number of previous targets.

% We start with the misdetection assignment

[Nm,Nt]=size(C); %Nm number of measurements, Nt number of potential targets (new and previous)
% Note that C is already shrinked to only consider the relevant columns and
% rows.

assignments= zeros(k_good,Nm);
costs= zeros(k_good,1)';

exp_C=exp(-C);

%gamma value of gamma (association) for each measurement
gamma= Nt-Nm+1:Nt; % We associate missed detections for previous Bernoullis as the initial solution (all measurements associated to new targets)
assignments(1,:)= gamma;
costs(1)=sum(C(sub2ind(size(C),1:Nm,gamma)));

logical_all=true(1,Nm);


for sol= 2:k_good
    for var= 1:Nm %We go through all previous measurements
        tempsamp= exp_C(var,:); % We take the row of costs for the current association variable

        %We fix current assignments except for the one indicated by var (we set their probability to zero so that they cannot be sampled)
        logical_all_but_j=logical_all;
        logical_all_but_j(var)=false;
        gamma_i=gamma(logical_all_but_j);
        tempsamp(gamma_i)= 0;

        idxold= find(tempsamp>0);
        tempsamp= tempsamp(idxold);

        random_n=rand(1,1);
        list_c=[0;cumsum(tempsamp(:))/sum(tempsamp)];


        chosen_bin=find(list_c<random_n);
        gamma(var)=chosen_bin(end);
        gamma(var)= idxold(gamma(var));


    end
    assignments(sol,:)= gamma;
    costs(sol)= sum(C(sub2ind(size(C),1:Nm,gamma)));
end

[C_2,I,~]= unique(assignments,'rows');
assignments= C_2;
costs= costs(I);

assignments=assignments';
costs=costs';
