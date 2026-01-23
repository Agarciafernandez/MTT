function X_estimate_agents=DPMBMStateEstimation(filter_agents,Nagents,type_estimator,existence_estimation_threshold1)

%Multi-target state estimation for distributed PMBM/PMB filtering
%X_estimate_sens is a cell returning the multi-target estimates for each
%agent

%Author: Ángel García-Fernández

X_estimate_agents=cell(Nagents,1);
for q=1:Nagents
    filter_upd=filter_agents{q}.filter_upd;

    switch type_estimator
        case 1
            X_estimate=PoissonMBMtarget_estimate1(filter_upd,existence_estimation_threshold1);
        case 2
            X_estimate=PoissonMBMtarget_estimate2(filter_upd);
        case 3
            X_estimate=PoissonMBMtarget_estimate3(filter_upd);
    end
    X_estimate_agents{q}=X_estimate;
end