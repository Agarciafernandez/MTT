function X_estimate=PoissonMBMtarget_estimate1(filter_upd,existence_estimation_threshold)

%Author: Angel F. Garcia Fernandez 


%Option 1: one picks global hypothesis with highest weight and then takes the
%estimates above a certain threshold  
X_estimate=[];
if(~isempty(filter_upd.globHyp))
    globHypWeight=filter_upd.globHypWeight;
    [m,index]=max(globHypWeight);

    HypMax=filter_upd.globHyp(index,:);

    for i=1:length(HypMax)
        hyp_i=HypMax(i);
        if(hyp_i>0)
            Existence=filter_upd.tracks{i}.eB(hyp_i);
            if(Existence>existence_estimation_threshold)
       
                X_estimate_i=filter_upd.tracks{i}.meanB{hyp_i};
                X_estimate=[X_estimate;X_estimate_i];
            end
        end
    end
end
