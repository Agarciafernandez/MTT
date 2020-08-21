function [X_estimate,t_b_estimate,length_estimate]=TPMBM_estimate1(filter_upd,existence_estimation_threshold)

%Author: Angel F. Garcia Fernandez


%Option 1: one picks global hypothesis with highest weight and then takes the
%estimates above a certain threshold
X_estimate=cell(0,1);
t_b_estimate=[];
length_estimate=[];

if(~isempty(filter_upd.globHyp))
    globHypWeight=filter_upd.globHypWeight;
    [m,index]=max(globHypWeight);
    
    HypMax=filter_upd.globHyp(index,:);
    
    index_output=1;
    
    for i=1:length(HypMax)
        hyp_i=HypMax(i);
        if(hyp_i>0)
            Existence=filter_upd.tracks{i}.eB(hyp_i);
            if(Existence>existence_estimation_threshold)
                
                X_estimate{index_output}=filter_upd.tracks{i}.meanB{hyp_i};
                t_b_estimate=[t_b_estimate,filter_upd.tracks{i}.t_b];
                length_estimate=[length_estimate,filter_upd.tracks{i}.length];
                index_output=index_output+1;
            end
        end
    end
end
