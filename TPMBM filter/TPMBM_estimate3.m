function [X_estimate,t_b_estimate,length_estimate]=TPMBM_estimate3(filter_upd)

%Author: Angel F. Garcia Fernandez

%We select the hypothesis (with fixed cardinality, i.e. considering MBM01 expansion) with highest
%weight

if(~isempty(filter_upd.globHyp))
    
    globHyp=filter_upd.globHyp;
    globHypWeight=filter_upd.globHypWeight;
    Nhyp=length(globHypWeight);
    N_tracks=length(filter_upd.tracks);
    
    
    weight_tot=globHypWeight;
    for j=1:Nhyp
        for i=1:N_tracks
            hyp=globHyp(j,i);
            if(hyp~=0)
                eB=filter_upd.tracks{i}.eB(hyp);
                if(eB>0.5)
                    weight_tot(j)=weight_tot(j)*eB;
                else
                    weight_tot(j)=weight_tot(j)*(1-eB);
                end
            end
        end
    end
    
    [m,index_hyp]=max(weight_tot);
    %We construct the estimate
    X_estimate=cell(0,1);
    t_b_estimate=[];
    length_estimate=[];
    index_output=1;
    
    for i=1:N_tracks
        hyp=globHyp(index_hyp,i);
        if(hyp~=0)
            eB=filter_upd.tracks{i}.eB(hyp);
            if(eB>0.5)
                
                X_estimate{index_output}=filter_upd.tracks{i}.meanB{hyp};
                t_b_estimate=[t_b_estimate,filter_upd.tracks{i}.t_b];
                length_estimate=[length_estimate,filter_upd.tracks{i}.length];
                
                index_output=index_output+1;
            end
        end
    end
    
else
    X_estimate=cell(0,1);
    t_b_estimate=[];
    length_estimate=[];
end