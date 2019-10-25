function filter_pred=PoissonMBMtarget_pred(filter_upd,F,Q,p_s,weights_b,means_b,covs_b)
%Author: Angel F. Garcia-Fernandez

%Prediction for Poisson component
filter_pred.weightPois=p_s*filter_upd.weightPois;
filter_pred.meanPois=F*filter_upd.meanPois;
Ncom=length(filter_upd.weightPois);

if(Ncom>0)
    for i=1:Ncom
        filter_pred.covPois(:,:,i)=F*filter_upd.covPois(:,:,i)*F'+Q;
    end
    
    %We add PHD of new born targets
    filter_pred.weightPois=cat(2,filter_pred.weightPois,weights_b);
    filter_pred.meanPois=cat(2,filter_pred.meanPois,means_b);
    filter_pred.covPois=cat(3,filter_pred.covPois,covs_b);
else
    filter_pred.weightPois=weights_b;
    filter_pred.meanPois=means_b;
    filter_pred.covPois=covs_b;
    
    
end


%Prediction for Bernoulli components
filter_pred.globHyp=filter_upd.globHyp;
filter_pred.globHypWeight=filter_upd.globHypWeight;


Ntracks=length(filter_upd.tracks);

if(Ntracks>0)
    for i=1:Ntracks
        Nhyp_i=length(filter_upd.tracks{i}.eB);
        filter_pred.tracks{i}.t_ini=filter_upd.tracks{i}.t_ini;
        
        for j=1:Nhyp_i %We go through all hypotheses
            
            filter_pred.tracks{i}.meanB{j}=F*filter_upd.tracks{i}.meanB{j};
            filter_pred.tracks{i}.covB{j}=F*filter_upd.tracks{i}.covB{j}*F'+Q;
            filter_pred.tracks{i}.eB(j)=p_s*filter_upd.tracks{i}.eB(j);
            filter_pred.tracks{i}.aHis{j}=filter_upd.tracks{i}.aHis{j};
            
        end
        
    end
else

    filter_pred.tracks=cell(0,1);
    filter_pred.globHyp=[];
    filter_pred.globHypWeight=[];
    
end
