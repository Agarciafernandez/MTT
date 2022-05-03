function filter_pred=PMBMtarget_pred_spawning(filter_upd,F,Q,p_s,weights_b,means_b,covs_b,F1,F2,p_s1,p_s2,d_spawning)
%Author: Angel F. Garcia-Fernandez

%Same code as PoissonMBMtarget_pred.m, but with recyling the spawned components


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

globHyp=filter_upd.globHyp;
globHypWeight=filter_upd.globHypWeight;


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
            
            
            %Recycling of spawning components
            indeces_hyp_j=globHyp(:,i)==j;
            weight_hyp_j=sum(globHypWeight(indeces_hyp_j));
            eBj=filter_upd.tracks{i}.eB(j);
            
            mean_j=filter_upd.tracks{i}.meanB{j};
            cov_j=filter_upd.tracks{i}.covB{j};
            
            
            vel=sqrt(mean_j(2)^2+mean_j(4)^2);
            unit_vector=[mean_j(2)/vel,0,mean_j(4)/vel,0]';
            
            d_offset1=d_spawning*[-unit_vector(3),0,unit_vector(1),0]';
            d_offset2=-d_offset1; %Signed changed w.r.t. spawning mode 1
            
            %Spawning mode 1
            weight1=p_s1*eBj*weight_hyp_j;
            mean1=d_offset1+F1*mean_j;
            cov1=F1*cov_j*F1'+Q;
            
            %Spawning mode 2
            weight2=p_s2*eBj*weight_hyp_j;
            mean2=d_offset2+F2*mean_j;
            cov2=F2*cov_j*F2'+Q;
            
            
            
            filter_pred.weightPois=[filter_pred.weightPois,weight1,weight2];
            filter_pred.meanPois=[filter_pred.meanPois,mean1,mean2];
            
            filter_pred.covPois=cat(3,filter_pred.covPois,cov1,cov2);

            
            
            
        end
        
    end
    
    
        
        
        
        
    
    
    
else

    filter_pred.tracks=cell(0,1);
    filter_pred.globHyp=[];
    filter_pred.globHypWeight=[];
    
end
