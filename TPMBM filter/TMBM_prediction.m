function filter_pred=TMBM_prediction(filter_upd,F,Q,p_s,means_b,P_ini,e_ini,Lscan,Ncom_b,k,Nx)
%Author: Angel F. Garcia-Fernandez

%Same as TPMBM_prediction.m, but removing Poisson part and adding new
%Bernoulli components according to the multi-Bernoulli birth model




%Prediction for Bernoulli components
filter_pred.globHyp=[filter_upd.globHyp,ones(size(filter_upd.globHyp,1),Ncom_b)];
filter_pred.globHypWeight=filter_upd.globHypWeight;


Ntracks=length(filter_upd.tracks);

for i=1:Ntracks
    Nhyp_i=length(filter_upd.tracks{i}.eB);
    filter_pred.tracks{i}.t_ini=filter_upd.tracks{i}.t_ini;
    %We predict the birth time and length of the trajectory
    filter_pred.tracks{i}.t_b=filter_upd.tracks{i}.t_b;
    filter_pred.tracks{i}.length=filter_upd.tracks{i}.length+1;
    
    
    for j=1:Nhyp_i %We go through all hypotheses
        
        %Prediction of the mean
        filter_pred.tracks{i}.meanB{j}=[filter_upd.tracks{i}.meanB{j};F*filter_upd.tracks{i}.meanB{j}(end-Nx+1:end)];
        
        %Prediction of the covariance matrix
        length_k_i=filter_upd.tracks{i}.length;
        
        if(Lscan>1)
            min_length=min([Lscan-1,length_k_i]);
            cov_i=filter_upd.tracks{i}.covB{j}(end-Nx*min_length+1:end,end-Nx*min_length+1:end);
            cov_i_pred=F*cov_i(end-Nx+1:end,end-Nx+1:end)*F'+Q;
            F_mode=zeros(Nx,size(cov_i,1));
            F_mode(end-3:end,end-3:end)=F;
            cross=F_mode*cov_i;
            cov_upd_i=[cov_i,cross';cross,cov_i_pred];
        else
            cov_i=filter_upd.tracks{i}.covB{j};
            cov_upd_i=F*cov_i*F'+Q;
        end
        
        
        filter_pred.tracks{i}.covB{j}=cov_upd_i;
        filter_pred.tracks{i}.eB(j)=p_s*filter_upd.tracks{i}.eB(j);
        filter_pred.tracks{i}.aHis{j}=filter_upd.tracks{i}.aHis{j};
        
    end
    
end


%We add the new tracks (Bernoulli components)



for j=1:Ncom_b
    filter_pred.tracks{Ntracks+j}.meanB{1}=means_b(:,j);
    filter_pred.tracks{Ntracks+j}.covB{1}=P_ini;
    filter_pred.tracks{Ntracks+j}.eB=e_ini(j);
    filter_pred.tracks{Ntracks+j}.t_ini=k+1;
    filter_pred.tracks{Ntracks+j}.aHis{1}=j;
    filter_pred.tracks{Ntracks+j}.t_b=k+1;
    filter_pred.tracks{Ntracks+j}.length=1;

    
end

