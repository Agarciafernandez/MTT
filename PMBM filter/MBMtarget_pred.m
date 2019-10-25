function filter_pred=MBMtarget_pred(filter_upd,F,Q,p_s,means_b,P_ini,e_ini,Ncom_b,k)
%Author: Angel F. Garcia-Fernandez

%Same as PoissonMBMtarget_pred.m, but removing Poisson part and adding new
%Bernoulli components according to the birth model



%Prediction for Bernoulli components
filter_pred.globHyp=[filter_upd.globHyp,ones(size(filter_upd.globHyp,1),Ncom_b)];
filter_pred.globHypWeight=filter_upd.globHypWeight;


Ntracks=length(filter_upd.tracks);

%Prediction of previous tracks (Bernoulli components)
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

%We add the new tracks (Bernoulli components)

for j=1:Ncom_b
    filter_pred.tracks{Ntracks+j}.meanB{1}=means_b(:,j);
    filter_pred.tracks{Ntracks+j}.covB{1}=P_ini;
    filter_pred.tracks{Ntracks+j}.eB=e_ini(j);
    filter_pred.tracks{Ntracks+j}.t_ini=k;
    filter_pred.tracks{Ntracks+j}.aHis{1}=j;
    
end
