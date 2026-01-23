function power_pmb=CalculatePowerPMB(filter_upd,omega)

%Author: Ángel García-Fernández


%This function approximates the power of a PMB as a PMB, see 
%Á. F. García-Fernández and G. Battistelli, "Distributed Poisson Multi-Bernoulli Filtering via Generalized Covariance Intersection," in IEEE Transactions on Signal Processing, vol. 74, pp. 246-257, 2026, doi: 10.1109/TSP.2026.3651805.


power_pmb.globHyp=filter_upd.globHyp;
power_pmb.globHypWeight=filter_upd.globHypWeight;


%We compute the power of the PPP intensity
power_pmb.meanPois=filter_upd.meanPois;
Ncom=length(filter_upd.weightPois);
kappa_list=zeros(1,Ncom);
for i=1:Ncom
     power_pmb.covPois(:,:,i)=filter_upd.covPois(:,:,i)/omega;
     kappa_list(i)=KappaPowerGaussian(filter_upd.covPois(:,:,i),omega);
end

power_pmb.weightPois=kappa_list.*(filter_upd.weightPois).^omega;

%Power for Bernoulli components

Ntracks=length(filter_upd.tracks);

if(Ntracks>0)
    for i=1:Ntracks
        Nhyp_i=length(filter_upd.tracks{i}.eB); 

        %Nhyp_i must be equal to one since this is PMB filtering
        if (Nhyp_i~=1)
            warning('CalculatePowerPMB function should only be used for PMB densities, not PMBM.')
        end
        
        %We do not change the t_ini metadata
        power_pmb.tracks{i}.t_ini=filter_upd.tracks{i}.t_ini;
        
        for j=1:Nhyp_i %We go through all hypotheses
            
            power_pmb.tracks{i}.meanB{j}=filter_upd.tracks{i}.meanB{j};
            power_pmb.tracks{i}.covB{j}=filter_upd.tracks{i}.covB{j}/omega;
            
            kappa=KappaPowerGaussian(filter_upd.tracks{i}.covB{j},omega);
            
            eB=filter_upd.tracks{i}.eB(j);

            eB=min(eB,1); %This is done to improve numerical accuracy

            power_pmb.tracks{i}.eB(j)=(eB^omega*kappa)/((1-eB)^omega+eB^omega*kappa);
            
            %We do not change the aHis metadata
            power_pmb.tracks{i}.aHis{j}=filter_upd.tracks{i}.aHis{j};
            
        end
        
    end
else

    power_pmb.tracks=cell(0,1);
    power_pmb.globHyp=[];
    power_pmb.globHypWeight=[];
    
end
