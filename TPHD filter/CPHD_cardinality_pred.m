function cardinality_pred=CPHD_cardinality_pred(cardinality_k,p_s,cardinality_birth)
%This function performs the prediction step for the cardinality of the CPHD
%filter.

%Cardinality distribution of the surviving targets
cardinality_surv_pred2=zeros(size(cardinality_k));
N_max=length(cardinality_surv_pred2)-1;
for j=0:N_max
    aux_v=j:N_max;
    %The use of gamma is a faster way to compute the factorial than
    %factorial gamma(n)=factorial(n-1)
    nchoosek_a=gamma(aux_v+1)./(gamma(j+1)*gamma(aux_v-j+1));  
    aux=(p_s^j)*nchoosek_a.*(1-p_s).^(aux_v-j).*cardinality_k(aux_v+1)';    
    cardinality_surv_pred2(j+1)=sum(aux);
end

%Cardinality distribution of the surviving targets
cardinality_surv_pred=zeros(size(cardinality_k));
N_max=length(cardinality_surv_pred)-1;
for j=0:N_max
    for p=j:N_max
        cardinality_surv_pred(j+1)=cardinality_surv_pred(j+1)+exp(sum(log(1:p))-sum(log(1:j))-sum(log(1:p-j))+j*log(p_s)+(p-j)*log(1-p_s))*cardinality_k(p+1);
    end
end

%Cardinality distribution with new born targets
cardinality_pred=zeros(size(cardinality_k));
for n=0:N_max    
    cardinality_pred(n+1)=sum(cardinality_surv_pred(1:n+1).*flipud(cardinality_birth(1:n+1)));
end

%Normalisation of the cardinality distribution (for robustness)
cardinality_pred=cardinality_pred/sum(cardinality_pred);