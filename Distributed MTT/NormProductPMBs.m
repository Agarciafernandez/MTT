function pmbm_fusion=NormProductPMBs(pmb_1,pmb_2,gating_threshold,Nhyp_max)

% This function computes the normalised product of two PMB densities
% The output is a PMBM density
% Á. F. García-Fernández and G. Battistelli, "Distributed Poisson Multi-Bernoulli Filtering via Generalized Covariance Intersection," in IEEE Transactions on Signal Processing, vol. 74, pp. 246-257, 2026, doi: 10.1109/TSP.2026.3651805.


% The Bernoulli components in pmb_2 are considered measurements that update
% the Bernoulli components in pmb_1

% Author: Ángel García-Fernández


Nx=size(pmb_1.meanPois,1);

%PPP intensity product

Np_1=length(pmb_1.weightPois);
Np_2=length(pmb_2.weightPois);

pmbm_fusion.weightPois=zeros(1,Np_1*Np_2);
pmbm_fusion.meanPois=zeros(Nx,Np_1*Np_2);
pmbm_fusion.covPois=zeros(Nx,Nx,Np_1*Np_2);


%We go through all PPP intensity components

for i=1:Np_1
    for j=1:Np_2

        mean_i=pmb_1.meanPois(:,i);
        mean_j=pmb_2.meanPois(:,j);

        cov_i=pmb_1.covPois(:,:,i);
        cov_j=pmb_2.covPois(:,:,j);


        [mean_prod,cov_prod,alpha_prod]=ProductTwoGaussians(mean_i,cov_i,mean_j,cov_j);
        pmbm_fusion.weightPois(i+Np_1*(j-1))=pmb_1.weightPois(i)*pmb_2.weightPois(j)*alpha_prod;

        pmbm_fusion.meanPois(:,i+Np_1*(j-1))=mean_prod;
        pmbm_fusion.covPois(:,:,i+Np_1*(j-1))=cov_prod;


    end
end


%New track generation in PMB1
Nprev_tracks=length(pmb_1.tracks); %Number of previous tracks in PMB1
Nnew_tracks=length(pmb_2.tracks); %New tracks in PMB1 (some may have existence probability equal to zero)

%Fusion of Bernoulli components in PMB1 with the PPP in PMB2
pmbm_fusion=NormProductPMBs_unassigned(pmb_1,pmb_2,pmbm_fusion,gating_threshold,Nprev_tracks,Nnew_tracks,Nx);

%Fusion of Bernoulli components in PMB2 with the PPP in PMB1
pmbm_fusion=NormProductPMBs_unassigned(pmb_2,pmb_1,pmbm_fusion,gating_threshold,0,Nprev_tracks,Nx);



%We update the tracks of PMB1 with the PMB2 components as measurements
%(assignment hypotheses
pmbm_fusion=NormProductPMBs_assigned(pmb_1,pmb_2,pmbm_fusion,gating_threshold,Nprev_tracks,Nnew_tracks);


%Update global hypotheses 


if(Nprev_tracks==0)
    %No measurement hypothesised to be received from a target yet (No
    %previous tracks)
    %We create one global hypothesis
    pmbm_fusion.globHypWeight=1;
    pmbm_fusion.globHyp=ones(1,Nnew_tracks);
    
    if(Nnew_tracks==0)
        pmbm_fusion.tracks=cell(0,1);
    end
        
else
    
    globWeightLog=[];
    globHyp=[];
    
    for p=1:length(pmb_1.globHypWeight)
        %We create cost matrix for each global hypothesis.
        cost_matrix_log=-Inf(Nprev_tracks+Nnew_tracks,Nnew_tracks);
        cost_misdetection=zeros(Nprev_tracks,1);
        
        for i=1:Nprev_tracks
            index_hyp=pmb_1.globHyp(p,i); %Hypothesis for track i in p global hypothesis
            Nhyp_i=length(pmb_1.tracks{i}.eB);
            %We generate the cost matrix for measurements
            
            if(index_hyp~=0)
                              
                index_max=length(pmbm_fusion.tracks{i}.weightBLog_k);
                indices=index_hyp+Nhyp_i*(1:Nnew_tracks);
                indices_c=indices<=index_max;
                indices_c=indices(indices_c);
                
                weights_log=pmbm_fusion.tracks{i}.weightBLog_k(indices_c);
                
                %We remove the weight of the misdetection hypothesis to use Murty (Later this weight is added).
                cost_matrix_log(i,1:length(indices_c))=weights_log-pmbm_fusion.tracks{i}.weightBLog_k(index_hyp);
                cost_misdetection(i)=pmbm_fusion.tracks{i}.weightBLog_k(index_hyp);
                              
                
            end
        end
        
        %New targets
        for i=Nprev_tracks+1:Nnew_tracks+Nprev_tracks
            weights_log=pmbm_fusion.tracks{i}.weightBLog_k;
            index=pmbm_fusion.tracks{i}.aHis{1};
            cost_matrix_log(i,index)=weights_log;
        end
        
                
        %We remove -Inf rows and columns for performing optimal assignment. We take them
        %into account for indexing later.
        %Columns that have only one value different from Inf are not fed
        %into Murty either.
               
        cos_matrix_inf=~isinf(cost_matrix_log);    
        indices_stay_column=find(sum(cos_matrix_inf,1)>1); %There is one after inequality

        indices_stay_row=find(sum(cos_matrix_inf(:,indices_stay_column),2)>0);    %There is a zero after inequality
        
        fixed_indices_stay_column=find(sum(cos_matrix_inf,1)==1); %These belong to the final hypotheses but are not fed into Murty's algorithms
        [fixed_indices_stay_row,~]=find(cos_matrix_inf(:,fixed_indices_stay_column));
        
      
        %cost_fixed is the overall cost of the hypotheses that always belong to the output 
        cost_fixed=sum(cost_matrix_log(fixed_indices_stay_row+size(cost_matrix_log,1)*(fixed_indices_stay_column'-1)));      
        cost_matrix_log_trimmed=cost_matrix_log(indices_stay_row,indices_stay_column);       
        
       
        
        globWeightLog_pred=log(pmb_1.globHypWeight(p));
        
        if(isempty(cost_matrix_log_trimmed))
            %One hypothesis: All targets are misdetected, no need for Murty                 
            opt_indices_trans=zeros(1,Nnew_tracks);
            opt_indices_trans(fixed_indices_stay_column)=fixed_indices_stay_row; 
            globWeightLog=[globWeightLog,sum(cost_misdetection)+cost_fixed+globWeightLog_pred];
      
        else
            %Number of new global hypotheses from this global hypothesis
            kbest=ceil(Nhyp_max*pmb_1.globHypWeight(p)); 
               
            
            %Option 1: We run Murty algorithm (making use of Hungarian algorithm to
            %solve the assginment problem). We call the function by using transpose and negative value
            %[opt_indices,nlcost]=murty(-cost_matrix_log_trimmed',kbest);
            
            %Option 2: We run MURTY algorithm (making use of Hungarian algorithm to
            %solve the assginment problem). We call the function by using transpose and negative value          
            %kBest2DAssign by David F. Crouse
            [opt_indices,~,nlcost]=kBest2DAssign(-cost_matrix_log_trimmed',kbest);
            
            %Option 3: Gibbs sampling            
            %display('We use Gibbs')
           %[opt_indices,nlcost]=Gibbs2DAssign(-cost_matrix_log_trimmed',kbest);

            
            opt_indices=opt_indices';
            nlcost=nlcost';
            
            
            %Optimal indices without removing Inf rows
            opt_indices_trans=Inf(size(opt_indices,1),Nnew_tracks);           
            for i=1:size(opt_indices,1)               
                opt_indices_trans(i,indices_stay_column)=indices_stay_row(opt_indices(i,:));
                
              %We add the single trajectory hypotheses that belong to all k
               %max
               opt_indices_trans(i,fixed_indices_stay_column)=fixed_indices_stay_row; 
            end
 
            globWeightLog=[globWeightLog,-nlcost+sum(cost_misdetection)+cost_fixed+globWeightLog_pred];
        end
        
        %We write the corresponding global hypothesis
        globHypProv=zeros(size(opt_indices_trans,1),Nprev_tracks+Nnew_tracks);
        
        for i=1:Nprev_tracks
            index_hyp=pmb_1.globHyp(p,i); %Hypothesis for track i in p global hypothesis
            Nhyp_i=length(pmb_1.tracks{i}.eB);            
            for j=1:size(opt_indices_trans,1)
                index_track= find(opt_indices_trans(j,:)==i);             
                if(isempty(index_track))
                    index=index_hyp;
                else
                    index=index_hyp+Nhyp_i*index_track;                   
                end
                globHypProv(j,i)=index;
            end
        end
        
        
        for i=Nprev_tracks+1:Nnew_tracks+Nprev_tracks
            
            for j=1:size(opt_indices_trans,1)
                index_track=find(opt_indices_trans(j,:)==i);
                if(isempty(index_track))
                    index=0;
                else
                    index=1;
                end       
                globHypProv(j,i)=index; %It is either 0 or 1
            end
        end
        
        globHyp=[globHyp;globHypProv];
        
   
    end
    pmbm_fusion.globHyp=globHyp;
    %Normalisation of weights of global hypotheses
    globWeight=exp(globWeightLog-max(max(globWeightLog)));
    globWeight=globWeight/sum(globWeight);
    pmbm_fusion.globHypWeight=globWeight;
      
end
