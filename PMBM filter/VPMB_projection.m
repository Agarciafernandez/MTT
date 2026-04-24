function filter_upd_pmb=VPMB_projection(filter_upd,max_iter,threshold)

% This function obtains the best fitting PMB to the updated
% PMBM based on coordinate descent KLD minimisation with auxiliary variables.
% This projection gives rise to the variational PMB filter in [1] based on the most likely assignment
% This is equivalent to the zero-temperature case in the variational PMB filter in Section III.B in [2].
% An alternative derivation based on the use of auxiliary variables and
% coordinate descent minimisation of the Kullback-Leibler divergence is
% provided in [3].

% [1] Y. Xia, K. Granström, L. Svensson, M. Fatemi, Á. F. García-Fernández and J. L. Williams,
% "Poisson Multi-Bernoulli Approximations for Multiple Extended Object Filtering,"
% in IEEE Transactions on Aerospace and Electronic Systems, vol. 58, no. 2, pp. 890-906, April 2022

% [2] J. L. Williams, "An Efficient, Variational Approximation of the Best Fitting Multi-Bernoulli Filter,"
% in IEEE Transactions on Signal Processing, vol. 63, no. 1, pp. 258-273, Jan.1, 2015

% [3]  A. F. García-Fernández, Y. Xia, “Variational PMB filter via coordinate descent Kullback-Leibler divergence minimisation,” 
% in 29th International Conference on Information Fusion, 2026.

%Author: Ángel García-Fernández


%Weights of global hypotheses
weights=filter_upd.globHypWeight;

%Number of tracks
Ntracks=length(filter_upd.tracks);


if(length(weights)>1 && Ntracks>0)

    %We initiate the algorithm with the standard PMB_projection
    filter_upd_pmb=PMB_projection(filter_upd);



    %Number of global hypotheses
    Nglob_hyp=length(filter_upd.globHypWeight);

    %Initial permutation matrix for each global hypothesis
    %Perm_hyp=repmat([1:Ntracks]',1,Nglob_hyp);

    cost_min_hyp=Inf;

    for n=1:max_iter
        
        %display(n)

        %Calculation of the best permutation for each global hypothesis
        [Perm_hyp,cost_min_hyp_u]=BestPermutation(filter_upd,filter_upd_pmb,Ntracks,Nglob_hyp,weights);


        %Calculation of the PMB projection with this permutation
        filter_upd_pmb=PMB_projection_perm(filter_upd,Perm_hyp);


        if(abs(cost_min_hyp_u-cost_min_hyp)<threshold)
            %We stop the iteration
            break
        end

        cost_min_hyp=cost_min_hyp_u;

        if(n==max_iter)
            display('The VPMB projection has not reached convergence')
        end

    end



else
    %There is only one global hypotheses. It is already in PMB form.
    filter_upd_pmb=filter_upd;

end
end

function [Perm_hyp,cost_min_hyp]=BestPermutation(filter_upd,filter_upd_pmb,Ntracks,Nglob_hyp,weights)

Perm_hyp=zeros(Ntracks,Nglob_hyp);

%We calculate all the KLDs between the original Bernoullis and the
%approximated Bernoullis
KLD_table=cell(Ntracks,Ntracks);

for i=1:Ntracks %i goes through the Bernoullis of q (filter_upd_pmb)

    eB_i=filter_upd_pmb.tracks{i}.eB(1);
    meanB_i=filter_upd_pmb.tracks{i}.meanB{1};
    covB_i=filter_upd_pmb.tracks{i}.covB{1};

    for j=1:Ntracks %i goes through the Bernoullis of the PMBM in hypothesis p

        Nhyp_j=length(filter_upd.tracks{j}.eB);

        KLD_table{i,j}=zeros(Nhyp_j,1);


        for index_hyp=1:Nhyp_j
            eB_j=filter_upd.tracks{j}.eB(index_hyp);
            meanB_j=filter_upd.tracks{j}.meanB{index_hyp};
            covB_j=filter_upd.tracks{j}.covB{index_hyp};

            KLD_table{i,j}(index_hyp)=KLD_Bernoulli(eB_j,meanB_j,covB_j,eB_i,meanB_i,covB_i);

        end

    end
end


%Minimum cost for each global hypothesis
cost_min_hyp_series = zeros(Nglob_hyp,1);

for p=1:Nglob_hyp
    %We create a cost matrix for each global hypothesis with the KLDs
    %between the Bernoulli components
    cost_matrix=zeros(Ntracks,Ntracks);

   
    for i=1:Ntracks %i goes through the Bernoullis of q (filter_upd_pmb)


        for j=1:Ntracks %i goes through the Bernoullis of the PMBM in hypothesis p


            index_hyp=filter_upd.globHyp(p,j); %Hypothesis for track j in p global hypothesis
            if(index_hyp>0)
                %The KLD first considers the Bernoullis of the PMBM and then the Bernoullis of the PMB
                cost_matrix(i,j)=KLD_table{i,j}(index_hyp);
            else

                %This is the KLD when the Bernoulli of the PMBM has
                %probability of existence 0
                eB_i=filter_upd_pmb.tracks{i}.eB(1);
                cost_matrix(i,j)=KLD_Bernoulli(0,[],[],eB_i,[],[]);
            end


        end
    end

    %We solve the 2D assignmnet problem
    [assignment,~,cost_min_hyp_series(p)] = assign2D(cost_matrix);

    %cost_min_hyp_series(p)=trace(cost_matrix);

    Perm_hyp(:,p)=assignment;

end
cost_min_hyp=weights*cost_min_hyp_series;

end


function KLD = KLD_Bernoulli(r1,mean1,P1,r2,mean2,P2)
% Kullback-Leibler divergence between Bernoullis
% Result given in
% M. Fontana, Á. F. García-Fernández and S. Maskell, "Data-Driven Clustering and Bernoulli Merging for the Poisson Multi-Bernoulli Mixture Filter,"
% in IEEE Transactions on Aerospace and Electronic Systems, vol. 59, no. 5, pp. 5287-5301, Oct. 2023


%We do this operation to improve numerical accuracy
r1=min(r1,1-10^(-8));
r2=min(r2,1-10^(-8));

if(r1==0)
    if(r2==1)
        KLD=Inf;
    else
        KLD=log(1/(1-r2));
    end
elseif(r1==1)
    if(r2==0)
        KLD=Inf;
    else
        N=length(mean1);
        sub=mean2-mean1;
        inv_P2=inv(P2);

        KLD_gaussian=1/2*(log(det(P2)/det(P1))+trace(inv_P2*P1)+sub'*inv_P2*sub-N);
        KLD=r1*log(r1/r2)+r1*KLD_gaussian;

    end
else
    %r1>0 and r1<1

    %r1>0
    if(r2==0 || r2==1)
        KLD=Inf;
    else
        N=length(mean1);
        sub=mean2-mean1;
        inv_P2=inv(P2);

        KLD_gaussian=1/2*(log(det(P2)/det(P1))+trace(inv_P2*P1)+sub'*inv_P2*sub-N);
        KLD=(1-r1)*log((1-r1)/(1-r2))+r1*log(r1/r2)+r1*KLD_gaussian;

    end

end

end
