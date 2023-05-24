function     [assignments, costs]=Gibbs2DAssign_nb(C,k_good,log_pdf_clutter)

% This function selects k_good good data association hypotheses for negative binomial
% clutter and a PMB prior using Gibbs sampling

% This function calculates the cost in log form


[Nm,Nt]=size(C); %Nm number of measurements, Nt number of potential targets (new and previous)
%Note that C was already shrinked to only consider the
%relevant columns and rows

assignments= zeros(k_good,Nm);
costs= zeros(k_good,1)';

%exp_C=exp(-C);

%gamma value of gamma (association) for each measurement
gamma= zeros(1,Nm); % We start with all measurements assigned to clutter
assignments(1,:)= gamma;

costs(1)=Calculate_log_cost(gamma,C,log_pdf_clutter);

logical_all=true(1,Nm);

for k= 2:k_good
    for j= 1:Nm %We go through all previous measurements

        tempsamp= [-C(j,:),0]; % We take the row of costs for the current association variable

        %We fix current assignments except for the one indicated by var (we set their probability to zero so that they cannot be sampled)
        logical_all_but_j=logical_all;
        logical_all_but_j(j)=false;
        gamma_fixed=gamma(logical_all_but_j);
        indices_assigned=gamma_fixed>0;
        gamma_fixed_assigned=gamma_fixed(indices_assigned);
        tempsamp(gamma_fixed_assigned)=-Inf;

        %Number of measurements associated to clutter in fixed data
        %associations
        Nc=Nm-length(gamma_fixed_assigned);
        log_pdf_clutter_Nc=log_pdf_clutter(Nc+1);  %index Nc+1 corresponds to Nc clutter measurements
        log_pdf_clutter_Nc_plus_1=log_pdf_clutter(Nc+2);  %index Nc+2 corresponds to Nc+1 clutter measurements



        tempsamp(end)=log_pdf_clutter_Nc_plus_1;
        %In the association vector for target assignments, we need to multiply the assigned values
        %to the density of the clutter with Nc+1 measurements (in log
        %scale, we sum)

        tempsamp(1:end-1)=tempsamp(1:end-1)+log_pdf_clutter_Nc;

        %We convert to natural scale (non-log scale)
        maximum=max(tempsamp);
        tempsamp=exp(tempsamp-maximum);
        tempsamp=tempsamp/sum(tempsamp);

        random_n=rand(1,1);
        list_c=[0;cumsum(tempsamp(:))];

        chosen_bin=find(list_c<random_n);
        gamma(j)=chosen_bin(end);


        if(gamma(j)==Nt+1)
            %It means gamma(j) is clutter and should be set to zero
            gamma(j)=0;
        end


    end
    assignments(k,:)= gamma;
    costs(k)= Calculate_log_cost(gamma,C,log_pdf_clutter);
end

[C_2,I,~]= unique(assignments,'rows');
assignments= C_2;
costs= costs(I);

assignments=assignments';
costs=costs';

end

function cost=Calculate_log_cost(gamma,C,log_pdf_clutter)

indices=find(gamma>0);

% Cost for the detected measurements
cost=sum(C(sub2ind(size(C),indices,gamma(indices))));
% Now we calculate the cost for the measurements associated to clutter

% Number of clutter measurements
Nc=sum(gamma==0);
log_density_c=log_pdf_clutter(Nc+1);
cost=cost-log_density_c; %The cost comes with a minus sign

end
