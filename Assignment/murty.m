function [assignments,costs]= murty(C,k_best)

%Author: Angel F. Garcia-Fernandez
%Implementation of [1] K. G. Murty, "An algorithm for ranking all the assignments in order of increasing cost", Oper. Res., vol. 16, no. 3, pp. 682–687, 1968.
%We consider non-square matrices.It just requires to consider an extra
%iteration in the foor loop
%Copyright (c) 2018, Angel F. Garcia-Fernandez


if k_best==0
    assignments=[];
    costs=[];
    return;
end

% Third party assignmentoptimal.c requires non-negative costs. We
% compensate the costs before returning the values
offset_cost=min(min(C));
C=C-offset_cost;


[n,m]= size(C);
assignments=zeros(k_best,n);
costs=zeros(1,k_best);

%Start with the optimal assignment
[assignment_jni,C_ini]=assignmentoptimal(C);
assignments(1,:)=assignment_jni';
costs(1)=C_ini;

if(k_best==1)
   costs=costs+(offset_cost.*sum(assignments>0,2))';
    return
end

%Initiation of lists Murty algorithm
N_list_max=max(n,10);

N_included=zeros(N_list_max,1);
N_excluded=zeros(N_list_max,1);
Included_rows=zeros(n,N_list_max);
Included_columns=zeros(n,N_list_max);
Excluded_rows=zeros(N_list_max,N_list_max);
Excluded_columns=zeros(N_list_max,N_list_max);
cost_nodes=zeros(N_list_max,1);
optimal_assignment_nodes=zeros(n,N_list_max);

%assignemnt_act denotes the assignment for which we perform the
%partitioning
assignment_act=assignment_jni';


for j=2:k_best
    
    %Positions where we can allocate the new nodes
    empty_indices=find(N_included+N_excluded==0);

    if(j==2)
        
        
        %First iteration
        for i=1:n-j+2 %In Murty's paper index goes up to n-1, as the cost matrix is square. IF it is not square, we go up to n.
                    
            Excluded_rows(1,i)=i;
            Excluded_columns(1,i)=assignment_act(i);
            
            Included_rows(1:i-1,i)=1:i-1;
            Included_columns(1:i-1,i)=assignment_act(1:i-1);
   
            N_included(i)=N_included(i)+i-1;
            N_excluded(i)=N_excluded(i)+1;
            
            %Calculation of the corresponding cost matrix
            C_prov=C;
            C_prov(i,assignment_act(i))=Inf;
            ind_row=true(n,1);
            ind_row(1:i-1)=false;
            ind_column=true(m,1);
            ind_column(assignment_act(1:i-1))=false;
            C_prov=C_prov(ind_row,ind_column);
            
            if(sum(sum(isinf(C_prov),2)==size(C_prov,2)))
                %A whole row is equal to Inf so the cost is Inf so we remove the node from the list                
                N_included(i)=0;
                N_excluded(i)=0;
                Excluded_rows(:,i)=0;
                Included_rows(:,i)=0;
                Excluded_columns(:,i)=0;
                Included_columns(:,i)=0;
                cost_nodes(i)=0;                
            else
          
              
                
                [opt_assign, C_ini] = assignmentoptimal(C_prov);
                cost_nodes(i)=sum(C(Included_rows(1:i-1,i)+size(C,1)*(Included_columns(1:i-1,i)-1)))+C_ini;
                
                %We should write the optimal assignment in terms of the original
                %C matrix (not C_prov)
                list_column=find(ind_column);
                optimal_assignment_nodes(Included_rows(1:i-1,i),i)=Included_columns(1:i-1,i);
                optimal_assignment_nodes(ind_row,i)=list_column(opt_assign); %To convert to the original matrix C
            end
            
        end
        
    else
            if(length(empty_indices)<N_new)
                %We need to allocated more memory to the variables
                N_new_list=10*N_new;
                N_list_max=N_list_max+N_new_list;
                
                N_included=[N_included;zeros(N_new_list,1)];
                N_excluded=[N_excluded;zeros(N_new_list,1)];
                
                Included_rows=[Included_rows,zeros(size(Included_rows,1),N_new_list)];
                Included_columns=[Included_columns,zeros(size(Included_columns,1),N_new_list)];
                
                 Excluded_rows_prov=Excluded_rows;
                Excluded_rows=zeros(N_list_max,N_list_max);
                Excluded_rows(1:size(Excluded_rows_prov,1),1:size(Excluded_rows_prov,2))=Excluded_rows_prov;
                
                Excluded_columns_prov=Excluded_columns;
                Excluded_columns=zeros(N_list_max,N_list_max);
                Excluded_columns(1:size(Excluded_columns_prov,1),1:size(Excluded_columns_prov,2))=Excluded_columns_prov;
                
                
                
                cost_nodes=[cost_nodes;zeros(N_new_list,1)];
                optimal_assignment_nodes=[optimal_assignment_nodes,zeros(size(optimal_assignment_nodes,1),N_new_list)];

                empty_indices=find(N_included+N_excluded==0);

            end
    
        
        %We perform a new partitioning
        for i=1:N_new
            
            
            empty_index_i= empty_indices(i);
            
            %We firstassignment_act_row_new(1:i-1) copy the parent node's excluded and included items
            Excluded_rows(1:N_excluded_act,empty_index_i)=Excluded_rows_act(1:N_excluded_act);
            Excluded_columns(1:N_excluded_act,empty_index_i)=Excluded_columns_act(1:N_excluded_act);
            Included_rows(:,empty_index_i)=Included_rows_act;
            Included_columns(:,empty_index_i)=Included_columns_act;
            
            %We add the new excluded node
            Excluded_rows(N_excluded_act+1,empty_index_i)=assignment_act_row_new(i);
            Excluded_columns(N_excluded_act+1,empty_index_i)=assignment_act_column_new(i);
            
            %We add the new included rows
            Included_rows(N_included_act+(1:i-1),empty_index_i)=assignment_act_row_new(1:i-1);
            Included_columns(N_included_act+(1:i-1),empty_index_i)=assignment_act_column_new(1:i-1);
            
            N_included(empty_index_i)=N_included_act+i-1;
            N_excluded(empty_index_i)=N_excluded_act+1;
            
            Included_rows_i=Included_rows(1:N_included(empty_index_i),empty_index_i);
            Included_columns_i=Included_columns(1:N_included(empty_index_i),empty_index_i);
            
            Excluded_rows_i=Excluded_rows(1:N_excluded(empty_index_i),empty_index_i);
            Excluded_columns_i=Excluded_columns(1:N_excluded(empty_index_i),empty_index_i);
            
            %Calculation of the corresponding cost matrix
            C_prov=C;
            C_prov(Excluded_rows_i++size(C,1)*(Excluded_columns_i-1))=Inf;
            ind_row=true(n,1);
            
            
            ind_row(Included_rows_i)=false;
            ind_column=true(m,1);
            ind_column(Included_columns_i)=false;
            C_prov=C_prov(ind_row,ind_column);
            
            
            if(sum(sum(isinf(C_prov),2)==size(C_prov,2)))
                %A whole row is equal to Inf so the cost is Inf so we remove the node from the list
                
                N_included(empty_index_i)=0;
                N_excluded(empty_index_i)=0;
                Excluded_rows(:,empty_index_i)=0;
                Included_rows(:,empty_index_i)=0;
                Excluded_columns(:,empty_index_i)=0;
                Included_columns(:,empty_index_i)=0;
                cost_nodes(empty_index_i)=0;
               
            else
             
                [opt_assign, C_ini] = assignmentoptimal(C_prov);
                
                
                cost_nodes(empty_index_i)=sum(C(Included_rows_i+size(C,1)*(Included_columns_i-1)))+C_ini;
                
                %We should write the optimal assignment in terms of the original
                %C matrix (not C_prov)
                list_column=find(ind_column);
                optimal_assignment_nodes(Included_rows_i,empty_index_i)=Included_columns(1:N_included(empty_index_i),empty_index_i);
                optimal_assignment_nodes(ind_row,empty_index_i)=list_column(opt_assign); %To convert to the original matrix C
                             
            end
            
        end
    
    end
    

    %We take the minimum
    list_nodes=find(cost_nodes>0);
    
    if(isempty(list_nodes))
        %There are no more possible assignments. We return the current
        %values
        assignments=assignments(1:j-1,:);
        costs=costs(1:j-1);  
        costs=costs+(offset_cost.*sum(assignments>0,2))';

        return;
   
    else
        
        [cost_min_act,index]=min(cost_nodes(list_nodes));
        index_node_act=list_nodes(index);
        
        assignments(j,:)=optimal_assignment_nodes(:,index_node_act);
        costs(j)=cost_min_act;
        
        %The node we consider for the next partitioning has these
        %characteristics
        assignment_act=assignments(j,:);
        
        
        Excluded_rows_act=Excluded_rows(:,index_node_act);
        Included_rows_act=Included_rows(:,index_node_act);
        Excluded_columns_act=Excluded_columns(:,index_node_act);
        Included_columns_act=Included_columns(:,index_node_act);
        N_included_act=N_included(index_node_act);
        N_excluded_act=N_excluded(index_node_act);
        
        %assignment_act_new (New nodes in assignment_act w.r.t. its node)
        N_new=n-N_included_act;
        assignment_act_row_new=find(Included_columns_act==0);
        assignment_act_column_new=assignment_act(assignment_act_row_new);
          
        
        %We remove this node from the list by setting
        N_included(index_node_act)=0;
        N_excluded(index_node_act)=0;
        Excluded_rows(:,index_node_act)=0;
        Included_rows(:,index_node_act)=0;
        Excluded_columns(:,index_node_act)=0;
        Included_columns(:,index_node_act)=0;
        cost_nodes(index_node_act)=0;
    end
    
end

%We compensate the costs
costs=costs+(offset_cost.*sum(assignments>0,2))';






