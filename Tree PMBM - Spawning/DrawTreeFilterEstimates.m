function DrawTreeFilterEstimates(X_truth,t_birth,t_death,father_spawning,X_estimate,...
    t_b_estimate,tree_index_estimate,genealogy_estimate,XLim,YLim,z,k)
%Author Ángel F. García-Fernández
%This function plots the output of the filter at each time step
%XLim and YLim set the limits for the representation

figure(4)
clf
axis([XLim(1) XLim(2),YLim(1) YLim(2)])
xlabel('x position (m)')
ylabel('y position (m)')
grid on
hold on

N_targets_tot=length(t_birth);

for i=1:N_targets_tot
    plot(X_truth(4*i-3,t_birth(i):t_death(i)-1),X_truth(4*i-1,t_birth(i):t_death(i)-1),'b','Linewidth',1.2)
    
    plot(X_truth(4*i-3,t_birth(i):10:t_death(i)-1),X_truth(4*i-1,t_birth(i):10:(t_death(i)-1)),'ob','Linewidth',1.2,'MarkerSize',3)
    
    if(father_spawning(i)==0)
        plot(X_truth(4*i-3,t_birth(i)),X_truth(4*i-1,t_birth(i)),'ob','Linewidth',1.2,'MarkerFaceColor','b','MarkerSize',3)
    else
        plot(X_truth(4*i-3,t_birth(i)),X_truth(4*i-1,t_birth(i)),'sr','Linewidth',1.2,'MarkerFaceColor','r','MarkerSize',4)     
    end
    
    text(X_truth(4*i-3,t_birth(i))-10,X_truth(4*i-1,t_birth(i))+10,num2str(t_birth(i)),'FontSize',8)
end

%Plot links for spawning parent-child
for i=1:N_targets_tot   
    if(father_spawning(i)>0)
        link_x=[X_truth(4*father_spawning(i)-3,t_birth(i)-1);...
            X_truth(4*i-3,t_birth(i))];
         link_y=[X_truth(4*father_spawning(i)-1,t_birth(i)-1);...
            X_truth(4*i-1,t_birth(i))];
        plot(link_x,link_y,'b','Linewidth',1.2)
    end
end

%We first plot all branches
for i=1:length(X_estimate)
     X_estimate_i=X_estimate{i};   
     plot(X_estimate_i(1:4:end),X_estimate_i(3:4:end),'-or')       
end


%Now, we need to link the branches that belong to the same tree according
%to the genealogy

unique_trees=unique(tree_index_estimate);
for i=1:length(unique_trees)   
    indices_branches=find(tree_index_estimate==unique_trees(i));
    
    for j=1:length(indices_branches)
        genealogy_ij=genealogy_estimate{indices_branches(j)};
        
        if(sum(genealogy_ij>1)>0)
            %This branch has spawned from another branch. If it does not
            %meet this condition, it is the main branch and we do not need
            %to do anything
            
            spawning_events=find(genealogy_ij>1);
            
            last_spawning_event=spawning_events(end);
            
            parent_branch_gen=genealogy_ij(1:last_spawning_event-1); %We take the genealogy variable cut just before the last spawning
            
            display(i)
            display(j)
            
            for p=1:j-1
                %We look for the parent branch, which should be one of the
                %previous branches
                genealogy_ip=genealogy_estimate{indices_branches(p)};

                if(parent_branch_gen==genealogy_ip(1:length(parent_branch_gen)))
                    %There is match. This is the parent branch
                    %We plot the link between branches
                    index_branch=indices_branches(p);
                    
                    X_estimate_parent=X_estimate{index_branch};
                    X_estimate_child=X_estimate{indices_branches(j)};
                    
                    %We need to calculate the state index in the parent
                    %branch (which is indexed w.r.t. the the previous last
                    %spawning, as the branch only contains information
                    %since the last spawning.
                    spawning_events_parent=find(parent_branch_gen>1);

                    if(isempty(spawning_events_parent))
                        index_state_parent=length(parent_branch_gen);
                    else
                        index_state_parent=length(parent_branch_gen)-spawning_events_parent(end)+1;
                    end
                    
                    
                    x_axis_link=[X_estimate_parent(4*index_state_parent-3),...
                        X_estimate_child(1)];
                    y_axis_link=[X_estimate_parent(4*index_state_parent-1),...
                              X_estimate_child(3)];

                    
                    
%                     x_axis_link=[X_truth(4*index_branch-3,t_birth(index_branch)+length(parent_branch_gen)-1),...
%                         X_truth(4*indices_branches(j)-3,t_birth(indices_branches(j)))];
%                     y_axis_link=[X_truth(4*index_branch-1,t_birth(index_branch)+length(parent_branch_gen)-1),...
%                               X_truth(4*indices_branches(j)-1,t_birth(indices_branches(j)))];

                    plot(x_axis_link,y_axis_link,'g','Linewidth',5)
                    
                   
                    break;
                end
                
                
                
               
                
                
                
            end
            
            
            
        end
        
        
    end
   
    
end






    
%We plot the measurements

plot(z(1,:),z(2,:),'oblack')