function filter_upd_pruned=TrPMBM_all_pruning(filter_upd, T_pruning,T_pruningPois,Nhyp_max,existence_threshold,T_alive)

%Author: Angel F. Garcia-Fernandez 
 
filter_upd_pruned=filter_upd;

%Pruning Poisson
index_remove=false(length(filter_upd.Pois),1);
for i=1:length(filter_upd.Pois)
    weightPois=filter_upd.Pois{i}.weightPois;
    if(weightPois<T_pruningPois)
        index_remove(i)=1;
    end
end
filter_upd_pruned.Pois(index_remove)=[];


%MB Pruning
weights=filter_upd.globHypWeight;
pruned=weights>T_pruning;
weights_pruned=weights(pruned);
globHyp=filter_upd.globHyp;
globHyp_pruned=globHyp(pruned,:);

%Pruning components so that there is a maximum number of components.
if length(weights_pruned)>Nhyp_max
    [~,pos]=sort(weights_pruned,'descend');
    index_pruned=pos(1:Nhyp_max);
else
    index_pruned=1:length(weights_pruned);
end
weights_pruned=weights_pruned(index_pruned)/sum(weights_pruned(index_pruned));
globHyp_pruned=globHyp_pruned(index_pruned,:);

%We remove Bernoulli components whose exitence probability is below a
%threshold, by setting the global hypotheses indices to zero (They will be
%deleted by the next algorithms.
Ntracks=length(filter_upd_pruned.tracks);
for i=1:Ntracks
    eB_i=filter_upd_pruned.tracks{i}.eB;
    remove=find(eB_i<existence_threshold);
    list_remove=ismember(globHyp_pruned(:,i),remove);  
    globHyp_pruned(list_remove,i)=0;
end

%We remove tracks(Bernoulli components) that do not take part in any global hypothesis
index=find(sum(globHyp_pruned,1)==0);
filter_upd_pruned.tracks(index)=[];
globHyp_pruned(:,index)=[];

%Pruning of tree_index
tree_index_pruned=filter_upd.tree_index;
tree_index_pruned(index)=[];




%We remove single-target hypotheses in each track (Bernoulli component) that do not belong to any global
%hypothesis
Ntracks=length(filter_upd_pruned.tracks);

for i=1:Ntracks
    Nhyp_i=length(filter_upd_pruned.tracks{i}.eB);
    indeces_valid=globHyp_pruned(:,i); %Track indices that appear in a global hypothesis
    %We remove single target hypothesis that do not appear in the global hypothesis
    index_remove=~ismember(1:Nhyp_i,indeces_valid);
    filter_upd_pruned.tracks{i}.meanB(index_remove)=[];
    filter_upd_pruned.tracks{i}.covB(index_remove)=[];
    filter_upd_pruned.tracks{i}.eB(index_remove)=[];
    filter_upd_pruned.tracks{i}.aHis(index_remove)=[];
    %filter_upd_pruned.tracks{i}.weightBLog(index_remove)=[];
    
    %Variables when we consider all trajectories
    filter_upd_pruned.tracks{i}.prob_length(index_remove)=[];
    filter_upd_pruned.tracks{i}.mean_past(index_remove)=[];
    
    
    %We need to change the indices of globHyp_pruned to account for the
    %removal of single target hypotheses
    
    if(sum(index_remove)>0)
        for j=1:size(globHyp_pruned,1)
            sub=sum(index_remove(1:globHyp_pruned(j,i)));
            globHyp_pruned(j,i)=globHyp_pruned(j,i)-sub;
        end
    end
end


% %We remove tracks (Bernoulli components) whose existence probability is very low in all its
% %single target hypotheses)
% Ntracks=length(filter_upd_pruned.tracks);
% index_remove=[];
% for i=1:Ntracks
%     eB_i=filter_upd_pruned.tracks{i}.eB;
%     remove=eB_i<existence_threshold;
%     if(sum(remove)==length(eB_i))
%         %We delete track (Bernoulli component)
%         index_remove=[index_remove,i];
%         
%     end
%     
% end
% filter_upd_pruned.tracks(index_remove)=[];
% %We adjust the global hypotheses
% globHyp_pruned(:,index_remove)=[];



%We also remove tracks (Bernoulli components) which are considered dead in
%all its single target hypotheses and its existence probaiblity is smaller
%than one. This means that the trajectory was never detected
Ntracks=length(filter_upd_pruned.tracks);
index_remove=[];
for i=1:Ntracks
    eB_i=filter_upd_pruned.tracks{i}.eB;
    prob_alive=zeros(length(eB_i),1);
    for j=1:length(prob_alive)
        prob_length_j=filter_upd_pruned.tracks{i}.prob_length{j};
        prob_alive(j)=prob_length_j(1);
    end
    dead_all=sum(prob_alive>T_alive)==0;
    if(dead_all)        
        remove=eB_i<1-10^(-5);
        if(sum(remove)==length(eB_i))
            %We delete track (Bernoulli component)
            index_remove=[index_remove,i];
        end
    end
    
end
filter_upd_pruned.tracks(index_remove)=[];
%We adjust the global hypotheses
globHyp_pruned(:,index_remove)=[];

%We adjust the tree_index variable
tree_index_pruned(index_remove)=[];






%When we eliminate a track (Bernoulli component), there can be
%duplicate global hypotheses that are merged into one by adding their
%weights and removing the duplicated ones.

if(~isempty(index_remove))
    [globHyp2,~,ic] = unique(globHyp_pruned,'rows');
    if(size(globHyp2,1)~=size(globHyp_pruned,1))
        %There are duplicate entries
        weights2=zeros(1,size(globHyp2,1));
        for i=1:size(globHyp2,1)
            indices_dupli=(ic==i);
            weights2(i)=sum(weights_pruned(indices_dupli));
        end
        globHyp_pruned=globHyp2;
        weights_pruned=weights2;
    end
end


%We store the result
filter_upd_pruned.globHyp=globHyp_pruned;
filter_upd_pruned.globHypWeight=weights_pruned/sum(weights_pruned);

filter_upd_pruned.tree_index=tree_index_pruned;

%Check everything is correct

if (length(tree_index_pruned)~=size(globHyp_pruned,2))
    display('error in pruning!!!!!!!')
end
