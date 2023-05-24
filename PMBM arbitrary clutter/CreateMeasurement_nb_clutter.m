function z=CreateMeasurement_nb_clutter(X_multi_k,t_birth,t_death,p_d,Area,k,H,chol_R,Nx,r_nb,p_nb)

%Author: Angel F. Garcia-Fernandez

index_targets=and(t_birth<=k,t_death>k);

detected_targets=rand(1,length(index_targets))<p_d;

detected_targets=and(detected_targets,index_targets);

N_detected_targets=sum(detected_targets);

N_clutter_measurements=nbinrnd(r_nb,p_nb);

X_multi_k_r=reshape(X_multi_k,Nx,size(X_multi_k,2)*size(X_multi_k,1)/Nx);

z_targets=H*X_multi_k_r(:,detected_targets)+chol_R*randn(size(H,1),N_detected_targets);
z_clutter=[Area(1)*rand(1,N_clutter_measurements)...
    ;Area(2)*rand(1,N_clutter_measurements)];


z=[z_targets,z_clutter];

