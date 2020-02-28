function DrawFilterEstimates(X_truth,t_birth,t_death,X_estimate,XLim,YLim,z,k)
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

for i=1:size(X_truth,1)/4
        plot(X_truth((i-1)*4+1,t_birth(i):t_death(i)-1),X_truth((i-1)*4+3,t_birth(i):t_death(i)-1),'b','LineWidth',1.2)
        text(X_truth((i-1)*4+1,t_birth(i)),X_truth((i-1)*4+3,t_birth(i)),num2str(i),'color','b')    
end

%We plot the state of the targets alive at the current time step
for i=1:size(X_truth,1)/4
    if(and(k>=t_birth(i),k<t_death(i)))
        plot(X_truth((i-1)*4+1,k),X_truth((i-1)*4+3,k),'ob','LineWidth',1.2)        
    end
end

X_binary=X_estimate~=0;

for k=1:size(X_binary,2)
    X_binario_k=X_binary(:,k);
    X_binario_m=reshape(X_binario_k,4,length(X_binario_k)/4);
    indices=find(sum(X_binario_m,1)>0);
    
    X_filtrado_k=X_estimate(:,k);
   
    for i=1:length(indices)     
        plot(X_filtrado_k(4*indices(i)-3),X_filtrado_k(4*indices(i)-1),'or')
    end

    
end
    
%We plot the measurements

plot(z(1,:),z(2,:),'oblack')