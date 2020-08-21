function DrawTrajectoryFilterEstimates(X_truth,t_birth,t_death,X_estimate,XLim,YLim,z,k)
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
    
        plot(X_truth((i-1)*4+1,t_birth(i):t_death(i)-1),X_truth((i-1)*4+3,t_birth(i):t_death(i)-1),'b','LineWidth',1)
        text(X_truth((i-1)*4+1,t_birth(i)),X_truth((i-1)*4+3,t_birth(i)),num2str(i),'color','b')    

        if(t_birth(i)<=k && k<=t_death(i)-1)
              plot(X_truth((i-1)*4+1,t_birth(i):k),X_truth((i-1)*4+3,t_birth(i):k),'b','LineWidth',2)
        end
        
    
end


for i=1:length(X_estimate)
    X_estimate_i=X_estimate{i};   
    plot(X_estimate_i(1:4:end),X_estimate_i(3:4:end),'-or')       
end
    
%We plot the measurements

plot(z(1,:),z(2,:),'oblack')