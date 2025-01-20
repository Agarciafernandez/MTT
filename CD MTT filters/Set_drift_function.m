function [f_drift,F_drift]=Set_drift_function(Gm0,beta0,Psi_aero,r0,h0)

%This function sets the drift function and the Jacobian for the re-entry
%tracking problem

%Author: Ángel F. García-Fernández


%Drift function
f_drift=@(x) f_drift_total(x,Gm0,beta0,Psi_aero,r0,h0);

%Gradient

F_drift=@(x) F_drift_total(x,Gm0,beta0,Psi_aero,r0,h0);



end

%Drift function
function output=f_drift_total(x,Gm0,beta0,Psi_aero,r0,h0)

distance=sqrt(x(1)^2+x(3)^2);
speed=sqrt(x(2)^2+x(4)^2);
g=-Gm0/distance^3;

drag=-beta0*exp(Psi_aero+(r0-distance)/h0)*speed;

%drag=0;

Matrix=[0 1 0 0;
    g drag 0 0;
    0 0 0 1
    0 0 g drag];
output=Matrix*x;

end

%Jacobian drift function
function output=F_drift_total(x,Gm0,beta0,Psi_aero,r0,h0)

distance=sqrt(x(1)^2+x(3)^2);
speed=sqrt(x(2)^2+x(4)^2);
exponential=exp(Psi_aero+(r0-distance)/h0);

%Derivative velocity w.r.t. position in x (and then y)
der_vel_pos_x=-Gm0/(distance^5)*(x(3)^2-2*x(1)^2)...
    +beta0*exponential*speed*x(1)*x(2)/(h0*distance);

der_vel_pos_y=-Gm0/(distance^5)*(x(1)^2-2*x(3)^2)...
    +beta0*exponential*speed*x(3)*x(4)/(h0*distance);





%Derivative velocity w.r.t. velocity in x (and then y)

der_vel_vel_x=-beta0*exponential*(x(2)^2/speed+speed);
der_vel_vel_y=-beta0*exponential*(x(4)^2/speed+speed);


output=[0 1 0 0;
    der_vel_pos_x der_vel_vel_x 0 0;
    0 0 0 1
    0 0 der_vel_pos_y der_vel_vel_y];

end