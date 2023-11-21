% 08/11/2022 SGN - Jules GOMEL
% MatLab code for Assignment 1 - Exercise 1

clear;clc;close all
format long g
%% Load Kernels

%% EX1 1)
mu=3.0359e-6;
dUdx = @(x) x - (1-mu)*(x+mu)./abs(x+mu).^3 - mu*(x+mu-1)./abs(x+mu-1).^3 ;

%L1=fzero(dUdx,[-mu+mu/10 1-mu-mu/10]);
%L3=fzero(dUdx,-10);
L2=fzero(dUdx,10);
%% EX1 2)

x0 = 1.008296144180133;
y0 = 0;
z0 = 0.001214294450297;
vx0 = 0;
vy0 = 0.010020975499502;
vz0 = 0;
xx0=[x0 y0 z0 vx0 vy0 vz0];
STM0=reshape(eye(6),1,36);
X0=[xx0 STM0]';

options=odeset('AbsTol',2.5e-14,'RelTol',2.5e-14,'Events',@ev_orbit);
[tt,xxSTM,tte,xxSTMe,ie]=ode113(@RTBP_rhs,[0 2*pi], X0, options, mu);

while(abs(xxSTMe(end,4))>1e-13 || abs(xxSTMe(end,6))>1e-13)
    dXdt=RTBP_rhs(1,xxSTMe,mu);
    phi = reshape(xxSTMe(end, 7:end), 6, 6)';

    XZdot_cross=[-xxSTMe(end,4) -xxSTMe(end,6)]';
    mat=[phi(4,1) phi(4,5);phi(6,1) phi(6,5)];
    mat2 = (1/dXdt(2)) * [dXdt(4);dXdt(6)]*[phi(2, 1) phi(2, 5)];

    delta_xydot=(mat-mat2)\XZdot_cross;

    X0(1)=X0(1)+delta_xydot(1);
    X0(5)=X0(5)+delta_xydot(2);

    [tt,xxSTM,tte,xxSTMe,ie]=ode113(@RTBP_rhs,[0 2*pi], X0, options, mu);
end

plot3(xxSTM(:,1),xxSTM(:,2),xxSTM(:,3))
%% EX1 3)
mu=3.0359e-6;
figure(1)
xlabel('x')
ylabel('y')
zlabel('z')

% Initialization 
x0 = 1.008296144180133;
y0 = 0;
z0 = 0.001214294450297;
vx0 = 0;
vy0 = 0.010020975499502;
vz0 = 0;
xx0=[x0 y0 z0 vx0 vy0 vz0];
STM0=reshape(eye(6),1,36);
X0=[xx0 STM0]';
for z=0.001214294450297:0.0001:0.0046
    hold on
    % keep the result from the previous iteration and change z
    X0(3)=z;
    % First guess orbit
    options=odeset('AbsTol',2.5e-14,'RelTol',2.5e-14,'Events',@ev_orbit);
    [tt,xxSTM,tte,xxSTMe,ie]=ode113(@RTBP_rhs,[0 2*pi], X0, options, mu);
    % Correction
    while(abs(xxSTMe(end,4))>1e-13 || abs(xxSTMe(end,6))>1e-13)
        dXdt=RTBP_rhs(1,xxSTMe,mu);
        phi = reshape(xxSTMe(end, 7:end), 6, 6)';
    
        XZdot_cross=[-xxSTMe(end,4) -xxSTMe(end,6)]';
        mat=[phi(4,1) phi(4,5);phi(6,1) phi(6,5)];
        mat2 = (1/dXdt(2)) * [dXdt(4);dXdt(6)]*[phi(2, 1) phi(2, 5)];
    
        delta_xydot=(mat-mat2)\XZdot_cross;
    
        X0(1)=X0(1)+delta_xydot(1);
        X0(5)=X0(5)+delta_xydot(2);
    
        [tt,xxSTM,tte,xxSTMe,ie]=ode113(@RTBP_rhs,[0 2*pi], X0, options, mu);
    end
    % Recontruction for plot
    plot3(xxSTM(:,1),xxSTM(:,2),xxSTM(:,3),'k',flip(xxSTM(:,1),2),-flip(xxSTM(:,2),2),flip(xxSTM(:,3),2),'k')
end
scatter3(L2,0,0,10,'red')
view(3)
%% Components of STM blocks
syms x y z mu;
r1=sqrt((x+mu)^2+y^2+z^2);
r2=sqrt((x+mu-1)^2+y^2+z^2);
U=0.5*(x^2+y^2)+(1-mu)/r1^3+mu/r2^3+0.5*mu*(1-mu);

Ux=diff(U,x);
Uy=diff(U,y);
Uz=diff(U,z);
omegazz=diff(Uz,z);
omegazy=diff(Uz,y);
omegazx=diff(Uz,x);

%% Functions

function RHS = RTBP_rhs(t,X,mu)
%RTBP_RHS Summary of this function goes here
%   Detailed explanation goes here

% State vector extraction

x=X(1);
y=X(2);
z=X(3);
u=X(4);
v=X(5);
w=X(6);

STM=(reshape(X(7:42),6,6))';

r1=((x+mu)^2+y^2+z^2)^(1/2);
r2=((x+mu-1)^2+y^2+z^2)^(1/2);

% RTBP dynamics

dXdt = [u ;
        v ;
        w;
        2*v+x-(1-mu)*(x+mu)/r1^3-mu*(x+mu-1)/r2^3;
       -2*u+y-(1-mu)*y/r1^3-mu*y/r2^3;
       -(1-mu)*z/r1^3-mu*z/r2^3];


% Variational approach 2D
omegaxx=1-(1-mu)/r1^3+3*(1-mu)*(x+mu)^2/r1^5-mu/r2^3+3*mu*(x-1+mu)^2/r2^5;
omegaxy=3*(1-mu)*(x+mu)*y/r1^5+3*mu*(x+mu-1)*y/r2^5;
omegayy=1-(1-mu)/r1^3+3*(1-mu)*y^2/r1^5-mu/r2^3+3*mu*y^2/r2^5;
omegazx=(3*mu*z*(2*mu + 2*x - 2))/(2*r2^5) - (3*z*(2*mu + 2*x)*(mu - 1))/(2*r1^5);
omegazy=(3*mu*y*z)/r2^5 - (3*y*z*(mu - 1))/r1^5;
omegazz=(mu-1)/r1^3-mu/r2^3-(3*z^2*(mu-1))/r1^5 + (3*mu*z^2)/r2^5;

A=[0 0 0 1 0 0;
   0 0 0 0 1 0;
   0 0 0 0 0 1;
   omegaxx omegaxy omegazx  0 2 0;
   omegaxy omegayy omegazy -2 0 0;
   omegazx omegazy omegazz  0 0 0];

dSTMdt=A*STM;

RHS=[dXdt;
    dSTMdt(1,1:6)';
    dSTMdt(2,1:6)';
    dSTMdt(3,1:6)';
    dSTMdt(4,1:6)';
    dSTMdt(5,1:6)';
    dSTMdt(6,1:6)';];

end

function [value,isterminal,direction] = ev_orbit(t,X,mu)
%EV_ORBIT2D event function for the x-axis crossing

value=X(2);
isterminal=1;
direction=-1; % stop when y decreases

end