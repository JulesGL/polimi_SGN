% 30/10/2022 SGN - Jules GOMEL
% MatLab code for Assignment 1 - Exercise 2

clear;clc;close all
format long g
%% EX2 1) - Done
clear;clc;close all
format long g
% Constants
mu=1.21506683e-2;
Re=6378e3;
hi=167e3;
Rm=1738e3;
hf=100e3;
DU=3.84405000e8;
% Parameters
alpha=1.5*pi;
beta=1.41;
delta=7;
ti=0;
r0=(Re+hi)/DU;
ri=r0;
rf=(Rm+hf)/DU;
v0=beta*sqrt((1-mu)/r0);
% Initial state
x0=r0*cos(alpha)-mu;
y0=r0*sin(alpha);
x0dot=-(v0-r0)*sin(alpha);
y0dot=(v0-r0)*cos(alpha);
xx0=[x0 y0 x0dot y0dot];
% First guess solution
options=odeset('AbsTol',2.5e-14,'RelTol',2.5e-14);
[tt,xx]=ode113(@PBRFBP_rhs,[ti ti+delta],xx0, options, mu);

XX=zeros(size(xx));
for i=1:length(tt)
    X1=(xx(i,1)+mu)*cos(tt(i))-xx(i,2)*sin(tt(i));
    Y1=(xx(i,1)+mu)*sin(tt(i))+xx(i,2)*cos(tt(i));
    X1dot=(xx(i,3)-xx(i,2))*cos(tt(i))-(xx(i,4)-xx(i,1))*sin(tt(i));
    Y1dot=(xx(i,3)-xx(i,2))*sin(tt(i))+(xx(i,4)-xx(i,1))*cos(tt(i));
    XX(i,:)=[X1 Y1 X1dot Y1dot];
end
% Rotating frame
figure
hold on
xlabel('x')
ylabel('y')
plot(xx(:,1),xx(:,2))
scatter(-mu,0,10,'red')
scatter(1-mu,0,10,'k')
hold off
% Earth-centered
theta=0:0.01:2*pi;
circle=(1-mu)*[cos(theta);sin(theta)];
figure
hold on
xlabel('x')
ylabel('y')
plot(XX(:,1),XX(:,2))
plot(circle(1,:),circle(2,:))
axis equal
scatter(0,0,10,'red')
hold off

%% EX2 2)a) - Done

f=@(X) abs(deltavi(X))+abs(deltavf(X));
ti0=ti;
X0=[x0 y0 x0dot y0dot ti0 ti0+delta];
LB(1:4)=-Inf;
LB(5:6)=-0.1;
UB(1:4)=Inf;
UB(5)=2*pi/0.925195985;
UB(6)=UB(5)+23*4.34811305;
options = optimset('Display','iter', 'LargeScale', 'off','Algorithm','active-set');
Xsol=fmincon(@(X) f(X),X0,[],[],[],[],LB,UB,@(X) con(X),options);
f=f(Xsol)
% Plot the optimized solution
xx0_op=Xsol(1:4);
ti_op=Xsol(5);
tf_op=Xsol(6);


options=odeset('AbsTol',2.5e-14,'RelTol',2.5e-14);
[tt_op,xx_op]=ode113(@PBRFBP_rhs,[ti_op tf_op],xx0_op, options, mu);

XX_op=zeros(size(xx_op));
for i=1:length(tt_op)
    X1=(xx_op(i,1)+mu)*cos(tt_op(i))-xx_op(i,2)*sin(tt_op(i));
    Y1=(xx_op(i,1)+mu)*sin(tt_op(i))+xx_op(i,2)*cos(tt_op(i));
    X1dot=(xx_op(i,3)-xx_op(i,2))*cos(tt_op(i))-(xx_op(i,4)-xx_op(i,1))*sin(tt_op(i));
    Y1dot=(xx_op(i,3)-xx_op(i,2))*sin(tt_op(i))+(xx_op(i,4)-xx_op(i,1))*cos(tt_op(i));
    XX_op(i,:)=[X1 Y1 X1dot Y1dot];
end


% Rotating frame
figure
hold on
xlabel('x')
ylabel('y')
plot(xx_op(:,1),xx_op(:,2))
scatter(-mu,0,10,'red')
scatter(1-mu,0,10,'k')
hold off
% Earth-centered
figure
hold on
theta=0:0.01:2*pi;
circle=(1-mu)*[cos(theta);sin(theta)];
axis equal
plot(circle(1,:),circle(2,:))
xlabel('x')
ylabel('y')
plot(XX_op(:,1),XX_op(:,2))
scatter(0,0,10,'red')
hold off

%% EX2 2)b) - Done
clear;clc;close all
format long g
% First guess
constants=constant();
X0=[constants.xx0 constants.ti constants.ti+constants.delta];
LB(1:4)=-Inf;
LB(5:6)=0;
UB(1:4)=Inf;
UB(5)=2*pi/0.925195985;
UB(6)=UB(5)+23*4.34811305;

options_grad = optimoptions(@fmincon,'Display','iter',"Algorithm","active-set",...
   "EnableFeasibilityMode",true,'MaxFunctionEvaluations',10000,'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);

Xsol=fmincon(@(X) grad_f(X),X0,[],[],[],[],LB,UB,@(X) grad_con(X),options_grad);
f=grad_f(Xsol)
% Plot the optimized solution
mu=constants.mu;
xx0_op=Xsol(1:4);
ti_op=Xsol(5);
tf_op=Xsol(6);


options=odeset('AbsTol',2.5e-14,'RelTol',2.5e-14);
[tt_op,xx_op]=ode113(@PBRFBP_rhs,[ti_op tf_op],xx0_op, options, mu);

XX_op=zeros(size(xx_op));
for i=1:length(tt_op)
    X1=(xx_op(i,1)+mu)*cos(tt_op(i))-xx_op(i,2)*sin(tt_op(i));
    Y1=(xx_op(i,1)+mu)*sin(tt_op(i))+xx_op(i,2)*cos(tt_op(i));
    X1dot=(xx_op(i,3)-xx_op(i,2))*cos(tt_op(i))-(xx_op(i,4)-xx_op(i,1))*sin(tt_op(i));
    Y1dot=(xx_op(i,3)-xx_op(i,2))*sin(tt_op(i))+(xx_op(i,4)-xx_op(i,1))*cos(tt_op(i));
    XX_op(i,:)=[X1 Y1 X1dot Y1dot];
end

% Rotating frame
figure
hold on
xlabel('x')
ylabel('y')
plot(xx_op(:,1),xx_op(:,2))
scatter(-mu,0,10,'red')
scatter(1-mu,0,10,'k')
hold off
% Earth-centered
figure
hold on
axis equal 
xlabel('x')
ylabel('y')
theta=0:0.01:2*pi;
circle=(1-mu)*[cos(theta);sin(theta)];
axis equal
plot(circle(1,:),circle(2,:))
plot(XX_op(:,1),XX_op(:,2))
scatter(0,0,10,'red')
hold off

%% EX2 3) - Done
clear;clc;close all
format long g
% Initialize first-guess and N
mu=1.21506683e-2;
Re=6378e3;
hi=167e3;
Rm=1738e3;
hf=100e3;
DU=3.84405000e8;
alpha=1.5*pi;
beta=1.41;
delta=7;
ti=0;
r0=(Re+hi)/DU;
ri=r0;
rf=(Rm+hf)/DU;
v0=beta*sqrt((1-mu)/r0);
options=odeset('AbsTol',2.5e-14,'RelTol',2.5e-14);
% Initial state
x0=r0*cos(alpha)-mu;
y0=r0*sin(alpha);
x0dot=-(v0-r0)*sin(alpha);
y0dot=(v0-r0)*cos(alpha);
xx0=[x0 y0 x0dot y0dot];
N=4;
% Time tab
tt0=time_tab(ti,ti+delta,N);
% First guess tab
xx1=xx0;
[~,xx]=ode113(@PBRFBP_rhs,[tt0(1) tt0(end)],xx0, options, mu);


y0=[xx1 xx(round(length(xx1)/4),:) xx(round(length(xx1)/2),:) xx(end,:) tt0(1) tt0(end)];
UB=zeros(1,18);
LB=zeros(1,18);
UB(1:16)=Inf;
LB(1:16)=-Inf;
UB(17)=2*pi/0.925195985;
UB(18)=Inf;
% Optimization !
%options_grad = optimoptions(@fmincon,'Display','iter',"Algorithm","active-set",'MaxFunctionEvaluations',100000,'SpecifyConstraintGradient',true,'CheckGradient',true,'MaxIterations',50000);
options_grad = optimoptions(@fmincon,'Display','iter',"Algorithm","active-set",'MaxFunctionEvaluations',100000,'SpecifyObjectiveGradient',true,'MaxIterations',50000);
ysol=fmincon(@(X) obj_fun(X),y0,[],[],[],[],LB,UB,@(X) con_ms(X),options_grad);

xx1_sol=ysol(1:4);
xx2_sol=ysol(5:8);
xx3_sol=ysol(9:12);
xx4_sol=ysol(13:16);
tt_sol=time_tab(ysol(17),ysol(18),4);


[tt2,traj1]=ode113(@PBRFBP_rhs_STM,[tt_sol(1) tt_sol(2)],[xx1_sol reshape(eye(4),1,16)], options, mu);
[tt3,traj2]=ode113(@PBRFBP_rhs_STM,[tt_sol(2) tt_sol(3)],[xx2_sol reshape(eye(4),1,16)], options, mu);
[tt4,traj3]=ode113(@PBRFBP_rhs_STM,[tt_sol(3) tt_sol(4)],[xx3_sol reshape(eye(4),1,16)], options, mu);

figure
hold on
xlabel('x')
ylabel('y')
scatter(-mu,0,10,'red')
scatter(1-mu,0,10,'k')
%scatter(xx1_sol(1),xx1_sol(2),10,'b','filled')
scatter(xx2_sol(1),xx2_sol(2),10,'b','filled')
scatter(xx3_sol(1),xx3_sol(2),10,'b','filled')
%scatter(xx4_sol(1),xx4_sol(2),10,'b','filled')
plot(traj1(:,1),traj1(:,2))
plot(traj2(:,1),traj2(:,2))
plot(traj3(:,1),traj3(:,2))
hold off
%% Functions

function dXdt = PBRFBP_rhs(t,X,mu)
%RTBP_RHS Dynamics
    constants=constant();
    rho=constants.rho;
    ms=constants.ms;
    omegas=constants.omegas;

% State vector extraction

x=X(1);
y=X(2);
u=X(3);
v=X(4);


% RTBP dynamics

dXdt = [u ;
        v ;
        2*v+x - (mu*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(3/2)) - (ms*cos(omegas*t))/rho^2 + ((2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(3/2)) - (ms*(2*x - 2*rho*cos(omegas*t)))/(2*((x - rho*cos(omegas*t))^2 + (y - rho*sin(omegas*t))^2)^(3/2));
       -2*u+y - (ms*sin(omegas*t))/rho^2 - (mu*y)/((mu + x - 1)^2 + y^2)^(3/2) - (ms*(2*y - 2*rho*sin(omegas*t)))/(2*((x - rho*cos(omegas*t))^2 + (y - rho*sin(omegas*t))^2)^(3/2)) + (y*(mu - 1))/((mu + x)^2 + y^2)^(3/2)];

end

function constants=constant()
    % Constants
    constants.rho=3.88811143e2;
    constants.ms=3.28900541e5;
    constants.omegas=-9.25195985e-1;
    constants.N=4;
    constants.mu=1.21506683e-2;
    constants.Re=6378e3;
    constants.hi=167e3;
    constants.Rm=1738e3;
    constants.hf=100e3;
    constants.DU=3.84405000e8;
    constants.options=odeset('AbsTol',2.5e-14,'RelTol',2.5e-14);
    % Parameters
    constants.alpha=1.5*pi;
    constants.beta=1.41;
    constants.delta=7;
    constants.ti=0;
    constants.r0=(constants.Re+constants.hi)/constants.DU;
    constants.ri=constants.r0;
    constants.rf=(constants.Rm+constants.hf)/constants.DU;
    constants.v0=constants.beta*sqrt((1-constants.mu)/constants.r0);
    % Initial state
    constants.x0=constants.r0*cos(constants.alpha)-constants.mu;
    constants.y0=constants.r0*sin(constants.alpha);
    constants.x0dot=-(constants.v0-constants.r0)*sin(constants.alpha);
    constants.y0dot=(constants.v0-constants.r0)*cos(constants.alpha);
    constants.xx0=[constants.x0 constants.y0 constants.x0dot constants.y0dot];
end

%% Functions 2)a)
function Xf = phi(X)
    constants=constant();
    mu=constants.mu;
    xi=X(1:4);
    ti=X(5);
    tf=X(6);
    options=odeset('AbsTol',2.5e-14,'RelTol',2.5e-14);
    [~,xf]=ode113(@PBRFBP_rhs,[ti tf],xi, options, mu);
    Xf=[xf(end,:) ti tf];
end

function delta_vi=deltavi(X)
    constants=constant();
    mu=constants.mu;
    ri=constants.ri;
    xi=X(1:4);
    delta_vi=sqrt((xi(3)-xi(2))^2+(xi(4)+xi(1)+mu)^2)-sqrt((1-mu)/ri);
end

function delta_vf=deltavf(X)
    constants=constant();
    mu=constants.mu;
    rf=constants.rf;
    Xf=phi(X);
    xf=Xf(1:4);
    delta_vf=sqrt((xf(3)-xf(2))^2+(xf(4)+xf(1)+mu-1)^2)-sqrt((mu)/rf);
end

function [c, ceq]=con(X)
    constants=constant();
    mu=constants.mu;
    ri=constants.ri;
    rf=constants.rf;
    xi=X(1:4);
    ti=X(5);
    tf=X(6);
    Xf=phi(X);
    xf=Xf(1:4);
    c=ti-tf;
    ceq=[(xi(1)+mu)^2+xi(2)^2-ri^2 ;
        (xi(1)+mu)*(xi(3)-xi(2))+xi(2)*(xi(4)+xi(1)+mu) ; 
        (xf(1)+mu-1)^2+xf(2)^2-rf^2 ;
        (xf(1)+mu-1)*(xf(3)-xf(2))+xf(2)*(xf(4)+xf(1)+mu-1)];
end

%% Functions 2)b)

function RHS = PBRFBP_rhs_STM(t,X,mu)

    constants=constant();
    rho=constants.rho;
    ms=constants.ms;
    omegas=constants.omegas;
    % State vector extraction
    
    x=X(1);
    y=X(2);
    u=X(3);
    v=X(4);
    
    STM=(reshape(X(5:20),4,4))';
        
    % RTBP dynamics
    
    dXdt = [u ;
        v ;
        2*v+x - (mu*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(3/2)) - (ms*cos(omegas*t))/rho^2 + ((2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(3/2)) - (ms*(2*x - 2*rho*cos(omegas*t)))/(2*((x - rho*cos(omegas*t))^2 + (y - rho*sin(omegas*t))^2)^(3/2));
       -2*u+y - (ms*sin(omegas*t))/rho^2 - (mu*y)/((mu + x - 1)^2 + y^2)^(3/2) - (ms*(2*y - 2*rho*sin(omegas*t)))/(2*((x - rho*cos(omegas*t))^2 + (y - rho*sin(omegas*t))^2)^(3/2)) + (y*(mu - 1))/((mu + x)^2 + y^2)^(3/2)];
    
    % Variational approach 2D
    omegaxx=(mu - 1)/((mu + x)^2 + y^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) - ms/((x - rho*cos(omegas*t))^2 + (y - rho*sin(omegas*t))^2)^(3/2) + (3*ms*(2*x - 2*rho*cos(omegas*t))^2)/(4*((x - rho*cos(omegas*t))^2 + (y - rho*sin(omegas*t))^2)^(5/2)) + (3*mu*(2*mu + 2*x - 2)^2)/(4*((mu + x - 1)^2 + y^2)^(5/2)) - (3*(2*mu + 2*x)^2*(mu - 1))/(4*((mu + x)^2 + y^2)^(5/2)) + 1;
    omegaxy=(3*ms*(2*x - 2*rho*cos(omegas*t))*(2*y - 2*rho*sin(omegas*t)))/(4*((x - rho*cos(omegas*t))^2 + (y - rho*sin(omegas*t))^2)^(5/2)) + (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2));
    omegayy=(mu - 1)/((mu + x)^2 + y^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) - ms/((x - rho*cos(omegas*t))^2 + (y - rho*sin(omegas*t))^2)^(3/2) + (3*ms*(2*y - 2*rho*sin(omegas*t))^2)/(4*((x - rho*cos(omegas*t))^2 + (y - rho*sin(omegas*t))^2)^(5/2)) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2)^(5/2) + (3*mu*y^2)/((mu + x - 1)^2 + y^2)^(5/2) + 1;
    
    A=[0 0 1 0;
       0 0 0 1;
       omegaxx omegaxy  0 2;
       omegaxy omegayy -2 0];
    
    dSTMdt=A*STM;
    
    RHS=[dXdt;
        dSTMdt(1,1:4)';
        dSTMdt(2,1:4)';
        dSTMdt(3,1:4)';
        dSTMdt(4,1:4)'];

end

function [Xf,STM]=phi_STM(X)
    constants=constant();
    mu=constants.mu;
    xi=X(1:4);
    ti=X(5);
    tf=X(6);
    xi_STM=[xi reshape(eye(4),1,16)];
    options=odeset('AbsTol',2.5e-14,'RelTol',2.5e-14);
    [~,xf_STM]=ode113(@PBRFBP_rhs_STM,[ti tf],xi_STM, options, mu);
    xf=xf_STM(end,1:4);
    STM=reshape(xf_STM(end,5:20),4,4)';
    Xf=[xf(end,:) ti tf];
end

function [f,gradf]=grad_f(X)
    constants=constant();
    mu=constants.mu;
    ri=constants.ri;
    rf=constants.rf;
    [Xf,STM]=phi_STM(X);
    xf=Xf(1:4);
    xi=X(1:4);
    f=abs(sqrt((xi(3)-xi(2))^2+(xi(4)+xi(1)+mu)^2)-sqrt((1-mu)/ri)) ...
    + abs(sqrt((xf(3)-xf(2))^2+(xf(4)+xf(1)+mu-1)^2)-sqrt((mu)/rf));
    if nargout>1
        gradf=zeros(6,1);
        gradf(1:4)=1/sqrt((X(3)-X(2))^2+(X(4)+X(1)+mu)^2)*...
            [X(4)+X(1)+mu,...
            X(2)-X(3), ...
            -X(2)+X(3), ...
            X(4)+X(1)+mu];
        gradf(5)=1/sqrt((Xf(3)-Xf(2))^2+(Xf(4)+Xf(1)+mu-1)^2)*...
            [Xf(4)+Xf(1)+mu-1, ...
            Xf(2)-Xf(3), ...
            -Xf(2)+Xf(3), ...
            Xf(4)+Xf(1)+mu-1]*(-STM)*PBRFBP_rhs(X(5),X,mu);
        gradf(6)=1/sqrt((Xf(3)-Xf(2))^2+(Xf(4)+Xf(1)+mu-1)^2)*...
            [Xf(4)+Xf(1)+mu-1, ...
            Xf(2)-Xf(3), ...
            -Xf(2)+Xf(3), ...
            Xf(4)+Xf(1)+mu-1]*PBRFBP_rhs(X(6),Xf,mu);
    end
end

function [c, ceq,gc,jceq]=grad_con(X)
    constants=constant();
    mu=constants.mu;
    ri=constants.ri;
    rf=constants.rf;
    xi=X(1:4);
    ti=X(5);
    tf=X(6);
    [Xf,STM]=phi_STM(X);
    xf=Xf(1:4);
    c=ti-tf;
    ceq=[(xi(1)+mu)^2+xi(2)^2-ri^2 ;
        (xi(1)+mu)*(xi(3)-xi(2))+xi(2)*(xi(4)+xi(1)+mu) ; 
        (xf(1)+mu-1)^2+xf(2)^2-rf^2 ;
        (xf(1)+mu-1)*(xf(3)-xf(2))+xf(2)*(xf(4)+xf(1)+mu-1)];
    if nargout>2
        % Extraction of the useful blocks of the STM

        dXidt=PBRFBP_rhs(X(5),X,mu);
        dxfdt2=PBRFBP_rhs(X(6),Xf,mu);
        dxfdt1=-STM*dXidt;

        dc3dxi=[2*(xf(1)+mu-1)    2*xf(2)    0       0]*STM;
        dc3dti=[2*(xf(1)+mu-1)    2*xf(2)    0       0]*dxfdt1;
        dc3dt2=[2*(xf(1)+mu-1)    2*xf(2)    0       0]*dxfdt2;

        dc4dxi=[xf(3)    xf(4)      xf(1)+mu-1 xf(2)]*STM;
        dc4dti=[xf(3)    xf(4)      xf(1)+mu-1 xf(2)]*dxfdt1;
        dc4dt2=[xf(3)    xf(4)      xf(1)+mu-1 xf(2)]*dxfdt2;

        gc=[0 0 0 0 1 -1]';
        jceq=[2*(xi(1)+mu)    2*xi(2)    0       0    0 0;
              xi(3)    xi(4)      xi(1)+mu xi(2) 0 0;
              dc3dxi dc3dti dc3dt2;
              dc4dxi dc4dti dc4dt2]';
              
        
    end
end

%% Functions 3)

function tt=time_tab(ti,tf,N)
    tt=zeros(1,N);
    for j=1:N
        tt(j)=ti+(j-1)*(tf-ti)/(N-1);
    end
end

function [J,grad_J]=obj_fun(Y)
    constants=constant();
    mu=constants.mu;
    ri=constants.ri;
    rf=constants.rf;
    J=norm(sqrt((Y(3)-Y(2))^2+(Y(4)+Y(1)+mu)^2)-sqrt((1-mu)/ri))+...
    norm(sqrt((Y(3+12)-Y(2+12))^2+(Y(4+12)+Y(1+12)+mu-1)^2)-sqrt((mu)/rf));
    if nargout>1
        grad_J=zeros(1,18);
        P1=(1/sqrt((Y(3)-Y(2))^2+(Y(4)+Y(1)+mu)^2))*...
        [Y(4)+Y(1)+mu -Y(3)+Y(2) Y(3)-Y(2) Y(4)+Y(1)+mu];
        PN=(1/sqrt((Y(3+12)-Y(2+12))^2+(Y(4+12)+Y(1+12)+mu-1)^2))*...
        [Y(4+12)+Y(1+12)+mu-1 -Y(3+12)+Y(2+12) Y(3+12)-Y(2+12) Y(4+12)+Y(1+12)+mu-1];
        grad_J(1:4)=P1;
        grad_J(13:16)=PN;
        grad_J=grad_J';
    end
end

function [c,ceq,Jc,Jceq]=con_ms(Y)
    constants=constant();
    mu=constants.mu;
    ri=constants.ri;
    rf=constants.rf;
    Re=constants.Re;
    Rm=constants.Rm;
    DU=constants.DU;
    N=constants.N;
    options=constants.options;

    tt0=time_tab(Y(17),Y(18),4);
    
    eta1_e=(Re/DU)^2-(Y(1)+mu)^2-Y(2)^2;
    eta1_m=(Rm/DU)^2-(Y(1)+mu-1)^2-Y(2)^2;
    eta2_e=(Re/DU)^2-(Y(1+4)+mu)^2-Y(2+4)^2;
    eta2_m=(Rm/DU)^2-(Y(1+4)+mu-1)^2-Y(2+4)^2;
    eta3_e=(Re/DU)^2-(Y(1+8)+mu)^2-Y(2+8)^2;
    eta3_m=(Rm/DU)^2-(Y(1+8)+mu-1)^2-Y(2+8)^2;
    eta4_e=(Re/DU)^2-(Y(1+12)+mu)^2-Y(2+12)^2;
    eta4_m=(Rm/DU)^2-(Y(1+12)+mu-1)^2-Y(2+12)^2;
    tau=Y(17)-Y(18);
    c(1)=eta1_e;
    c(2)=eta1_m;
    c(3)=eta2_e;
    c(4)=eta2_m;
    c(5)=eta3_e;
    c(6)=eta3_m;
    c(7)=eta4_e;
    c(8)=eta4_m;
    c(9)=tau;

    psi1_1=(Y(1)+mu)^2+Y(2)^2-ri^2;
    psi1_2=(Y(1)+mu)*(Y(3)-Y(2))+Y(2)*(Y(4)+Y(1)+mu);
    psiN_1=(Y(3*4+1)+mu-1)^2+Y(3*4+2)^2-rf^2;
    psiN_2=(Y(3*4+1)+mu-1)*(Y(3*4+4)-Y(3*4+2))+Y(3*4+2)*(Y(3*4+4)+Y(3*4+1)+mu-1);

    xx1=Y(1:4);
    xx2=Y(5:8);
    xx3=Y(9:12);
    xx4=Y(13:16);
    STM0=reshape(eye(4),1,16);

    [~,xx2_STM4]=ode113(@PBRFBP_rhs_STM,[tt0(1) tt0(2)],[xx1 STM0], options, mu);
    [~,xx3_STM4]=ode113(@PBRFBP_rhs_STM,[tt0(2) tt0(3)],[xx2 STM0], options, mu);
    [~,xx4_STM4]=ode113(@PBRFBP_rhs_STM,[tt0(3) tt0(4)],[xx3 STM0], options, mu);

    xx2_flow=xx2_STM4(end,1:4);
    xx3_flow=xx3_STM4(end,1:4);
    xx4_flow=xx4_STM4(end,1:4);

    zeta_1=xx2-xx2_flow;
    zeta_2=xx3-xx3_flow;
    zeta_3=xx4-xx4_flow;
    ceq=[zeta_1 zeta_2 zeta_3 psi1_1 psi1_2 psiN_1 psiN_2 ];
    
    if nargout>2
        Jc=zeros(9,18);
        S1=[-2*(Y(1)+mu)   -2*Y(2) 0 0;...
            -2*(Y(1)+mu-1) -2*Y(2)  0 0];
        S2=[-2*(Y(1+4)+mu)   -2*Y(2+4) 0 0;...
            -2*(Y(1+4)+mu-1) -2*Y(2+4)  0 0];
        S3=[-2*(Y(1+8)+mu)   -2*Y(2+8) 0 0;...
            -2*(Y(1+8)+mu-1) -2*Y(2+8)  0 0];
        S4=[-2*(Y(1+12)+mu)   -2*Y(2+12) 0 0;...
            -2*(Y(1+12)+mu-1) -2*Y(2+12)  0 0];
        St=[1 -1];
        Jc(1:2,1:4)=S1;
        Jc(3:4,5:8)=S2;
        Jc(5:6,9:12)=S3;
        Jc(7:8,13:16)=S4;
        Jc(9,17:18)=St;
        Jc=Jc';

        Jceq=zeros(16,18);
        R1=[2*(Y(1)+mu) 2*Y(2) 0 0;
            Y(3) Y(4) Y(1)+mu Y(2)];
        RN=[2*(Y(1+12)+mu-1) 2*Y(2+12) 0 0;
            Y(4+12) Y(4+12) Y(1+12)+mu-1 Y(2+12)];

        phi2=reshape(xx2_STM4(end,5:20),4,4)';
        phi3=reshape(xx3_STM4(end,5:20),4,4)';
        phi4=reshape(xx4_STM4(end,5:20),4,4)';

        rhs1=PBRFBP_rhs(tt0(1),xx1,mu);
        rhs2=PBRFBP_rhs(tt0(2),xx2,mu);
        rhs3=PBRFBP_rhs(tt0(3),xx3,mu);

        rhs1_flow=PBRFBP_rhs(tt0(2),xx2_flow,mu);
        rhs2_flow=PBRFBP_rhs(tt0(3),xx3_flow,mu);
        rhs3_flow=PBRFBP_rhs(tt0(4),xx4_flow,mu);
        %Q1j, QNj
        Q11=-((N-1)/(N-1))*phi2*rhs1+((N-2)/(N-1))*rhs1_flow;
        Q1N=-((1-1)/(N-1))*phi2*rhs1+(1/(N-1))*rhs1_flow;
        Q21=-((N-2)/(N-1))*phi3*rhs2+((N-3)/(N-1))*rhs2_flow;
        Q2N=-((2-1)/(N-1))*phi3*rhs2+(2/(N-1))*rhs2_flow;
        Q31=-((N-3)/(N-1))*phi4*rhs3+((N-4)/(N-1))*rhs3_flow;
        Q3N=-((3-1)/(N-1))*phi4*rhs3+(3/(N-1))*rhs3_flow;

        Jceq(1:4,1:4)=phi2;
        Jceq(1:4,5:8)=-eye(4);
        Jceq(1:4,17)=Q11;
        Jceq(1:4,18)=Q1N;

        Jceq(5:8,5:8)=phi3;
        Jceq(5:8,9:12)=-eye(4);
        Jceq(5:8,17)=Q21;
        Jceq(5:8,18)=Q2N;

        Jceq(9:12,9:12)=phi4;
        Jceq(9:12,13:16)=-eye(4);
        Jceq(9:12,17)=Q31;
        Jceq(9:12,18)=Q3N;
        
        Jceq(13:14,1:4)=R1;
        Jceq(15:16,13:16)=RN;
        Jceq=Jceq';
    end
end