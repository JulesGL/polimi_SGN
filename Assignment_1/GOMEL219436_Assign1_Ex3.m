% 30/10/2022 MSAS - Jules GOMEL
% MatLab code for Assignment 1

clear;clc;close all
cspice_kclear();
%% Load Kernels 
cspice_furnsh('kernels\naif0012.tls');
cspice_furnsh('kernels\pck00010.tpc');
cspice_furnsh('kernels\gm_de431.tpc');
cspice_furnsh('kernels\de425s.bsp');
lambda0_tf0=[51.548026115666690,48.626493624983230,-21.554596093166367,3.109557803551133e+02,-6.576266243768765e+02,4.120921760392176e+02,3.183284637742068,8.479878419506089e+03];
%% - Random Loop for generation of initial costates
clear;clc;close all
format long g
constants=constant();
exitflag=0;
while(exitflag<=0)
    int_pow_r1=round(-4+6*rand);
    int_pow_v1=round(-4+6*rand);
    int_pow_m=round(-4+6*rand);
    lambda0_tf0=[(-1+2*rand(1,3)).*10.^int_pow_r1  ...
        (-1+2*rand(1,3)).*10.^int_pow_v1 ...
        (rand).*10.^int_pow_m (constants.t0_et/constants.t_unit+0.5)];
    options_solve = optimoptions(@fsolve,'Display','iter');
    [lambda0_tf,fval,exitflag,output]=fsolve(@(X) f(X),lambda0_tf0,options_solve);
end

%% EX3 2) - Problem Solving
clear;clc;close all
format long g

constants=constant();
lambda0_tf0=[-4.849834917525271,6.814345119673250,-4.914356420569379,0.628569652137633,-0.512950062550021,0.858527246374456,0.349983765984809,23.102828170446930];
options_solve = optimoptions(@fsolve,'Display','iter','FunctionTolerance',0.0010/constants.l_unit,'MaxIterations',100000,'MaxFunctionEvaluations',100000);
[lambda0_tf,fval,exitflag,output]=fsolve(@(X) f(X),lambda0_tf0,options_solve);
%% Ex3 2) - Results display

% Time of ending computation
TOF=lambda0_tf(end)-constants.t0_et/constants.t_unit;
tformat='YYYY-MM-DD-HR:MN:SC.#### UTC';
t_end_et=lambda0_tf(end)*constants.t_unit;
t_end=cspice_timout(t_end_et,tformat);

% Integrating with the solution of the solver

x_lambda_0=[constants.x0 lambda0_tf(1:7)];
[~,xf_lambdaf]=ode113(@cs_dyn,[constants.t0_et/constants.t_unit lambda0_tf(8)],x_lambda_0,constants.options);
xf=xf_lambdaf(end,1:6);
xf_mars=cspice_spkezr('MARS',lambda0_tf(8)*constants.t_unit,'ECLIPJ2000','NONE','SUN');
xf_mars=[xf_mars(1:3)'/constants.l_unit xf_mars(4:6)'/constants.v_unit];

% Error in final state in km
delta_r=norm(xf(1:3)-xf_mars(1:3))*constants.l_unit
delta_v=norm(xf(4:6)-xf_mars(4:6))*1000*constants.v_unit

%% EX3 - 3)


%% Functions

function constants=constant()
    % Units
    constants.t_unit=60*60*24*365; % s/day
    constants.l_unit=149597870.700; % km/AU
    constants.v_unit=constants.l_unit/constants.t_unit; % km/s->AU/day
    constants.a_unit=149597870700/(constants.t_unit^2); % m/s²->AU/day2

    % Time constants
    
    constants.t0='2022-08-03-12:45:20.000 UTC';
    constants.t0_et=cspice_str2et(constants.t0);
    constants.Isp=3000/constants.t_unit;
    
    % Mass
    constants.m0=1500;
    
    % Length constants
    constants.mu=cspice_bodvrd('SUN','GM',15)*(constants.t_unit^2)/(constants.l_unit^3);
    constants.earth_x0=(cspice_spkezr('Earth',constants.t0_et,'ECLIPJ2000','NONE','SUN')');
    constants.x0=[constants.earth_x0(1:3)/constants.l_unit constants.earth_x0(4:6)/constants.v_unit constants.m0];

    % Force constants
    constants.g0=9.80665/constants.a_unit;
    constants.N=4;
    %constants.N=3;
    constants.Tmax_1=150e-3;
    constants.Tmax=constants.N*constants.Tmax_1/constants.a_unit;
    

    % Other
    constants.options=odeset('AbsTol',2.5e-14,'RelTol',2.5e-14);
    
    
    
end

function costate_dyn=cs_dyn(~,X_lambda)
    constants=constant();
    mu=constants.mu;
    Tmax=constants.Tmax;
    Isp=constants.Isp;
    g0=constants.g0;
    costate_dyn=zeros(size(X_lambda));
    % Equation of motion
    costate_dyn(1:3)=X_lambda(4:6);
    costate_dyn(4:6)=-(mu/(norm(X_lambda(1:3)))^3)*X_lambda(1:3)-(Tmax/(X_lambda(7)*norm(X_lambda(11:13))))*X_lambda(11:13);
    costate_dyn(7)=-Tmax/(g0*Isp*constants.m0);
    % Costate dynamics
    costate_dyn(8:10)=-(3*mu/(norm(X_lambda(1:3))^5))*(X_lambda(11:13)*X_lambda(1:3)')*X_lambda(1:3)+(mu/norm(X_lambda(1:3))^3)*X_lambda(11:13);
    costate_dyn(11:13)=-X_lambda(8:10);
    costate_dyn(14)=-(Tmax/X_lambda(7)^2)*norm(X_lambda(11:13));
end

function f=f(lambda0_tf)
    constants=constant();
    mu=constants.mu;
    x_lambda_0=[constants.x0 lambda0_tf(1:7)];
    [~,xf_lambdaf]=ode113(@cs_dyn,[constants.t0_et/constants.t_unit lambda0_tf(8)],x_lambda_0,constants.options);
    xf_mars=cspice_spkezr('MARS',lambda0_tf(8)*constants.t_unit,'ECLIPJ2000','NONE','SUN');
    xf_mars=[xf_mars(1:3)'/constants.l_unit xf_mars(4:6)'/constants.v_unit];
    f(1:6)=xf_lambdaf(end,1:6)-xf_mars;
    f(7)=xf_lambdaf(end,14);
    rhs_tf=cs_dyn(0,xf_lambdaf(end,:));
    f(8)=1+xf_lambdaf(end,8:13)*(rhs_tf(1:6)'-[xf_mars(4:6) -(mu/(norm(xf_mars(1:3))^3))*xf_mars(1:3)]');
end


