% Exercise 2 :  Batch filters
%
% MatLab code for Assignment 2 - Exercise 2 : state estimation and batch
% filters 
%
% SGN - 12/2022 - Jules GOMEL ☺
% AY 2022-23 -- Prof. F. Topputo and P. Di Lizia; TA: A. Morselli and M. Maestrini

%% Clear Constants & Kernels
clc; clear; close all
cspice_kclear();
format long g

%% Constants & Kernels
clc; clear; close all
cspice_kclear();
format long g

cspice_furnsh('assignment02.tm');
addpath('sgp4')

% Initialize data
t0='2022-11-11-19:08:49.824 UTC';
t0_vis_win='2022-11-12-04:30:00.000 UTC';
tf='2022-11-14-16:30:00.000 UTC';
t0_et=cspice_str2et(t0);
t0_vw_et=cspice_str2et(t0_vis_win);
tf_et=cspice_str2et(tf);

% Time span we need
% Visibility windows
t_span=t0_vw_et:60:tf_et;
% Integration
t_span_itg=t0_et:60:tf_et;

% Index offset to compute visibility window
t_span_0=t0_et:60:(t0_vw_et-60);
offset_t_span=length(t_span_0);

R_earth=cspice_bodvrd('EARTH','RADII',15);

% Satellite
r0_mean=[6054.30795817484,-3072.03883303992,-133.115352431876];
v0_mean=[4.64750094824087,9.18608475681236,-0.62056520749034];
x0_mean=[r0_mean,v0_mean]';

% KOUROU
coord_K=[5.25144,-52.80466,-14.67];
r_K=(R_earth'+[0 0 coord_K(3)]).*[cosd(coord_K(1))*cosd(coord_K(2)),cosd(coord_K(1))*sind(coord_K(2)),sin(coord_K(2))]+[0 0 coord_K(3)];
sigma_AzEl_K=100e-3;
sigma_r_K=0.01;
min_El_K=10;
R_K=diag([0.1;0.1;0.01].^2);
W_K=inv(R_K);
% PERTH
coord_P=[-31.80252,115.88516,22.16];
r_P=(R_earth'+[0 0 coord_P(3)]).*[cosd(coord_P(1))*cosd(coord_P(2)),cosd(coord_P(1))*sind(coord_P(2)),sin(coord_P(2))];
sigma_AzEl_P=100e-3;
sigma_r_P=0.01;
min_El_P=5;
R_P=R_K;
W_P=inv(R_P);
%% 1- Integrate trajectory

options=odeset('AbsTol',2.5e-14,'RelTol',2.5e-14);
[tt,xx]=ode113(@Kepler,t_span_itg,x0_mean,options);

Az_K=t_span;
Az_P=t_span;
El_K=t_span;
El_P=t_span;

for i=(offset_t_span+1):length(tt)
    ROT_ECI2TOPO_K = cspice_sxform('J2000', 'KOUROU_TOPO', tt(i));
    ROT_ECI2TOPO_P = cspice_sxform('J2000', 'PERTH_TOPO', tt(i));
    % Compute station position in ECI
    rv_K_eci = cspice_spkezr('KOUROU', tt(i), 'J2000', 'NONE', 'EARTH');
    rv_P_eci = cspice_spkezr('PERTH', tt(i), 'J2000', 'NONE', 'EARTH');
    % Compute station-satellite vector in ECI
    rv_K_sat_eci = xx(i,:)' - rv_K_eci;
    rv_P_sat_eci = xx(i,:)' - rv_P_eci;
    
    % Convert state into topocentric frame
    rv_K_sat_topo = ROT_ECI2TOPO_K*rv_K_sat_eci;
    rv_P_sat_topo = ROT_ECI2TOPO_P*rv_P_sat_eci;
    
    % Compute range, azimuth and elevation using cspice_xfmsta
    rll_K_sat = cspice_xfmsta(rv_K_sat_topo,'RECTANGULAR','LATITUDINAL','EARTH');
    rll_P_sat = cspice_xfmsta(rv_P_sat_topo,'RECTANGULAR','LATITUDINAL','EARTH');
    
    Az_K(i-offset_t_span) = rll_K_sat(2)*cspice_dpr();   % [rad]
    El_K(i-offset_t_span) = rll_K_sat(3)*cspice_dpr(); % [rad]
    
    Az_P(i-offset_t_span) = rll_P_sat(2)*cspice_dpr();   % [rad]
    El_P(i-offset_t_span) = rll_P_sat(3)*cspice_dpr(); % [rad]
end
%% 1- Plot & Compute visibility windows
figure
hold on
title('Az El - Visibility Windows KOUROU')
plot(t_span,Az_K)
plot(t_span,El_K)
plot(t_span,min_El_K*ones(size(t_span)),'red')
xlabel('Ephemeris time(s)')
ylabel('Angle (°)')
legend('Az','El')
hold off

figure
hold on
title('Az El - Visibility Windows PERTH')
plot(t_span,Az_P)
plot(t_span,El_P)
plot(t_span,min_El_P*ones(size(t_span)),'red')
xlabel('Ephemeris time(s)')
ylabel('Angle (°)')
legend('Az','El')
hold off

v_win_K = El_K>=min_El_K;
v_win_P = El_P>=min_El_P;

% Conjugate the two visibility windows
v_win = max(v_win_P,v_win_K);

v_win_K=double(v_win_K);
v_win_P=double(v_win_P);

v_win_K(not(v_win_K))=nan;
v_win_P(not(v_win_P))=nan;

%% Time of visibility windows
clc
tformat = 'YYYY-MON-DD-HR:MN:SC.####::UTC';

t_beg_K=[t_span(66),t_span(1771),t_span(2871)];
t_end_K=[t_span(556),t_span(2006),t_span(3545)];
t_beg_P=[t_span(588),t_span(2229)];
t_end_P=[t_span(1505),t_span(2837)];

for k=1:3
    k
    cspice_timout(t_beg_K(k),tformat)
    cspice_timout(t_end_K(k),tformat)
end

for k=1:2
        cspice_timout(t_beg_P(k),tformat)
    cspice_timout(t_end_P(k),tformat)
end

%% 2- a) Retrieve sat position from TLE
close all; clc

typerun    = 'u';  % user-provided inputs to SGP4 Matlab function
opsmode    = 'a';  % afspc approach ('air force space command')
whichconst =  72;  % WGS72 constants (radius, gravitational parameter)

%TLE from upper stage of Ariane 5
longstr1 = '1 87654U 22110B   22316.00967942  .00000002  00000-0  32024-3 0  9990';
longstr2 = '2 87654   3.6309 137.1541 8138191 196.1632  96.6141  1.26411866   834';

satrec = twoline2rv(longstr1, longstr2, typerun,'e', opsmode, whichconst);

% Get TLE epoch
[year,mon,day,hr,min,sec] = invjday(satrec.jdsatepoch, satrec.jdsatepochf);
tle_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
sat_epoch_et = cspice_str2et(tle_epoch_str);

fprintf('Satellite num ID: %d\n', satrec.satnum);
fprintf('TLE reference epoch: UTC %s\n', tle_epoch_str);

% Initialize & Propagate the trajectory using sgp4 :
m   = length(t_span) ;
ECI.x = zeros(6,m);
TEME.a = [0;0;0];

% daily corrections in rad 
arcsec2rad = pi / (180*3600);
[dPsi,dEps] = Read_EOP('EOP-Last5Years.txt',t_span);
ddpsi = dPsi * arcsec2rad;
ddeps = dEps * arcsec2rad;

%SGP4
for i = 1:m
    % Evaluate TLE at the corresponding epoch, with respect to satellite
    % epoch
    [satrec,TEME.x,TEME.v] = sgp4(satrec, (t_span(i) - sat_epoch_et)/60);
    
    % Centuries from TDT 2000 January 1 00:00:00.000
    ttt = cspice_unitim(t_span(i), 'ET', 'TDT')/cspice_jyear()/100;

    [ECI.x, ECI.v, ~] = teme2eci(TEME.x, TEME.v, TEME.a, ttt, ddpsi(i), ddeps(i));
    ECI.xx(:,i) = [ECI.x; ECI.v];
end

ECI_xx0 = ECI.xx(:,1) ;

%% 2- b) Simulate KOUROU measurements
clc
% SPK → station state in J2000
ECI.x_sta = cspice_spkezr('KOUROU', t_span, 'J2000', 'NONE', 'EARTH');

% Delta vector between sat and station
ECI.x_sta_exp = ECI.xx - ECI.x_sta;
 
% Rotation from J2000 to topo
ROT_ECI2TOPO = cspice_sxform('J2000','KOUROU_TOPO', t_span);
 
% Convert delta vector to topo
TOPO.x_sta_exp = zeros(size(ECI.x_sta_exp));
for j = 1:size(ECI.x_sta_exp,2)
    TOPO.x_sta_exp(:,j) = ROT_ECI2TOPO(:,:,j)*ECI.x_sta_exp(:,j);
end

% Perfect measurements :
rll_sat_K = cspice_xfmsta(TOPO.x_sta_exp,'RECTANGULAR','LATITUDINAL','EARTH');
measurements_perf_K.Az = rll_sat_K(2,:)*cspice_dpr().*v_win_K;
measurements_perf_K.El = rll_sat_K(3,:)*cspice_dpr().*v_win_K;
measurements_perf_K.Range = rll_sat_K(1,:).*v_win_K;

measurements_perf_K.All = [measurements_perf_K.Az;...
                                     measurements_perf_K.El;...
                                     measurements_perf_K.Range];

m=length(measurements_perf_K.All);
% Generate the noise :
noise_K = mvnrnd([0;0;0], R_K, m)';

% Add the noise :
measurements_noise_K.All = measurements_perf_K.All + noise_K;

% Recover each noisy measurements :
measurements_noise_K.Az    = measurements_noise_K.All(1,:);
measurements_noise_K.El    = measurements_noise_K.All(2,:);
measurements_noise_K.Range = measurements_noise_K.All(3,:);

% Plot example of measurement

figure
hold on
title('Elevation measurement with noise - KOUROU')
plot(t_span,measurements_noise_K.El,'red')
xlabel('Ephemeris time(s)')
ylabel('Elevation angles (°)')
hold off

%% 2- b) Simulate PERTH measurements

% SPK → station state in J2000
ECI.x_sta = cspice_spkezr('PERTH', t_span, 'J2000', 'NONE', 'EARTH');

% Delta vector between sat and station
ECI.x_sta_exp = ECI.xx - ECI.x_sta;
 
% Rotation from J2000 to topo
ROT_ECI2TOPO = cspice_sxform('J2000','PERTH_TOPO', t_span);
 
% Convert delta vector to topo
TOPO.x_sta_exp = zeros(size(ECI.x_sta_exp));
for j = 1:size(ECI.x_sta_exp,2)
    TOPO.x_sta_exp(:,j) = ROT_ECI2TOPO(:,:,j)*ECI.x_sta_exp(:,j);
end

% Perfect measurements :
rll_sat_P = cspice_xfmsta(TOPO.x_sta_exp,'RECTANGULAR','LATITUDINAL','EARTH');
measurements_perf_P.Az = rll_sat_P(2,:)*cspice_dpr().*v_win_P;
measurements_perf_P.El = rll_sat_P(3,:)*cspice_dpr().*v_win_P;
measurements_perf_P.Range = rll_sat_P(1,:).*v_win_P;

measurements_perf_P.All = [measurements_perf_P.Az;...
                                     measurements_perf_P.El;...
                                     measurements_perf_P.Range];

m=length(measurements_perf_P.All);
% Generate the noise :
noise_P = mvnrnd([0;0;0], R_P, m)';

% Add the noise :
measurements_noise_P.All = measurements_perf_P.All + noise_P;

% Recover each noisy measurements :
measurements_noise_P.Az    = measurements_noise_P.All(1,:);
measurements_noise_P.El    = measurements_noise_P.All(2,:);
measurements_noise_P.Range = measurements_noise_P.All(3,:);

%% 3- a) LS - PERTH only - Keplerian

% Encapsulate the variables in the function handle :
fun = @(x) obj_fun3a(x,t_span,options,v_win_P,measurements_noise_P,W_P) ;

% Optimization options :
opt = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Display','Iter');

% Execute least squares :
[x0_a, resnorm, residual,exitflag,~,~,jac] = lsqnonlin(fun,ECI_xx0,[],[],opt);

% Compute resulting covariance :
Jac_a = full(jac);
P_ls_a = resnorm / (length(residual) - length(ECI_xx0)) .* inv(Jac_a'*Jac_a);

sqrt(trace(P_ls_a(1:3,1:3)))
sqrt(trace(P_ls_a(4:6,4:6)))
%% 3- b) LS - PERTH and KOUROU - Keplerian 

% Encapsulate the variables in the function handle :
fun = @(x) obj_fun3b(x,t_span,options,v_win_P,measurements_noise_P,W_P, ...
    v_win_K,measurements_noise_K,W_K) ;

% Optimization options :
opt = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Display','Iter');

% Execute least squares :
[x0_b, resnorm, residual,exitflag,~,~,jac] = lsqnonlin(fun,ECI_xx0,[],[],opt);

% Compute resulting covariance :
Jac_b = full(jac);
P_ls_b = resnorm / (length(residual) - length(ECI_xx0)) .* inv(Jac_b'*Jac_b);

sqrt(trace(P_ls_b(1:3,1:3)))
sqrt(trace(P_ls_b(4:6,4:6)))
%% 3- c) LS - PERTH and KOUROU - Perturbed

% Encapsulate the variables in the function handle :
fun = @(x) obj_fun3c(x,t_span,options,v_win_P,measurements_noise_P,W_P, ...
    v_win_K,measurements_noise_K,W_K) ;

% Optimization options :
opt = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Display','Iter');

% Execute least squares :
[x0_c, resnorm, residual,exitflag,~,~,jac] = lsqnonlin(fun,ECI_xx0,[],[],opt);

% Compute resulting covariance :
Jac_c = full(jac);
P_ls_c = resnorm / (length(residual) - length(ECI_xx0)) .* inv(Jac_c'*Jac_c);

sqrt(trace(P_ls_c(1:3,1:3)))
sqrt(trace(P_ls_c(4:6,4:6)))
%% Functions

function RHS=Kepler(~,X)
    GM=cspice_bodvrd('EARTH','GM',1);
    r=norm(X(1:3));
    dXdt(1:3)=X(4:6);
    dXdt(4:6)=-GM/r^3*X(1:3);
    RHS=dXdt';
end

function [dPsi,dEps] = Read_EOP(finput,tspan)

% Generate a file ID
fid = fopen(finput);

% Skip the 34 first lines
for i = 1:34
    fgetl(fid);
end

% Generate the start time (first day at midnight)
t0_et = tspan(1);
t0_str = cspice_timout(t0_et,'YYYY-MM-DD-HR:MN:SC.####::UTC');
t_start = cspice_str2et( strcat(t0_str(1:10),'T00:00:00') );

% Get the first line and its date
line_first = fgetl(fid);
year = str2double(line_first(1:4));
mon  = str2double(line_first(6:7));
day  = str2double(line_first(9:10));

% Generate initial time 
t_init_et = cspice_str2et( sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,0,0,0]) );

% While intial time is different from the first day of measurements, skip
% the lines:
while t_start ~= t_init_et
    line = fgetl(fid);
    year = str2double(line(1:4));
    mon  = str2double(line(6:7));
    day  = str2double(line(9:10));
    t_init_et = cspice_str2et(sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,0,0,0]));
end

% Get the line of each measurement day:
line_19 = line;
line_20 = fgetl(fid);
line_21 = fgetl(fid);
line_22 = fgetl(fid);

% Skip the break lines between the 22nd and the 23rd of November
for i = 1:4
    fgetl(fid);
end

% Get the last line
line_23 = fgetl(fid);

% Initialize the vector of corrections
n = length(tspan);
dPsi = zeros(n,1);
dEps = zeros(n,1);

% Associate a correction to each time span value depending on their
% corresponding day
for i = 1:n
    t = tspan(i);
    if t < cspice_str2et('2021-11-20T00:00:00')
        dPsi(i) = str2double(line_19(60:68));
        dEps(i) = str2double(line_19(70:78));
    elseif t < cspice_str2et('2021-11-21T00:00:00')
        dPsi(i) = str2double(line_20(60:68));
        dEps(i) = str2double(line_20(70:78));
    elseif t < cspice_str2et('2021-11-22T00:00:00')
        dPsi(i) = str2double(line_21(60:68));
        dEps(i) = str2double(line_21(70:78));
    elseif t < cspice_str2et('2021-11-23T00:00:00')
        dPsi(i) = str2double(line_22(60:68));
        dEps(i) = str2double(line_22(70:78));
    else
        dPsi(i) = str2double(line_23(60:68));
        dEps(i) = str2double(line_23(70:78));
    end
end
end

function dxdt = Kepler_J2(t,xx)

GM=cspice_bodvrd('EARTH','GM',1);
RE=mean(cspice_bodvrd('EARTH','RADII',15));
% State :
x = xx(1) ;
y = xx(2) ;
z = xx(3) ;
u = xx(4) ;
v = xx(5) ;
w = xx(6) ;

r = norm(xx(1:3)) ;

pos_ECI = xx(1:3) ;
rotm = cspice_pxform('J2000', 'ITRF93', t);
pos_ECEF = rotm * pos_ECI ;
r_ECEF   = norm(pos_ECEF) ;

J2 = 0.00108263;
Coeff_J2 = (1.5*J2*GM*RE^2/r_ECEF^4);

a_J2_ECEF = [Coeff_J2*(pos_ECEF(1)/r_ECEF)*(5*pos_ECEF(3)^2/r_ECEF^2 - 1); 
             Coeff_J2*(pos_ECEF(2)/r_ECEF)*(5*pos_ECEF(3)^2/r_ECEF^2 - 1);
             Coeff_J2*(pos_ECEF(3)/r_ECEF)*(5*pos_ECEF(3)^2/r_ECEF^2 - 3)] ;

a_J2_ECI = rotm' * a_J2_ECEF ;

dxdt = [u;
        v;
        w;
        -GM/r^3*x +  a_J2_ECI(1)  ;
        -GM/r^3*y +  a_J2_ECI(2)  ;
        -GM/r^3*z +  a_J2_ECI(3) ];       
end

function residual = obj_fun3a(ECI_xx0,tspan,options,v_win_P,measurements_noise_P,W_P)

m = length(tspan) ;
residual = [];

% Propagation 
[~,ECI_xx] = ode113(@Kepler, tspan, ECI_xx0, options);

for i = 1:m
    % Initialisation :
    Diff_meas = [] ;
    Weights = [] ;

    if ~isnan(v_win_P(i)) %

            ECI_x_sta = cspice_spkezr('PERTH', tspan(i), 'J2000', 'NONE', 'EARTH') ;
            
            ECI_x_sta_exp = ECI_xx(i,:)' - ECI_x_sta ;
             
            ROT_ECI2TOPO = cspice_sxform('J2000','PERTH_TOPO', tspan(i));
    
            TOPO_x_sta_exp = ROT_ECI2TOPO * ECI_x_sta_exp ;
    
            % Recover the predicted measurements :
            rll_sat = cspice_xfmsta(TOPO_x_sta_exp,'RECTANGULAR','LATITUDINAL','EARTH');
            Range = rll_sat(1,:);
            Az    = rll_sat(2,:)*cspice_dpr();
            El    = rll_sat(3,:)*cspice_dpr();

            % Get difference vector :
            Diff_meas = [Diff_meas                               ;
                angdiff(Az, measurements_noise_P.Az(i)) ;
                El    - measurements_noise_P.El(i)      ;
                Range - measurements_noise_P.Range(i)  ];
    
            % Add the weight matrix
            Weights = blkdiag(Weights,W_P);
    end 
    % Compute the weighted :
    Diff_meas_weighted = Weights*Diff_meas;
    residual = [residual; Diff_meas_weighted];
end 
end 

function residual =obj_fun3b(ECI_xx0,tspan,options,v_win_P,measurements_noise_P,W_P, ...
    v_win_K,measurements_noise_K,W_K)

m = length(tspan) ;
residual = [];

% Propagation 
[~,ECI_xx] = ode113(@Kepler, tspan, ECI_xx0, options);

% PERTH
for i = 1:m
    % Initialisation :
    Diff_meas = [] ;
    Weights = [] ;

    if ~isnan(v_win_P(i)) %

            ECI_x_sta = cspice_spkezr('PERTH', tspan(i), 'J2000', 'NONE', 'EARTH') ;
            
            ECI_x_sta_exp = ECI_xx(i,:)' - ECI_x_sta ;
             
            ROT_ECI2TOPO = cspice_sxform('J2000','PERTH_TOPO', tspan(i));
    
            TOPO_x_sta_exp = ROT_ECI2TOPO * ECI_x_sta_exp ;
    
            % Recover the predicted measurements :
            rll_sat = cspice_xfmsta(TOPO_x_sta_exp,'RECTANGULAR','LATITUDINAL','EARTH');
            Range = rll_sat(1,:);
            Az    = rll_sat(2,:)*cspice_dpr();
            El    = rll_sat(3,:)*cspice_dpr();

            % Get difference vector :
            Diff_meas = [Diff_meas                               ;
                angdiff(Az, measurements_noise_P.Az(i)) ;
                El    - measurements_noise_P.El(i)      ;
                Range - measurements_noise_P.Range(i)  ];
    
            % Add the weight matrix
            Weights = blkdiag(Weights,W_P);
    end 
    % Compute the weighted :
    Diff_meas_weighted = Weights*Diff_meas;
    residual = [residual; Diff_meas_weighted];
end 

% KOUROU
for i = 1:m
    % Initialisation :
    Diff_meas = [] ;
    Weights = [] ;

    if ~isnan(v_win_K(i)) %

            ECI_x_sta = cspice_spkezr('KOUROU', tspan(i), 'J2000', 'NONE', 'EARTH') ;
            
            ECI_x_sta_exp = ECI_xx(i,:)' - ECI_x_sta ;
             
            ROT_ECI2TOPO = cspice_sxform('J2000','KOUROU_TOPO', tspan(i));
    
            TOPO_x_sta_exp = ROT_ECI2TOPO * ECI_x_sta_exp ;
    
            % Recover the predicted measurements :
            rll_sat = cspice_xfmsta(TOPO_x_sta_exp,'RECTANGULAR','LATITUDINAL','EARTH');
            Range = rll_sat(1,:);
            Az    = rll_sat(2,:)*cspice_dpr();
            El    = rll_sat(3,:)*cspice_dpr();

            % Get difference vector :
            Diff_meas = [Diff_meas                               ;
                angdiff(Az, measurements_noise_K.Az(i)) ;
                El    - measurements_noise_K.El(i)      ;
                Range - measurements_noise_K.Range(i)  ];
    
            % Add the weight matrix
            Weights = blkdiag(Weights,W_K);
    end 
    % Compute the weighted :
    Diff_meas_weighted = Weights*Diff_meas;
    residual = [residual; Diff_meas_weighted];

end
end

function residual =obj_fun3c(ECI_xx0,tspan,options,v_win_P,measurements_noise_P,W_P, ...
    v_win_K,measurements_noise_K,W_K)

m = length(tspan) ;
residual = [];

% Propagation 
[~,ECI_xx] = ode113(@Kepler_J2, tspan, ECI_xx0, options);

% PERTH
for i = 1:m
    % Initialisation :
    Diff_meas = [] ;
    Weights = [] ;

    if ~isnan(v_win_P(i)) %

            ECI_x_sta = cspice_spkezr('PERTH', tspan(i), 'J2000', 'NONE', 'EARTH') ;
            
            ECI_x_sta_exp = ECI_xx(i,:)' - ECI_x_sta ;
             
            ROT_ECI2TOPO = cspice_sxform('J2000','PERTH_TOPO', tspan(i));
    
            TOPO_x_sta_exp = ROT_ECI2TOPO * ECI_x_sta_exp ;
    
            % Recover the predicted measurements :
            rll_sat = cspice_xfmsta(TOPO_x_sta_exp,'RECTANGULAR','LATITUDINAL','EARTH');
            Range = rll_sat(1,:);
            Az    = rll_sat(2,:)*cspice_dpr();
            El    = rll_sat(3,:)*cspice_dpr();

            % Get difference vector :
            Diff_meas = [Diff_meas                               ;
                angdiff(Az, measurements_noise_P.Az(i)) ;
                El    - measurements_noise_P.El(i)      ;
                Range - measurements_noise_P.Range(i)  ];
    
            % Add the weight matrix
            Weights = blkdiag(Weights,W_P);
    end 
    % Compute the weighted :
    Diff_meas_weighted = Weights*Diff_meas;
    residual = [residual; Diff_meas_weighted];
end 

% KOUROU
for i = 1:m
    % Initialisation :
    Diff_meas = [] ;
    Weights = [] ;

    if ~isnan(v_win_K(i)) %

            ECI_x_sta = cspice_spkezr('KOUROU', tspan(i), 'J2000', 'NONE', 'EARTH') ;
            
            ECI_x_sta_exp = ECI_xx(i,:)' - ECI_x_sta ;
             
            ROT_ECI2TOPO = cspice_sxform('J2000','KOUROU_TOPO', tspan(i));
    
            TOPO_x_sta_exp = ROT_ECI2TOPO * ECI_x_sta_exp ;
    
            % Recover the predicted measurements :
            rll_sat = cspice_xfmsta(TOPO_x_sta_exp,'RECTANGULAR','LATITUDINAL','EARTH');
            Range = rll_sat(1,:);
            Az    = rll_sat(2,:)*cspice_dpr();
            El    = rll_sat(3,:)*cspice_dpr();

            % Get difference vector :
            Diff_meas = [Diff_meas                               ;
                angdiff(Az, measurements_noise_K.Az(i)) ;
                El    - measurements_noise_K.El(i)      ;
                Range - measurements_noise_K.Range(i)  ];
    
            % Add the weight matrix
            Weights = blkdiag(Weights,W_K);
    end 
    % Compute the weighted :
    Diff_meas_weighted = Weights*Diff_meas;
    residual = [residual; Diff_meas_weighted];

end
end