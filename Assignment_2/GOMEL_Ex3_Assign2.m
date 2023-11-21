% Jules GOMEL 
% SGN - Assignement 2 - Exercise 3
% 23/12/2022 - Polimi

%% Clear Constants & Kernels
clc; clear; close all
cspice_kclear();
format long g

%% Constants & Kernels
clc; clear; close all
cspice_kclear();
format long g

cspice_furnsh('assignment02.tm');
addpath('simulator_ex3')

% Parameters of stereo camera (pix and meters)
cam.f=30;
cam.d=54;
cam.b=1;
cam.u0=960;
cam.v0=600;
cam.p0=[cam.u0;cam.v0];
Ss=[1920 1200];
sigma_u_square=10;
sigma_v_square=10;
sigma_d_square=10;
cam.R=diag([10 10 10]);

% LVLH frame to camera frame
cam.Cframe=[1 0 0;
         0 0 -1;
         0 1 0];

% Initial relative state at t0 in LVLH frame
t0='2023-04-01-14:55:12.023 UTC';
t0_et=cspice_str2et(t0);
tf_et=t0_et+24*60*60;
t_span=t0_et:tf_et;
r0=[15.792658268071492,-59.044939772661586,3.227106250277039];
v0=[-0.053960274403210,-0.053969644762889,-0.089140748762173];
x0=[r0,v0];
P0=diag([10 10 10 0.1 0.1 0.1]);
q0=[0.674156352338764,0.223585877389611,0.465489474399161,0.528055032413102]';
omega0=[-0.001262427155865,0.001204540074343,-0.000039180139156]';

% Parameters of SGN-I

l=10;
h=5;
d=3;
rho=1420;
m=l*h*d*rho;
J=m/12*diag([d^2+h^2,l^2+h^2,l^2+d^2]);
p=[l/2 -d/2 -h/2;
   l/2  d/2 -h/2;
   l/2  d/2  h/2;
   l/2 -d/2  h/2;
  -l/2 -d/2 -h/2;
  -l/2  d/2 -h/2;
  -l/2  d/2  h/2;
  -l/2 -d/2  h/2];


% mean motion of GEO orbit computation:
mu=398600;
R=42241.08;
n=sqrt(mu/R^3);

% Weights for UKF :

m=length(t_span);
n_dim   = 6 ;
alpha = 1e-3 ;
beta = 2 ;
k = 0 ;
lambda = alpha^2 * (n_dim + k) - n_dim ;

W0_m = lambda/(n_dim + lambda) ;
W0_c = W0_m + alpha^2 + 3    ;
Wi_m = 1/(2*(n_dim+lambda))    ;
Wi_c = Wi_m                  ;

W_M = [W0_m; Wi_m*ones(2*n_dim,1)];
W_C = [W0_c; Wi_c*ones(2*n_dim,1)];

%% 1- Integrate quat & state

r00=[12,-60,0];
v00=[1e-4,-2*n*r00(1),-1.2e-3];
options=odeset('AbsTol',2.5e-14,'RelTol',2.5e-14);
[tt_q,q_omega]=ode113(@euler,t_span,[q0;omega0],options,J);
[tt_x,xx]=ode113(@CW,t_span,[r00,v00],options,n);
qq=q_omega(:,1:4);

%% 1- Plot LVLH ref trajectory
figure
hold on 
xlabel('X-coordinate')
ylabel('Y-coordinate')
zlabel('Z-coordinate')
plot3(xx(:,1),xx(:,2),xx(:,3))
view(45,45)

%% 1- Simulate measurements 

for i=1:length(t_span)
    meas_array(i).meas=meas_sim(n,xx(i,1:3)',qq(i,:)',-t0_et+t_span(i),t0_et,cam);
end

%% 1- Number of visible vertices 
length_visibles=zeros(size(meas_array));
visible_v=[];

for i=1:length(t_span)
    meas=meas_array(i).meas;
    length_visibles(i)=length(meas.visible);
    for j=1:length(meas.visible)
        visible_v=[visible_v meas.y(:,j)];
    end
end

figure
hold on
xlabel('Ephemeris time (s)')
ylabel('Visible vertices')
grid on
plot(t_span,length_visibles)

%% 2- Plot visible vertices

figure
hold on
xlabel('Horizontal pixel')
ylabel('Vertical pixel')
plot(visible_v(1,:),visible_v(2,:),'.')
plot([0 1920 1920 0 0],[0 0 1200 1200 0],'k')
hold off

%% 3- UKF

% UKF struct will contain mean state and cov fields
UKF.x = zeros(n_dim,m) ;
UKF.P     = zeros(n_dim,n_dim,m);
UKF.x(:,1) = x0  ;
UKF.P(:,:,1)   = P0      ;

% Update step
for k =1:m-1    

    % Previous step
    x_k_1 = UKF.x(:,k) ;
    Pk_1= UKF.P(:,:,k)   ;
    t_k_1= t_span(k)   ;
    t_k= t_span(k+1) ;

    sqrt_Pk_1=sqrtm((n_dim+lambda)*Pk_1);

    % Frame matrices 
    C_TI=quat2dcm(qq(k,:));
    C_LI=dir_cos_mat_LVLH(n,k);

    % Generation of sigma points
    Sigma_k_1 = zeros(n_dim,13) ;
    Sigma_k_1(:,1) = x_k_1 ;

    for i = 2:(n_dim+1)
        Sigma_k_1(:,i)   = x_k_1 + sqrt_Pk_1(:,i-1);
    end 
    for i = 2:(n_dim+1)
        Sigma_k_1(:,i+n_dim) = x_k_1 - sqrt_Pk_1(:,i-1);
    end 

    % Use measurements array

    meas = meas_array(k).meas;
    y_k  = meas.y ;
    visible = meas.visible ;
    sz = size(y_k) ;
    width = sz(2) ;
    
    % Propagation of sigma points
    Sigma_k=zeros(n_dim,13) ;
    Gamma_k=zeros(3,width,13) ;
    Rk=10*eye(3*width) ;

    for i=1:(2*n_dim+1)
    
        [tt, xx_int] = ode113(@(t,x) CW(t,x,n), [t_k_1 t_k], Sigma_k_1(:,i), options); %integration
        Sigma_k(:,i) = xx_int(end,1:6) ;
        ind = 1 ;
        
        for j = visible 
            % Play with frames
            p_Body_target = p(j,1:3) ;
            p_Inertial_target = C_TI' * p_Body_target';
            p_LVLH_target = C_LI * p_Inertial_target;
            p_LVLH_target_chaser_centered = p_LVLH_target - Sigma_k(1:3,i) ; 
            p_Camera_frame_chaser_centered = cam.Cframe * p_LVLH_target_chaser_centered ;

            X = p_Camera_frame_chaser_centered(1) ;
            Y = p_Camera_frame_chaser_centered(2) ;
            Z = p_Camera_frame_chaser_centered(3) ;

            Gamma_k(1:3,ind,i) = [ cam.p0(1) - cam.d*cam.f*Y/Z ; cam.p0(2) + cam.d*cam.f*X/Z ; cam.b*cam.d*cam.f/Z] ;

            ind = ind + 1 ;

        end 
    end 

    % A priori state and cov

    x_k_minus = zeros(n_dim,1) ;
    y_k_minus = zeros(3,width) ;

    for i=1:(2*n_dim+1)
        x_k_minus = x_k_minus + W_M(i)*Sigma_k(:,i) ;
        y_k_minus = y_k_minus + W_M(i)*Gamma_k(:,:,i) ;
    end 

    P_k_minus = zeros(n_dim,n_dim) ;
    P_yy_k    = zeros(3*width,3*width) + Rk ;
    P_xy_k    = zeros(n_dim,3*width) ;

    for i=1:(2*n_dim+1)

        Dx = (Sigma_k(:,i) - x_k_minus) ;
        Dy = (Gamma_k(:,:,i) - y_k_minus) ;
        Dy = reshape(Dy,3*width,1) ;

        P_k_minus   = P_k_minus   + W_C(i)*(Dx*Dx') ;

        P_yy_k = P_yy_k + W_C(i)*(Dy*Dy') ;

        P_xy_k = P_xy_k + W_C(i)*(Dx*Dy') ;

    end 
    
    % Kalman gain

    K_k = P_xy_k/P_yy_k ;

    % A posteriori

    Dy2 = (y_k - y_k_minus)  ;
    Dy2 = reshape(Dy2,3*width,1) ;

    x_k_plus = x_k_minus + K_k * Dy2       ;
    P_k_plus     = P_k_minus - K_k * P_yy_k * K_k' ;

    UKF.x(:,k+1) = x_k_plus  ;
    UKF.P(:,:,k+1)   = P_k_plus      ;

end 

%% 3- EKF

n_dim=6;
EKF.x=zeros(n_dim,m);
EKF.P=zeros(n_dim,n_dim,m);
EKF.x(:,1)=x0;
EKF.P(:,:,1)=P0;

% From CW equations (h=1s)
STM_ = STM(n,1);

for k =1:m-1 

    x_k_1=EKF.x(:,k);

    Pk_1=EKF.P(:,:,k);

    t_k_1=t_span(k);
    t_k=t_span(k+1);

    [tt, xx_]=ode113(@(t,x) CW(t,x,n), [t_k_1 t_k], x_k_1, options); %integration
    x_k_minus=xx_(end,1:6) ;

    P_k_minus=STM_*Pk_1*STM_' ;

    meas=meas_array(k).meas;
    y_k=meas.y ;
    visible=meas.visible ;

    sz = size(y_k) ;
    width = sz(2) ;
    Rk = 10 * eye(3*width) ;

    y_k_minus = zeros(3,width);
    ind = 1; 

    H_k = zeros(3*width,6);

    % Frame matrices 

    C_TI = quat2dcm(qq(k,:));
    C_LI =dir_cos_mat_LVLH(n,k);

    for j = visible 

        p_Body_target = p(j,1:3) ;
        p_Inertial_target = C_TI' * p_Body_target';
        p_LVLH_target = C_LI * p_Inertial_target;
        p_LVLH_target_chaser_centered = p_LVLH_target - x_k_minus(1:3)';
        p_Camera_frame_chaser_centered = cam.Cframe * p_LVLH_target_chaser_centered ;
    
        X = p_Camera_frame_chaser_centered(1) ;
        Y = p_Camera_frame_chaser_centered(2) ;
        Z = p_Camera_frame_chaser_centered(3) ;
    
        y_k_minus(1:3,ind) = [ cam.p0(1) - cam.d*cam.f*Y/Z ; cam.p0(2) + cam.d*cam.f*X/Z ; cam.b*cam.d*cam.f/Z] ;

        % Jacobian with chain's rule as in the report

        JAC1 = [0                 -cam.d*cam.f/Z     cam.d*cam.f*Y/Z^2      ;
                cam.d*cam.f/Z      0                -cam.d*cam.f*X/Z^2      ;
                0                  0                -cam.d*cam.b*cam.f/Z^2 ];

        JAC2 = [ -cam.Cframe zeros(3,3) ];
        H_k(3*ind-2:3*ind,:) = JAC1*JAC2 ;
        ind = ind + 1 ;

    end 

    %Kalman gain 

    K_k = P_k_minus * H_k' / (H_k*P_k_minus*H_k' + Rk);

    %A posteriori

    Dy2 = (y_k - y_k_minus)  ;
    Dy2 = reshape(Dy2,3*width,1)     ;

    x_k_plus = x_k_minus' + K_k * Dy2 ;
    P_k_plus = (eye(6) -  K_k*H_k)*P_k_minus  ;

    EKF.x(:,k+1) = x_k_plus  ;
    EKF.P(:,:,k+1)   = P_k_plus      ;

end 

%% UKF and EKF plot

figure() ;
hold on 
plot3(xx(:,1), xx(:,2), xx(:,3),'r:','LineWidth',2)
plot3(EKF.x(1,1:m), EKF.x(2,1:m), EKF.x(3,1:m),'k','LineWidth',1.5)
legend('Reference','EKF')
view(140,10)
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
axis equal;

figure() ;
hold on 
plot3(xx(:,1), xx(:,2), xx(:,3),'r:','LineWidth',2)
plot3(UKF.x(1,:), UKF.x(2,:), UKF.x(3,:),'k','LineWidth',1.5)
legend('Reference','UKF')
view(140,10)
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
axis equal;




%% Functions

function dq_omegadt=euler(~,q_omega,J)
    q=q_omega(1:4);
    omega=q_omega(5:7);
    A=[0 -omega(1) -omega(2) -omega(3);
       omega(1) 0 omega(3) -omega(2);
       omega(2) -omega(3) 0 omega(1);
       omega(3) omega(2) -omega(1) 0];
    dq_omegadt(1:4)=.5*A*q;
    dq_omegadt(5:7)=-J\cross(omega,J*omega);
    dq_omegadt=dq_omegadt';
end


function dXdt=CW(~,X,n)
    dXdt=[X(4);
          X(5);
          X(6);
          3*n^2*X(1)+2*n*X(5);
          -2*n*X(4);
          -n^2*X(3)];
end

function Phi=STM(n,t)
    Phi=[4-3*cos(n*t)     0 0        1/n*sin(n*t)      2/n*(1-cos(n*t))       0;
         6*(sin(n*t)-n*t) 1 0        2/n*(cos(n*t)-1)  1/n*(4*sin(n*t)-3*n*t) 0;
         0                0 cos(n*t) 0                 0                      1/n*sin(n*t);
         3*n*sin(n*t)     0 0        cos(n*t)          2*sin(n*t)             0;
         6*n*(cos(n*t)-1) 0 0        -2*sin(n*t)       4*cos(n*t)-3           0;
         0                0 -n*sin(n*t) 0              0                      cos(n*t)];
end

function C_LI=dir_cos_mat_LVLH(n,t)
    C_LI=[cos(n*t) sin(n*t) 0;
         -sin(n*t) cos(n*t) 0;
          0        0        1];
end
