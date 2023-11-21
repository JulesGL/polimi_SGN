clear
clc
close all
cspice_kclear();

%% build Camera Model:
foc=30; %20 mm
dens=54; %pix/mm
b=1;
p0=[1920/2;1200/2]; % center pixel location;
Cframe=[1,0,0;
    0,cos(-pi/2),sin(-pi/2);
    0,-sin(-pi/2),cos(-pi/2)];
R=10;

Cam.f=foc;
Cam.d=dens;
Cam.p0=p0;
Cam.b=b;
Cam.Cframe=Cframe;
Cam.R=R;

%% guess state:
cspice_furnsh('assignment02.tm');

% mean motion of GEO orbit computation:
mu=398600;
R=42241.08;
n=sqrt(mu/R^3);

%% propagate trajectories:

% initial time definition:
t0='1 April 2023 14:55:12.023 UTC';
epoch0=cspice_str2et(t0);

% random measurement time delta from t0;
tmeas = rand(1)*1000; 

%relative state after tmeas seconds from epoch0:
x0=12;
y0=-60;
z0=0;

% random attitude of the target after tmeas from epoch0:
q0=rand(4,1);
q0=q0/norm(q0); % quaternions are only rotations if unit norm!

% measurement extraction:

[meas]=meas_sim_pvt(n,[x0;y0;z0],q0,tmeas,epoch0,Cam);

%%
fprintf('\n Number of visible vertices: %d\n', length(meas.visible))
% 
figure(100)
scatter3(meas.y(1,:),meas.y(2,:),meas.y(3,:)); hold on; xlim([0,1920]);ylim([0,1200]);grid minor;
for i=1:length(meas.visible)
    text(meas.y(1,i),meas.y(2,i),meas.y(3,i),sprintf('%d',meas.visible(i)));
end
xlabel('$p_x$ [pix]','interpreter','latex'); ylabel('$p_y$ [pix]','interpreter','latex'); zlabel('d [pix]','interpreter','latex');
set(gca,'TickLabelInterpreter','latex');
