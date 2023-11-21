% Exercise 1 :  Uncertainty propagation
%
% MatLab code for Assignment 2 - Exercise 1 : Different methods of
% uncertainty propagation (LinCov, UT, MC...)
%
% SGN - 12/2022 - Jules GOMEL ☺
% AY 2022-23 -- Prof. F. Topputo and P. Di Lizia; TA: A. Morselli and M. Maestrini

%% Clear Constants & Kernels
clc; clear all; close all
cspice_kclear();
format long g

%% Constants & Kernels
cspice_furnsh('assignment02.tm');

t0 = '2022-11-11-19:08:49.824 UTC';
t0_et = cspice_str2et(t0);
r0_mean = [6054.30795817484, -3072.03883303992, -133.115352431876];
v0_mean = [4.64750094824087, 9.18608475681236, -0.62056520749034];
x0_mean = [r0_mean, v0_mean]';
P0 = [5.6e-3  3.5e-3 -7.1e-4 0      0      0;
      3.5e-3  9.7e-3  7.6e-4 0      0      0;
     -7.1e-4  7.6e-4  8.1e-4 0      0      0;
      0      0      0      2.8e-7  0      0;
      0      0      0      0      2.7e-7  0;
      0      0      0      0      0      9.6e-8];

GM = cspice_bodvrd('EARTH', 'GM', 5);
R_earth = cspice_bodvrd('EARTH', 'RADII', 5);

% Time to be studied
% options=odeset('AbsTol',2.5e-14,'RelTol',2.5e-14);
% [tt,xx]=ode113(@Kepler,[0 1000000],[x0_mean' reshape(eye(6),1,36)]',options);
% max=norm(xx(1,1:3));
% min=norm(xx(1,1:3));
% for i=1:length(tt)
%     if norm(xx(i,1:3))<min 
%         min=norm(xx(i,1:3));
%     else 
%         if norm(xx(i,1:3))>max
%             max=norm(xx(i,1:3));
%         end
%     end
% end
% a=.5*(min+max);
% T=sqrt(4*pi^2/GM*a^3);

T = 68346.6610065473;
t_tab = [0, .5, 1, 1.5, 2, 2.5, 3, 3.5, 4] * T + t0_et;

%% Ex1 - 1) LinCov
tic
X_lincov = x0_mean .* ones(6, 8);
P_lincov = P0 .* ones(6, 6, 8);

options = odeset('AbsTol', 2.5e-12, 'RelTol', 2.5e-12);
[tt, xx_STM] = ode113(@Kepler, t_tab, [x0_mean' reshape(eye(6), 1, 36)]', options);

for i = 1:8
    X_lincov(:, i) = xx_STM(i + 1, 1:6);
    STM_i = (reshape(xx_STM(i + 1, 7:42), 6, 6))';
    P_lincov(:, :, i) = STM_i * P0 * STM_i';
end
toc

%% Ex1 - 1) UT
tic
[X_UT, P_UT] = UT(x0_mean, P0, t_tab);
toc
%% Ex1 - 1) Comparison
[x, y, z] = sphere;
x = x * R_earth(1);
y = y * R_earth(1);
z = z * R_earth(3);

% Plotting standard deviation
figure 
hold on
xlabel('Time (odd → pericenter, even → apocenter)')
ylabel('Ranging standard deviation (km)')

for i = 1:8
    scatter(i, sqrt(trace(P_UT(1:3, 1:3, i))), 'blue', "square")
    scatter(i, sqrt(trace(P_lincov(1:3, 1:3, i))), 'red', "o")
end

legend('UT', 'LinCov')
hold off

% 3D plot
figure   
for i = 1:8
    plot3(X_UT(1, i), X_UT(2, i), X_UT(3, i), "square")
    hold on
    plot3(X_lincov(1, i), X_lincov(2, i), X_lincov(3, i), "o")
end

surf(x, y, z)
legend('UT1', 'LinCov1', 'UT2', 'LinCov2', 'UT3', 'LinCov3', 'UT4', 'LinCov4', ...
    'UT5', 'LinCov5', 'UT6', 'LinCov6', 'UT7', 'LinCov7', 'UT8', 'LinCov8')

xlabel('x-axis (km)')
ylabel('y-axis (km)')
zlabel('z-axis (km)')
axis equal
view(40, 35)

% Range plot
figure
hold on
for i = 1:8
    scatter(i, norm(X_lincov(1:3, i)), 'blue', "square")
    scatter(i, norm(X_UT(1:3, i)), 'red', "o")
end

legend('LinCov', 'UT')
xlabel('Time (odd → pericenter, even → apocenter)')
ylabel('Range (km)')
hold off

% Velocity plot
figure
hold on
for i = 1:8
    scatter(i, norm(X_lincov(4:6, i)), 'blue', "square")
    scatter(i, norm(X_UT(4:6, i)), 'red', "o")
end

legend('LinCov', 'UT')
xlabel('Time (odd → pericenter, even → apocenter)')
ylabel('Velocity (km/s)')
hold off

%% Ex1 - 2) Monte-Carlo
tic
X_MC = x0_mean .* ones(6, 8);
P_MC = P0 .* ones(6, 6, 8);

n = 100;
X_samples = mvnrnd(x0_mean', P0, n);
Y_samples = zeros(6, 100, 8);

for i = 1:n
    options = odeset('AbsTol', 2.5e-8, 'RelTol', 2.5e-8);
    [~, YY_STM] = ode113(@Kepler, t_tab, [X_samples(i, :) reshape(eye(6), 1, 36)]', options);

    for k = 1:8
        Y_samples(:, i, k) = YY_STM(k+1, 1:6);
    end
end

for k = 1:8
    X_MC(:, k) = 1/n * (sum(Y_samples(:, :, k), 2));
    P_MC(:, :, k) = 1/(n-1) * (Y_samples(:, :, k) - X_MC(:, k)) * ((Y_samples(:, :, k) - X_MC(:, k))');
end
toc

%% Ex1 - Covariance Matrix comparison
figure 
hold on
for i = 1:8
    scatter(i, sqrt(trace(P_UT(1:3, 1:3, i))), 'blue', 'square')
    scatter(i, sqrt(trace(P_lincov(1:3, 1:3, i))), 'red', 'o')
    scatter(i, sqrt(trace(P_MC(1:3, 1:3, i))), 'green', '*')
end

xlabel('Time (odd → pericenter, even → apocenter)')
ylabel('Position standard deviation (km)')
legend('UT', 'LinCov', 'MC')
hold off

figure 
hold on
for i = 1:8
    scatter(i, sqrt(trace(P_UT(4:6, 4:6, i))), 'blue', 'square')
    scatter(i, sqrt(trace(P_lincov(4:6, 4:6, i))), 'red', 'o')
    scatter(i, sqrt(trace(P_MC(4:6, 4:6, i))), 'green', '*')
end

xlabel('Time (odd → pericenter, even → apocenter)')
ylabel('Velocity standard deviation (km/s)')
legend('UT', 'LinCov', 'MC')
hold off

%% Ex1 - Subplot MC - Apogee 

% Apogee
figure 

subplot(2,2,1)
for i = 1:100
    scatter3(Y_samples(1, i, 1), Y_samples(2, i, 1), Y_samples(3, i, 1))
    hold on
end
error_ellipse(P_MC(1:3, 1:3, 1), X_MC(1:3, 1));
xlabel('X-coordinate in ECI [km]')
ylabel('Y-coordinate in ECI [km]')
zlabel('Z-coordinate in ECI [km]')
view(80, 60)
title('First apogee')

subplot(2,2,2)
for i = 1:100
    scatter3(Y_samples(1, i, 3), Y_samples(2, i, 3), Y_samples(3, i, 3))
    hold on
end
error_ellipse(P_MC(1:3, 1:3, 3), X_MC(1:3, 3));
xlabel('X-coordinate in ECI [km]')
ylabel('Y-coordinate in ECI [km]')
zlabel('Z-coordinate in ECI [km]')
view(80, 60)
title('Second apogee')

subplot(2,2,3)
for i = 1:100
    scatter3(Y_samples(1, i, 5), Y_samples(2, i, 5), Y_samples(3, i, 5))
    hold on
end
error_ellipse(P_MC(1:3, 1:3, 5), X_MC(1:3, 5));
xlabel('X-coordinate in ECI [km]')
ylabel('Y-coordinate in ECI [km]')
zlabel('Z-coordinate in ECI [km]')
view(80, 60)
title('Third apogee')

subplot(2,2,4)
for i = 1:100
    scatter3(Y_samples(1, i, 7), Y_samples(2, i, 7), Y_samples(3, i, 7))
    hold on
end
error_ellipse(P_MC(1:3, 1:3, 7), X_MC(1:3, 7));
xlabel('X-coordinate in ECI [km]')
ylabel('Y-coordinate in ECI [km]')
zlabel('Z-coordinate in ECI [km]')
view(80, 60)
title('Fourth apogee')

%% Ex1 - Subplot MC - Perigee
figure 

subplot(2,2,1)
for i = 1:100
    scatter3(Y_samples(1, i, 2), Y_samples(2, i, 2), Y_samples(3, i, 2))
    hold on
end
error_ellipse(P_MC(1:3, 1:3, 2), X_MC(1:3, 2));
xlabel('X-coordinate in ECI [km]')
ylabel('Y-coordinate in ECI [km]')
zlabel('Z-coordinate in ECI [km]')
view(80,60)
title('First perigee')

subplot(2,2,2)
for i = 1:100
    scatter3(Y_samples(1, i, 4), Y_samples(2, i, 4), Y_samples(3, i, 4))
    hold on
end
error_ellipse(P_MC(1:3, 1:3, 4), X_MC(1:3, 4));
xlabel('X-coordinate in ECI [km]')
ylabel('Y-coordinate in ECI [km]')
zlabel('Z-coordinate in ECI [km]')
view(80,60)
title('Second perigee')

subplot(2,2,3)
for i = 1:100
    scatter3(Y_samples(1, i, 6), Y_samples(2, i, 6), Y_samples(3, i, 6))
    hold on
end
error_ellipse(P_MC(1:3, 1:3, 6), X_MC(1:3, 6));
xlabel('X-coordinate in ECI [km]')
ylabel('Y-coordinate in ECI [km]')
zlabel('Z-coordinate in ECI [km]')
view(80,60)
title('Third perigee')

subplot(2,2,4)
for i = 1:100
    scatter3(Y_samples(1, i, 8), Y_samples(2, i, 8), Y_samples(3, i, 8))
    hold on
end
error_ellipse(P_MC(1:3, 1:3, 8), X_MC(1:3, 8));
xlabel('X-coordinate in ECI [km]')
ylabel('Y-coordinate in ECI [km]')
zlabel('Z-coordinate in ECI [km]')
view(80,60)
title('Fourth perigee')

%% Functions

function RHS = Kepler(~, X)
    GM = cspice_bodvrd('EARTH', 'GM', 1);
    r = norm(X(1:3));
    dXdt(1:3) = X(4:6);
    dXdt(4:6) = -GM/r^3 * X(1:3);
    STM = (reshape(X(7:42), 6, 6))';
    
    omegaxx = 3 * GM/r^5 * X(1)^2 - GM/r^3;
    omegaxy = 3 * GM/r^5 * X(2) * X(1);
    omegayy = 3 * GM/r^5 * X(2)^2 - GM/r^3;
    omegazx = 3 * GM/r^5 * X(3) * X(1);
    omegazy = 3 * GM/r^5 * X(3) * X(2);
    omegazz = 3 * GM/r^5 * X(3)^2 - GM/r^3;
    A = [0 0 0 1 0 0;
         0 0 0 0 1 0;
         0 0 0 0 0 1;
         omegaxx omegaxy omegazx  0 0 0;
         omegaxy omegayy omegazy  0 0 0;
         omegazx omegazy omegazz  0 0 0];
    dSTMdt = A * STM;

    RHS = [dXdt';
           dSTMdt(1,1:6)';
           dSTMdt(2,1:6)';
           dSTMdt(3,1:6)';
           dSTMdt(4,1:6)';
           dSTMdt(5,1:6)';
           dSTMdt(6,1:6)';];
end

% UT

function sigma_p = sigma_points(n, X0, P0)
    sigma_p = zeros(6, 2*n+1);
    sigma_p(:, 1) = X0;
    k = 0;
    alpha = 1e-3;
    lambda = alpha^2 * (n+k) - n;
    Pn = sqrtm((n+lambda) * P0);
    
    for i = 1:n
        sigma_p(:, i+1) = X0 + Pn(:, i);
        sigma_p(:, i+1+n) = X0 - Pn(:, i);
    end
end

function [X_UT, P_UT] = UT(x0_mean, P0, t_tab)
    n = 6;
    k = 0;
    beta = 2;
    alpha = 1e-3;
    lambda = alpha^2 * (n+k) - n;
    X_UT = zeros(6, 8);
    P_UT = zeros(6, 6, 8);
    x_p = sigma_points(n, x0_mean, P0);
    y_p = zeros(6, 13, k);
    
    for j = 1:size(x_p, 2)
        options = odeset('AbsTol', 2.5e-14, 'RelTol', 2.5e-14);
        [~, dX_p] = ode113(@Kepler, t_tab, [x_p(:, j)' reshape(eye(6), 1, 36)]', options);
        
        for k = 1:8
            y_p(:, j, k) = dX_p(k+1, 1:6)';
        end
    end
    
    W0m = lambda / (n+lambda);
    W0c = W0m + (1 - alpha^2 + beta);
    Wim = 1 / (2 * (n+lambda));
    Wic = Wim;
    
    for k = 1:8
        Y_mean_UT = y_p(:, 1, k) * W0m + sum(y_p(:, 2:end, k), 2) * Wim;
        P_mean_UT = zeros(6);
        
        for j = 1:(2*n+1)
            if j == 1
                P_mean_UT = P_mean_UT + W0c * (y_p(:, j, k) - Y_mean_UT) * (y_p(:, j, k) - Y_mean_UT)';
            else
                P_mean_UT = P_mean_UT + Wic * (y_p(:, j, k) - Y_mean_UT) * (y_p(:, j, k) - Y_mean_UT)';
            end
        end
        
        X_UT(:, k) = Y_mean_UT;
        P_UT(:, :, k) = P_mean_UT;
    end
end


%% File exchange - error ellispoid - cf References 
% AJ Johnson (2023). error_ellipse (https://www.mathworks.com/matlabcentral/fileexchange/4705-error_ellipse)
% MATLAB Central File Exchange. Retrieved November 21, 2023. 

function h=error_ellipse(varargin)
% ERROR_ELLIPSE - plot an error ellipse, or ellipsoid, defining confidence region
%    ERROR_ELLIPSE(C22) - Given a 2x2 covariance matrix, plot the
%    associated error ellipse, at the origin. It returns a graphics handle
%    of the ellipse that was drawn.
%
%    ERROR_ELLIPSE(C33) - Given a 3x3 covariance matrix, plot the
%    associated error ellipsoid, at the origin, as well as its projections
%    onto the three axes. Returns a vector of 4 graphics handles, for the
%    three ellipses (in the X-Y, Y-Z, and Z-X planes, respectively) and for
%    the ellipsoid.
%
%    ERROR_ELLIPSE(C,MU) - Plot the ellipse, or ellipsoid, centered at MU,
%    a vector whose length should match that of C (which is 2x2 or 3x3).
%
%    ERROR_ELLIPSE(...,'Property1',Value1,'Name2',Value2,...) sets the
%    values of specified properties, including:
%      'C' - Alternate method of specifying the covariance matrix
%      'mu' - Alternate method of specifying the ellipse (-oid) center
%      'conf' - A value betwen 0 and 1 specifying the confidence interval.
%        the default is 0.5 which is the 50% error ellipse.
%      'scale' - Allow the plot the be scaled to difference units.
%      'style' - A plotting style used to format ellipses.
%      'clip' - specifies a clipping radius. Portions of the ellipse, -oid,
%        outside the radius will not be shown.
%
%    NOTES: C must be positive definite for this function to work properly.
default_properties = struct(...
  'C', [], ... % The covaraince matrix (required)
  'mu', [], ... % Center of ellipse (optional)
  'conf', 0.5, ... % Percent confidence/100
  'scale', 1, ... % Scale factor, e.g. 1e-3 to plot m as km
  'style', '', ...  % Plot style
  'clip', inf); % Clipping radius
if length(varargin) >= 1 & isnumeric(varargin{1})
  default_properties.C = varargin{1};
  varargin(1) = [];
end
if length(varargin) >= 1 & isnumeric(varargin{1})
  default_properties.mu = varargin{1};
  varargin(1) = [];
end
if length(varargin) >= 1 & isnumeric(varargin{1})
  default_properties.conf = varargin{1};
  varargin(1) = [];
end
if length(varargin) >= 1 & isnumeric(varargin{1})
  default_properties.scale = varargin{1};
  varargin(1) = [];
end
if length(varargin) >= 1 & ~ischar(varargin{1})
  error('Invalid parameter/value pair arguments.') 
end
prop = getopt(default_properties, varargin{:});
C = prop.C;
if isempty(prop.mu)
  mu = zeros(length(C),1);
else
  mu = prop.mu;
end
conf = prop.conf;
scale = prop.scale;
style = prop.style;
if conf <= 0 | conf >= 1
  error('conf parameter must be in range 0 to 1, exclusive')
end
[r,c] = size(C);
if r ~= c | (r ~= 2 & r ~= 3)
  error(['Don''t know what to do with ',num2str(r),'x',num2str(c),' matrix'])
end
x0=mu(1);
y0=mu(2);
% Compute quantile for the desired percentile
k = sqrt(qchisq(conf,r)); % r is the number of dimensions (degrees of freedom)
hold_state = get(gca,'nextplot');
if r==3 & c==3
  z0=mu(3);
  
  % Make the matrix has positive eigenvalues - else it's not a valid covariance matrix!
  if any(eig(C) <=0)
    error('The covariance matrix must be positive definite (it has non-positive eigenvalues)')
  end
  % C is 3x3; extract the 2x2 matricies, and plot the associated error
  % ellipses. They are drawn in space, around the ellipsoid; it may be
  % preferable to draw them on the axes.
  Cxy = C(1:2,1:2);
  Cyz = C(2:3,2:3);
  Czx = C([3 1],[3 1]);
  [x,y,z] = getpoints(Cxy,prop.clip);
  h1=plot3(x0+k*x,y0+k*y,z0+k*z,prop.style);hold on
  [y,z,x] = getpoints(Cyz,prop.clip);
  h2=plot3(x0+k*x,y0+k*y,z0+k*z,prop.style);hold on
  [z,x,y] = getpoints(Czx,prop.clip);
  h3=plot3(x0+k*x,y0+k*y,z0+k*z,prop.style);hold on
  
  [eigvec,eigval] = eig(C);
  [X,Y,Z] = ellipsoid(0,0,0,1,1,1);
  XYZ = [X(:),Y(:),Z(:)]*sqrt(eigval)*eigvec';
  
  X(:) = scale*(k*XYZ(:,1)+x0);
  Y(:) = scale*(k*XYZ(:,2)+y0);
  Z(:) = scale*(k*XYZ(:,3)+z0);
  h4=surf(X,Y,Z);
  colormap gray
  alpha(0.3)
  camlight
  if nargout
    h=[h1 h2 h3 h4];
  end
elseif r==2 & c==2
  % Make the matrix has positive eigenvalues - else it's not a valid covariance matrix!
  if any(eig(C) <=0)
    error('The covariance matrix must be positive definite (it has non-positive eigenvalues)')
  end
  [x,y,z] = getpoints(C,prop.clip);
  h1=plot(scale*(x0+k*x),scale*(y0+k*y),prop.style);
  set(h1,'zdata',z+1)
  if nargout
    h=h1;
  end
else
  error('C (covaraince matrix) must be specified as a 2x2 or 3x3 matrix)')
end
%axis equal
set(gca,'nextplot',hold_state);
end
%---------------------------------------------------------------
% getpoints - Generate x and y points that define an ellipse, given a 2x2
%   covariance matrix, C. z, if requested, is all zeros with same shape as
%   x and y.
function [x,y,z] = getpoints(C,clipping_radius)
n=100; % Number of points around ellipse
p=0:pi/n:2*pi; % angles around a circle
[eigvec,eigval] = eig(C); % Compute eigen-stuff
xy = [cos(p'),sin(p')] * sqrt(eigval) * eigvec'; % Transformation
x = xy(:,1);
y = xy(:,2);
z = zeros(size(x));
% Clip data to a bounding radius
if nargin >= 2
  r = sqrt(sum(xy.^2,2)); % Euclidian distance (distance from center)
  x(r > clipping_radius) = nan;
  y(r > clipping_radius) = nan;
  z(r > clipping_radius) = nan;
end
end
%---------------------------------------------------------------
function x=qchisq(P,n)
% QCHISQ(P,N) - quantile of the chi-square distribution.
if nargin<2
  n=1;
end
s0 = P==0;
s1 = P==1;
s = P>0 & P<1;
x = 0.5*ones(size(P));
x(s0) = -inf;
x(s1) = inf;
x(~(s0|s1|s))=nan;
for ii=1:14
  dx = -(pchisq(x(s),n)-P(s))./dchisq(x(s),n);
  x(s) = x(s)+dx;
  if all(abs(dx) < 1e-6)
    break;
  end
end
end
%---------------------------------------------------------------
function F=pchisq(x,n)
% PCHISQ(X,N) - Probability function of the chi-square distribution.
if nargin<2
  n=1;
end
F=zeros(size(x));
if rem(n,2) == 0
  s = x>0;
  k = 0;
  for jj = 0:n/2-1;
    k = k + (x(s)/2).^jj/factorial(jj);
  end
  F(s) = 1-exp(-x(s)/2).*k;
else
  for ii=1:numel(x)
    if x(ii) > 0
      F(ii) = quadl(@dchisq,0,x(ii),1e-6,0,n);
    else
      F(ii) = 0;
    end
  end
end
end
%---------------------------------------------------------------
function f=dchisq(x,n)
% DCHISQ(X,N) - Density function of the chi-square distribution.
if nargin<2
  n=1;
end
f=zeros(size(x));
s = x>=0;
f(s) = x(s).^(n/2-1).*exp(-x(s)/2)./(2^(n/2)*gamma(n/2));
end
%---------------------------------------------------------------
function properties = getopt(properties,varargin)
%GETOPT - Process paired optional arguments as 'prop1',val1,'prop2',val2,...
%
%   getopt(properties,varargin) returns a modified properties structure,
%   given an initial properties structure, and a list of paired arguments.
%   Each argumnet pair should be of the form property_name,val where
%   property_name is the name of one of the field in properties, and val is
%   the value to be assigned to that structure field.
%
%   No validation of the values is performed.
%
% EXAMPLE:
%   properties = struct('zoom',1.0,'aspect',1.0,'gamma',1.0,'file',[],'bg',[]);
%   properties = getopt(properties,'aspect',0.76,'file','mydata.dat')
% would return:
%   properties = 
%         zoom: 1
%       aspect: 0.7600
%        gamma: 1
%         file: 'mydata.dat'
%           bg: []
%
% Typical usage in a function:
%   properties = getopt(properties,varargin{:})
% Process the properties (optional input arguments)
prop_names = fieldnames(properties);
TargetField = [];
for ii=1:length(varargin)
  arg = varargin{ii};
  if isempty(TargetField)
    if ~ischar(arg)
      error('Propery names must be character strings');
    end
    f = find(strcmp(prop_names, arg));
    if length(f) == 0
      error('%s ',['invalid property ''',arg,'''; must be one of:'],prop_names{:});
    end
    TargetField = arg;
  else
    % properties.(TargetField) = arg; % Ver 6.5 and later only
    properties = setfield(properties, TargetField, arg); % Ver 6.1 friendly
    TargetField = '';
  end
end
if ~isempty(TargetField)
  error('Property names and values must be specified in pairs.');
end
end
