%% Model cablirate
clear
% Load data
load("empVolatilitySurfaceData.mat");%data
settings = calibrationSettings;

% Get values
kappa= settings.parameters0(1); 
theta= settings.parameters0(2); 
eta= settings.parameters0(3); 
rho= settings.parameters0(4); 
V0= settings.parameters0(5);

n = size(data.IVolSurf,1)*size(data.IVolSurf,2);

% Whether we want to calcualte standard error or not
settings.std = false;

param0 = [V0,theta,kappa, eta, rho];
fun = @(param)loss(data, settings, param);
%Optimization
[param_final, fFinal, exitFlag] = fminsearch(fun, param0, settings.calibrOptions);

%reference result

%For kappa, theta, eta, rho, and V0
%Estimated values: 10.4398    0.0235673      1.73195     -0.57519   0.00332752

%Standard errors: 1.0823  0.00034564     0.14861   0.0090438   0.0039123

%t -values: 9.6462      68.1855      11.6543     -63.6002      0.85053
%%In-sample MSE: 0.00011392

%% STD
% Calculate f value when plug in final parameters
f_value = loss(data, settings, param_final)*n;

% Indicate we want to calculate standard error
settings.std = true;
f = @(param) loss(data, settings, param);

p = settings.numberOfVariables;

% Jacobian matrix
J = jacobianest(f,param_final); % Jacobian

% sigma
sigma2 = f_value/(n - p);
sigma = sigma2*inv(J'*J);

% Standard errors
se = sqrt(diag(sigma))'; % SE with Jacobian

disp(['[V0, theta, kappa, eta, rho]']);
disp(['Starting values: ', num2str(param0)]);
disp(['Optimized values: ', num2str(param_final)]);
disp(['Standard error: ', num2str(se)]);
disp(['t-values: ', num2str(param_final./se)]);
disp(['In-sample MSE: ', num2str(f_value/n)]);
