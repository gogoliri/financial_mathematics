%%Option pricing & Calculate MSE
function res = loss(data, settings, param)
% Get the variable
V0= param(1); 
theta= param(2); 
kappa= param(3); 
eta= param(4); 
rho= param(5);
parameters = {V0; theta; kappa; eta; rho};

% Get values from data and settings
K = data.K;
T = data.T;
r = data.r;
S0 = data.S0;
IvolMarket  = data.IVolSurf;
n = settings.n;
model = settings.model;
N = size(data.IVolSurf,1)*size(data.IVolSurf,2);
std = settings.std;

% Get the constraint
minV0 = settings.minV0; maxV0 = settings.maxV0;
minTheta = settings.minTheta; maxTheta = settings.maxTheta;
minKappa = settings.minKappa; maxKappa = settings.maxKappa;
minEta = settings.minXi; maxEta = settings.maxXi;
minRho = settings.minRho; maxRho = settings.maxRho;

% Initialize Ct and IV
Ct = zeros(length(T), length(K));
IV = zeros(length(T), length(K));

% Calculate option price and implied volatility
for i = 1:length(T)
    Ct(i,:) = max(S0(1)*CallPricingFFT(model, n, 1 ,K'./S0, T(i), r,0, parameters{:}),0);%if FFT return negative value
    IV(i,:) = blsimpv(S0, K, r, T(i), Ct(i,:));

end
 
% Sum of square error
error = (IV - IvolMarket).^2;

% Check constraint
if V0 < minV0 || V0 > maxV0 ...
   || theta < minTheta || theta > maxTheta ...
   || kappa < minKappa || kappa > maxKappa ...
   || eta < minEta || eta > maxEta ...
   || rho < minRho || rho > maxRho
   
    inBounds = false;
else 
    inBounds = true;
end

% If calculate standard error, return IV, else return MSE
if std 
    res = IV(:);
else
    % Visualization
    drawnow    
    surf(IV)
    hold on
    surf(IvolMarket)
    hold off
    
    %Output
    if inBounds
        if sum(sum(isnan(error))) > 0
            error(isnan(error)) = 1;
        end
        res = sum(sum(error))/N;
    else
        % If not in bounds, give a *large* value
        res = 1E10;
    end

end

