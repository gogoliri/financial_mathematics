%% Option pricing and stock simulation
clear
%Settings
S0 = 1;
V0 = 0.04;
K = 0.98;
H = 1.1;
T = 3;
r = 0.01;
kappaQ = 6;
thetaQ = 0.05;
eta = 0.5;
rho = -0.7;
dt = 1/252;
n = 10000;
L = T/dt;

%% MC simulation of stock price
% Initizalize
S = zeros(L, 1);
V = zeros(L, 1);
S_a = zeros(L, 1);
V_a = zeros(L, 1);

S(1) = S0;
V(1) = V0;
S_a(1) = S0;
V_a(1) = V0;

payoff = zeros(n, 1);
payoff_a = zeros(n, 1);

% Simulations
for i = 1:n
    % Initialize random factors
    epsilon1 = randn(L, 1);
    e12 = randn(L, 1);
    epsilon2 = rho*epsilon1 + sqrt(1-rho^2)*e12;
    
    % Antithetic version of noise
    epsilon1_a = -epsilon1;
    e12_a = -e12;
    epsilon2_a = rho*epsilon1_a + sqrt(1-rho^2)*e12_a;
    
    % Stock simulation
    for j = 2:L
        S(j) = S(j-1)*exp((r - 0.5*V(j-1))*dt + sqrt(V(j-1)*dt)*epsilon1(j-1));
        V(j) = max(0, V(j-1) + ...
        kappaQ*(thetaQ - V(j-1))*dt + ...
        eta*sqrt(V(j-1)*dt)*epsilon2(j-1) + ...
        0.25*(eta^2)*dt*(epsilon2(j-1)^2-1));
        
        % Antithetic version
        S_a(j) = S_a(j-1)*exp((r - 0.5*V_a(j-1))*dt + sqrt(V_a(j-1)*dt)*epsilon1_a(j-1));
        V_a(j) = max(0, V_a(j-1) + ...
        kappaQ*(thetaQ - V_a(j-1))*dt + ...
        eta*sqrt(V_a(j-1)*dt)*epsilon2_a(j-1) + ...
        0.25*(eta^2)*dt*(epsilon2_a(j-1)^2-1)); 
    end
    
    %Payout
    if any(S >= H)
        payoff(i) = max(0, S(end) - K);
    else 
        payoff(i) = 0;
    end
    
    if any(S_a >= H)
        payoff_a(i) = max(0, S_a(end) - K);
    else 
        payoff_a(i) = 0;
    end   
end

% Plot the last simulation of S and S_a
plot(S)
hold on
plot(S_a)
yline(K, "r", "Strike K")
yline(H, "b", "Barrier H")
hold off

% Call option price
C_MC = exp(-r*T)*mean((payoff + payoff_a)./2, 1);
%Implied Volatility
IV = blsimpv(S0, K, r, T, C_MC);

disp(['MC Heston price: ', num2str(C_MC)]);
disp(['Implied volatility : ', num2str(IV)]);
