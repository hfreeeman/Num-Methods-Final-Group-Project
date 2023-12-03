
% load "EnvironmentalForcing.mat"
load ('EnvironmentalForcing.mat');

% generate initial conditions and constant parameters
betaMaxVec = [0.5, 1, 2]; % max rate infection under ideal conditions (1/day)
mu_LMinVec = [5, 6, 7]; % min length of latent period (min number of days latent)
mu_I = 10; % rate infection clears (number of days infectious)
e = 0.001; % rate of introduction from external sources
Ap = 5000; % norm factor for population ('final' plant surface area in cm^2)
Pi = 1.33 * 30 * (-0.35968 + 0.10789 * 15 - 0.00214 * 15 * 15) * 30; % inital size of the population in cm^2 (equal to the model after 30 days w/ constant temp of 15C)
Si = Pi/Ap; % initial size of suseptible population (normalized)
Li = 0.01 * Si; % initial fraction of the population that's latent
Ii = 0; % initial fraction of population that's infectious
Ri = mu_I * Ii; % initial fraction of population that's recovered
Bi = 1; % initial size of berry population in cm^2 (berries grow starting 1st day of simulation)

time = linspace(0,61,1465);

Tbeta = zeros(size(T));
for i = 1:length(T)
   if T(i) < 35 && T(i) > 0
       Tbeta(i) = 0.0002417.*(T(i).^2.08737) .* (35 - T(i)).^0.72859;
   else
       Tbeta(i) = 0;
   end
end



% call SLIRP function as function handler 
func = @(S, L, I, R, P, B, params, tSpan, T) SLIRP(S, L, I, R, P, B, params, tSpan, T );

%slopes = SLIRP(1, initial_conditions, beta(1), mu_L(1), T(1), tspan(1), mu_I, e, Ap )
% call rk4 function
% function takes output of SLIRP, tspan, and inital conditions for the 6
% ODEs as inputs.
% returns value of 6 dependent variables.

% plot the dependent variables versus time
% plot total population P and Berry Population B after normalizing by A
% (P/A and B/A)
figure
plot(tspan, T)
grid on
title("forcing air temperature as a function of days")
xlim([0, 61])
xlabel("time [days]")
ylabel("temperature [c]")

% for changing betaMax value

for i = 1:length(betaMaxVec)

% initialize parameters 
initial_cond = [Si; Li; Ii; Ri; Pi; Bi];
params = [betaMaxVec(i); mu_LMinVec(2); mu_I; e; Ap];

% call rk4 for each value of betaMax 
[S, I, L, R, P, B] = RK4(func, tspan, initial_cond, params, T);


figure
plot(tspan, P/Ap, tspan, B/Ap, '--', tspan, S, '-.', tspan, L, tspan, I, '--',  tspan, R * .5, '-.', 'LineWidth', 1.5 )
legend("Total Population", "Berry Population", "Susceptible", "Latent", "Infected", "Removed", Location="best")
ylim([0, 1.6])
xlim([0, 61])
grid on
ylabel("Population (fraction of initial)")
xlabel("time [days]")
txt = sprintf("SLIRP Model System for Beta_m_a_x = %.1f & Mu_Lmin = %.1f", betaMaxVec(i), mu_LMinVec(2));
title(txt)

end

% for changing mu_lmin values 
for i = 1:length(betaMaxVec)

% initialize parameters 
initial_cond = [Si; Li; Ii; Ri; Pi; Bi];
params = [betaMaxVec(2); mu_LMinVec(i); mu_I; e; Ap];

% call rk4 for each value of betaMax 
[S, I, L, R, P, B] = RK4(func, tspan, initial_cond, params, T);


figure
plot(tspan, P/Ap, tspan, B/Ap, '--', tspan, S, '-.', tspan, L, tspan, I, '--',  tspan, R * .5, '-.', 'LineWidth', 1.5 )
legend("Total Population", "Berry Population", "Susceptible", "Latent", "Infected", "Removed", Location="best")
ylim([0, 1.6])
xlim([0, 61])
grid on
ylabel("Population (fraction of initial)")
xlabel("time [days]")
txt = sprintf("SLIRP Model System for Beta_m_a_x = %.1f & Mu_Lmin = %.1f", betaMaxVec(2), mu_LMinVec(i));
title(txt)

end


%% functions 

function [slopes] = SLIRP(S, L, I, R, P, B, params, tSpan, T )
    
% unpack inital varibles

% unpack and solve for initial parameters
    beta = params(1);
    mu_L = params(2);
    mu_I = params(3);
    e = params(4);
    A_p = params(5);
    P_b = B;

    % find needed values 
T_B = zeros(size(T));
for i = 1:size(T)
   if T(i) < 35 && T(i) > 0
       T_B(i) = 0.0002417.*(T(i)^2.08737) .* (35 - T(i))^0.72859;
   else
       T_B(i) = 0;
   end
end


    t_day = 30 + tSpan;
    T_E = -0.35968 + 0.10789.*T - 0.00214.*T.^2;
    dP_l_dt = (1.33 .* t_day) .* T_E;
    dP_b_dt = (0.1724 .* P_b - 0.0000212 .* (P_b .^2)) .* T_E;

    % Calculate derivatives
    dL_dt = beta * S * I - (mu_L^-1) * L + e;
    dI_dt = (mu_L^-1) * L - (mu_I^-1) * I;
    dR_dt = (mu_L^-1) * I;
    dP_dt = dP_b_dt + dP_l_dt;
    dS_dt = -beta * S * I + (dP_dt) / A_p;
    

    % Return slopes
    slopes = [dS_dt; dL_dt; dI_dt; dR_dt; dP_dt; dP_b_dt];

end

function [S, I, L, R, P, B] = RK4(func, tSpan, initial_cond, params, temp)

% unpack inital varibles
    S(1) = initial_cond(1);
    I(1) = initial_cond(2);
    L(1) = initial_cond(3);
    R(1) = initial_cond(4);
    P(1) = initial_cond(5);
    B(1) = initial_cond(6);

% unpack and solve for initial parameters
    betaMax = params(1);
    mu_LMin = params(2);
    mu_I = params(3);
    e = params(4);
    Ap = params(5);
    t1 = tSpan(1);
    t2 = tSpan(2);
    h = t2 - t1;

% Tbeta using eq. 7
Tbeta = zeros(size(temp));
for i = 1:length(temp)
   if temp(i) < 35 && temp(i) > 0
       Tbeta(i) = 0.0002417.*(temp(i).^2.08737) .* (35 - temp(i)).^0.72859;
   else
       Tbeta(i) = 0;
   end
end

% beta using eq. 8
beta = betaMax * Tbeta;

% mu_L using eq. 9
mu_L = zeros(size(tSpan));
j = 1;
for i = 1:length(tSpan) 
    mu_L(i) = sum(Tbeta(j:i));
    while (mu_L(i) > mu_LMin) 
        j = j + 1;
        mu_L(i) = sum(Tbeta(j:i));
    end
end



for i = 2:length(tSpan)

    params1 = [beta(i), mu_L(i), mu_I, e, Ap, tSpan(i)];
    params2 = [beta(i), mu_L(i), mu_I, e, Ap, tSpan(i) + 0.5 * h];

 % First pass
    [slopes1] = func(S(i-1), L(i-1), I(i-1), R(i-1), P(i-1), B(i-1), params1(1:5), params1(6), temp(i));
    S1 = S(i-1) + 0.5 * slopes1(1) * h;
    L1 = L(i-1) + 0.5 * slopes1(2) * h;
    I1 = I(i-1) + 0.5 * slopes1(3) * h;
    R1 = R(i-1) + 0.5 * slopes1(4) * h;
    P1 = P(i-1) + 0.5 * slopes1(5) * h;
    B1 = B(i-1) + 0.5 * slopes1(6) * h;

 % Second pass
    [slopes2] = func(S1, L1, I1, R1, P1, B1, params2(1:5), params2(6), temp(i));
    S2 = S(i-1) + 0.5 * slopes2(1) * h;
    L2 = L(i-1) + 0.5 * slopes2(2) * h;
    I2 = I(i-1) + 0.5 * slopes2(3) * h;
    R2 = R(i-1) + 0.5 * slopes2(4) * h;
    P2 = P(i-1) + 0.5 * slopes2(5) * h;
    B2 = B(i-1) + 0.5 * slopes2(6) * h;

 % third pass
    [slopes3] = func(S2, L2, I2, R2, P2, B2, params2(1:5), params2(6), temp(i));
    S3 = S(i-1) + 0.5 * slopes3(1) * h;
    L3 = L(i-1) + 0.5 * slopes3(2) * h;
    I3 = I(i-1) + 0.5 * slopes3(3) * h;
    R3 = R(i-1) + 0.5 * slopes3(4) * h;
    P3 = P(i-1) + 0.5 * slopes3(5) * h;
    B3 = B(i-1) + 0.5 * slopes3(6) * h;

% fourth pass 
    [slopes4] = func(S3, L3, I3, R3, P3, B3, params2(1:5), params2(6), temp(i) );

% update SLIRP
    S(i) = S(i-1) + (slopes1(1) + 2 * slopes2(1) + 2 * slopes3(1) + slopes4(1)) * h / 6;
    L(i) = L(i-1) + (slopes1(2) + 2 * slopes2(2) + 2 * slopes3(2) + slopes4(2)) * h / 6;
    I(i) = I(i-1) + (slopes1(3) + 2 * slopes2(3) + 2 * slopes3(3) + slopes4(3)) * h / 6;
    R(i) = R(i-1) + (slopes1(4) + 2 * slopes2(4) + 2 * slopes3(4) + slopes4(4)) * h / 6;
    P(i) = P(i-1) + (slopes1(5) + 2 * slopes2(5) + 2 * slopes3(5) + slopes4(5)) * h / 6;
    B(i) = B(i-1) + (slopes1(6) + 2 * slopes2(6) + 2 * slopes3(6) + slopes4(6)) * h / 6;

end
end





