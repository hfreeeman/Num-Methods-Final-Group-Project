
% load "EnvironmentalForcing.mat"
load ('EnvironmentalForcing.mat');

% generate initial conditions and constant parameters
betaMax = 1 % max rate infection under ideal conditions (1/day)
mu_LMin = 6 % min length of latent period (min number of days latent)
mu_I = 10 % rate infection clears (number of days infectious)
e = 0.001 % rate of introduction from external sources
Ap = 5000 % norm factor for population ('final' plant surface area in cm^2)
Pi = 1.33 * 30 * (-0.35968 + 0.10789 * 15 - 0.00214 * 15 * 15) * 30 % inital size of the population in cm^2 (equal to the model after 30 days w/ constant temp of 15C)
Si = Pi/Ap % initial size of suseptible population (normalized)
Li = 0.01 * Si % initial fraction of the population that's latent
Ii = 0 % initial fraction of population that's infectious
Ri = muI * Ii % initial fraction of population that's recovered
Bi = 1 % initial size of berry population in cm^2 (berries grow starting 1st day of simulation)

% compute Tbeta using equation 7
if T <= 0 %C
    Tbeta = 0;
elseif T < 35 %C
    Tbeta = 0.000241 * (T^(2.06737)) * (35 - T)^0.72859;
else
    Tbeta = 0;
end

% compute beta using equation 8
beta = betaMax * Tbeta;

% compute muL using equation 9
j = 1;
for i:T 
    mu_L[i] = sum(Tbeta(j:i));
    while (mu2(i) < mu2min) % MIGHT BE < NOT > CHECK
        j = j + 1;
        mu2(i) = sum(Tbeta(j:i));
    end
end

% call SLIRP function as function handler 
(6 ODEs) = @SLIRP % DON'T KNOW HOW TO DO IT
% t and y unknown inputs
% beta, mu_L, T, tspan, mu_I, e, A inputs (soome vectors and some scalars)

% call rk4 function
% function takes output of SLIRP, tspan, and inital conditions for the 6
% ODEs as inputs.
% returns value of 6 dependent variables.

% plot the dependent variables versus time
% plot total population P and Berry Population B after normalizing by A
% (P/A and B/A)

