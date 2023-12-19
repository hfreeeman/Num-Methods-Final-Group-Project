%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% SLIRPE function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SLIRPE_model (pronounced like the 7-11 slushy drink) combines all the 
% SLIR equations and the plant growth equations with an external source 
% equation in one system. The dependent variables are stored in y as:
% y(1) = B (amount of population, surface area, that is berries)
% y(2) = P (total population, surface area, including berries and leaves: P = S+L+I+R)
% y(3) = S (susceptible population)
% y(4) = L (latent population)
% y(5) = I (infectious population)
% y(6) = R (recovered/removed population)
% y(7) = E (amount of new infections from external sources)
% y(8) = F (size of the spreading population, e.g. sporulating for a fungus) 
%
% and the parameters in a cell array:
% p{1} = beta_max (max rate of colony growth/new infections)
% p{2} = mu_i (inverse length of the infectious period in days)
% p{3} = T (array of temperature in C)
% p{4} = day (array of times in units of days)
% p{5} = A (total plant surface area at reference time)
% p{6} = Windspd (windspeed)
% p{7} = Winddir (wind direction)
% p{8} = eta     (release fraction scale factor)
% p{9} = kappa   (release fraction scale factor)  
% p{10}= xi      (release fraction offset)
% p{11}= Gamma   (spore production multiple)
% p{12}= alpha   (spore production 2nd factor)
%
% and the time input (idx) should be an integer for the iteration number
% Note that E is not calculated by the function (only integrated in time)
function [dy_dt] = SLIRPE_model(idx, y0, e, mu_L, p)

 %assign parameters
 [beta_max, mu_I, T, day, A_p, WindSpd, Winddir, eta, kappa, xi, Gamma, alpha] = p{:};
    
 % Assign variables
    B = y0(1);
    P = y0(2);
    S = y0(3);
    L = y0(4);
    I = y0(5);
    R = y0(6);
    E = y0(7);
    F = y0(8);

 T_used = T(idx);
 day_used = day(idx);
 mu_L_used = mu_L(idx);
 m_used = WindSpd(idx);


    beta = beta_max*Sall_temp_effect(T_used); %pathogen growth rate

    dPl_dt = (1.33 .* (day_used + 30)) .* Sall_temp_effect(T_used);

    dy_dt = zeros(1,8);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dy_dt(1) = (0.1724 .* B - 0.0000212 .* (B^2)) .* Sall_temp_effect(T_used); % B
    dy_dt(2) = dy_dt(1) + dPl_dt; % P
    dy_dt(3) = -(beta*S*I) + (dy_dt(2) / A_p); % S
    dy_dt(4) = beta*S*I - (mu_L_used) * L + e; % L
    dy_dt(5) = (mu_L_used) * L - mu_I*I; % I
    dy_dt(6) = (mu_I) * I; % R
    dy_dt(7) = e; % E
    if I == 0
    dy_dt(8) = 0; % F
    else 
    dy_dt(8) = Gamma * exp(alpha * I * A_p / 10000) - F *(exp(kappa * m_used + xi) / (eta *(1 + exp(kappa*m_used + xi))));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
