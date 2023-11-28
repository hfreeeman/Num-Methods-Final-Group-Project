function slopes = SLIRP(time_index, variables, beta_max, mu_L, mu_I, e, T, day, A_p)
    % Unpack variables
    S = variables(1);
    L = variables(2);
    I = variables(3);
    R = variables(4);
    P = variables(5);

    % find needed values 
T_B = zeros(size(T))
for i = 1:size(T)
   if T(i) < 35 && T(i) > 0
       T_B(i) = 0.0002417.*(T^2.08737) .* (35 - T)^0.72859;
   else
       T_B(i) = 0;
   end
end
    T_E = -0.35968 + 0.10789.*T - 0.00214.*T.^2;
    beta = beta_max .* T_B;
    dP_l_dt = (1.33 .* day) .* T_E;
    dP_b_dt = (0.1724 .* P - 0.0000212 .* (P^2)) .* T_E .* (time_index);

    % Calculate derivatives
    dS_dt = -beta .* S .* I + (dP_b_dt + dP_l_dt) ./ A_p;
    dL_dt = beta .* S .* I - (mu_L^-1) .* L + e;
    dI_dt = (mu_L^-1) .* L - (mu_I^-1) .* I;
    dR_dt = (mu_I^-1) .* I;
    dP_dt = dP_b_dt + dP_l_dt;

    % Return slopes
    slopes = [dS_dt; dL_dt; dI_dt; dR_dt; dP_dt];
end