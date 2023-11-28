function [t, y] = RungeKutta(func, tspan, y0, h)
    % Runge-Kutta method for solving systems of ODEs
    
    % Number of steps
    numSteps = round((tspan(2) - tspan(1)) / h);
    
    % Initialize vectors
    t = linspace(tspan(1), tspan(2), numSteps + 1);
    y = zeros(length(y0), numSteps + 1);
    y(:, 1) = y0;
    
    % Runge-Kutta method for systems
    for i = 1:numSteps
        k1 = h * func(t(i), y(:, i));
        k2 = h * func(t(i) + h/2, y(:, i) + k1/2);
        k3 = h * func(t(i) + h/2, y(:, i) + k2/2);
        k4 = h * func(t(i) + h, y(:, i) + k3);
        
        y(:, i + 1) = y(:, i) + (k1 + 2*k2 + 2*k3 + k4) / 6;
    end
end