function [t, x] = rungeKutta4(f, tspan, x0, u, b, param, nRK)
    % Runge-Kutta 4th order method for solving ODEs with additional parameters
    % and a specified number of steps per sample interval.
    %
    % Inputs:
    % f     - function handle representing the state derivative, f(t, x, u, b, param)
    % tspan - vector [t0 tf] where t0 is the initial time and tf is the final time
    % x0    - initial state (column vector)
    % u     - input (column vector or matrix where each column is u at a time step)
    % b     - parameter (column vector)
    % param - additional parameters (column vector)
    % nRK   - number of Runge-Kutta steps per sample interval
    %
    % Outputs:
    % t     - time vector
    % x     - solution matrix, each row corresponds to a state at a time step

    % Initial time and final time
    t0 = tspan(1);
    tf = tspan(2);
    
    % Time vector for the input u
    time_u = linspace(t0, tf, size(u, 2));

    % Time step size
    h = (tf - t0) / nRK;
    
    % Preallocate arrays for efficiency
    t = linspace(t0, tf, nRK+1)';
    x = zeros(length(x0), nRK+1);
    
    % Set initial condition
    x(:, 1) = x0;
    
    % RK4 algorithm
    for i = 1:nRK
        % Current time and state
        ti = t(i);
        xi = x(:, i);
        
        % Compute the RK4 slopes
        k1 = f(ti, xi, u, b, param);
        k2 = f(ti + h/2, xi + h/2 * k1, u, b, param);
        k3 = f(ti + h/2, xi + h/2 * k2, u, b, param);
        k4 = f(ti + h, xi + h * k3, u, b, param);
        
        % Update the state
        x(:, i+1) = xi + h/6 * (k1 + 2*k2 + 2*k3 + k4);
    end
end