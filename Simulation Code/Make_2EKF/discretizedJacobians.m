function [Fx, Fb] = discretizedJacobians(func,t,xk,uk,bk,param,dt,choice)

nx = length(xk);
nu = length(uk);
nb = length(bk);
ny = 16;
if choice ==1
    n_outputs = nx;
else
    n_outputs = ny;
end

pert = 1e-5;

% Calculate Jacobian for xk for the chosen function
Fxcont = zeros(n_outputs,nx);
for ii = 1:nx
    xk_plus = xk;
    xk_minus = xk;

    xk_plus(ii) = xk_plus(ii) + pert;
    xk_minus(ii) = xk_minus(ii) - pert;

    if choice ==1 
        output_plus = func(t, xk_plus, uk, bk, param);
        output_minus = func(t, xk_minus, uk, bk, param);
    else
         output_plus = func(xk_plus, uk, bk, param);
        output_minus = func(xk_minus, uk, bk, param);
    end

    Fxcont(:, ii) = (output_plus - output_minus) / (2 * pert);
end

% Calculate Jacobian for bk for the chosen function
Fbcont = zeros(n_outputs,nb);
for ii = 1:nb
    bk_plus = bk;
    bk_minus = bk;

    bk_plus(ii) = bk_plus(ii) + pert;
    bk_minus(ii) = bk_minus(ii) - pert;

    if choice ==1
        output_plus = func(t, xk, uk, bk_plus, param);
        output_minus = func(t, xk, uk, bk_minus, param);
    else
        output_plus = func(xk, uk, bk_plus, param);
        output_minus = func(xk, uk, bk_minus, param);
    end

    Fbcont(:, ii) = (output_plus - output_minus) / (2 * pert);
end

% Discretization using Euler's method
if choice == 1
    Fx = eye(n_outputs) + Fxcont * dt; % Discretize Fk
    Fb = Fbcont * dt;                  % Discretize Bk
else
    Fx = Fxcont; %Discretize Hk
    Fb = Fbcont; % Discreteize Dk
end