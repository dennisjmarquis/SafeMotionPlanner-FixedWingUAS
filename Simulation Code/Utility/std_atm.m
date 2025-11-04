


function [DENS]=std_atm(h)  % height in m
Rho0 = 1.225;
T0   = 288.15;
R   = 8.31432;
M   = 0.0289644;
g0   = 9.80665;
L   = 0.0065;
T   = T0-h*L;
Rho = Rho0*(T/T0)^(g0*M/R/L-1);
DENS=Rho;
end