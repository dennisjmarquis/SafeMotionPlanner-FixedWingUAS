load ../Utility/param_CZ150_pNoise.mat

%Unit conversions
m2f = 3.2808;
kg2s = 0.0685218;

% Aircraft parameters
cbar = 1.05/m2f; %ft to m
bSpan = 6.97/m2f; %ft to m
S = 7.32/m2f^2; %ft^2 to m^2
mass = 0.3357/kg2s; %mass in slugs to kg, measured for CZ150, matches with Simmons2023
g = 9.81; %m/s^2
Ix = 0.403/m2f^2/kg2s; %slug ft^2 to kg m^2
Iy = 0.317/m2f^2/kg2s; %slug ft^2 to kg m^2
Iz = 0.591/m2f^2/kg2s; %slug ft^2 to kg m^2
Ixz = 0.049/m2f^2/kg2s; %slug ft^2 to kg m^2

rho = 1.1526;     %air density (kg/m^3) at 630 m altitude 

Va = 21; %Nominal airspeed of controller (also average airspeed measured during training was 22.2m/s)

I = diag([Ix Iy Iz]);

Q_va = 1/2/mass*rho*Va^2*S*stdC(1:3).^2';

Q_omega = inv(I)*1/2*rho*Va^2*diag([bSpan cbar bSpan])*stdC(4:6).^2';