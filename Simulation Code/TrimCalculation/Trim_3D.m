function  out =  Trim_3D(x)


global H0 Vtrim k1 gammatrim

psi0 = 0;
g0 = 9.81;

% Vt - trim velocity
% phid - desired phi
% gm - desired gamma (flight path angle)

p     = x(1);
q     = x(2);
r     = x(3);
u     = x(4);
v     = x(5);
w     = x(6);
phi   = x(7);
theta = x(8);
psi   = x(9);
h     = x(10);
dE    = x(11);
dA    = x(12);
dR    = x(13);
dT    = x(14);


Va      =   sqrt(u^2 + v^2 + w^2);
alpha   =   atan2(w,u);
beta    =   asin(v/Va);

%Unit conversions
m2f = 3.2808;
kg2s = 0.0685218;

cbar = 1.05/m2f; %ft to m
bSpan = 6.97/m2f; %ft to m
Sref = cbar*bSpan;   %reference area (meters^2)
mass = 0.3357/kg2s; %mass in slugs to kg
Ix = 0.403/m2f^2/kg2s; %slug ft^2 to kg m^2
Iy = 0.317/m2f^2/kg2s; %slug ft^2 to kg m^2
Iz = 0.591/m2f^2/kg2s; %slug ft^2 to kg m^2
Ixz = 0.049/m2f^2/kg2s; %slug ft^2 to kg m^2

Rho = std_atm(h+H0);


load param_CZ150.mat param;

% Get actuator deflections
[deltaE, deltaA, deltaR] = model_actuator_CZ150(dE,dA,dR);
%Get rotation rate from throttle model
n = model_throttle_CZ150(dT);

%Propeller Model from Simmons2023
D = 16/12/m2f; %radius of propeller (16 inch diameter)
J = n*D/Va;
J0 = 2;
Jc = J-J0;

%Nondimensional angular rates
phat = p*bSpan/(2*Va);
qhat = q*cbar/(2*Va);
rhat = r*bSpan/(2*Va);

CX = [alpha^2 Jc Jc^2 deltaE*alpha 1]*param(1:5);
CY = [beta phat rhat deltaA deltaR 1]*param(6:11);
CZ = [alpha qhat deltaE 1]*param(12:15);

CL = [beta phat rhat deltaA 1]*param(16:20);
CM = [alpha alpha^3 qhat deltaE qhat*deltaE 1]*param(21:26);
CN = [beta phat rhat deltaA deltaR 1]*param(27:32);

Fx      =   1/2*Rho*Sref*Va^2*CX;
Fy      =   1/2*Rho*Sref*Va^2*CY;
Fz      =   1/2*Rho*Sref*Va^2*CZ;
L       =   1/2*Rho*Sref*Va^2*bSpan*CL;
M       =   1/2*Rho*Sref*Va^2*cbar*CM;
N       =   1/2*Rho*Sref*Va^2*bSpan*CN;

% EQUATIONS OF MOTION

udot = -q*w + r*v + 1/mass*Fx - g0*sin(theta);
vdot = -r*u + p*w + 1/mass*Fy + g0*cos(theta)*sin(phi);
wdot = -p*v + q*u + 1/mass*Fz + g0*cos(theta)*cos(phi);
pdot = (Iz*L+Ixz*N+(Iy*Iz-Iz^2-Ixz^2)*q*r+Ixz*(Ix-Iy+Iz)*p*q)/(Ix*Iz-Ixz^2);
qdot = (M+(Iz-Ix)*p*r+Ixz*(r^2-p^2))/Iy;
rdot = (Ix*N+Ixz*L+(Ix^2-Ix*Iy+Ixz^2)*p*q-Ixz*(Ix-Iy+Iz)*q*r)/(Ix*Iz-Ixz^2);
phidot = p + (q.*sin(phi)+r.*cos(phi)).*tan(theta);
thetadot = q.*cos(phi) - r.*sin(phi);
psidot = (q.*sin(phi) + r.*cos(phi)).*sec(theta);
xFdot = u.*cos(theta).*cos(psi-psi0) + ...
    v.*(sin(phi).*sin(theta).*cos(psi-psi0) - cos(phi).*sin(psi-psi0)) + ...
    w.*(cos(phi).*sin(theta).*cos(psi-psi0) + sin(phi).*sin(psi-psi0));
yFdot = u.*cos(theta).*sin(psi-psi0) + ...
    v.*(sin(phi).*sin(theta).*sin(psi-psi0) + cos(phi).*cos(psi-psi0)) + ...
    w.*(cos(phi).*sin(theta).*sin(psi-psi0) - sin(phi).*cos(psi-psi0));
hdot = u.*sin(theta) - v.*sin(phi).*cos(theta) - w.*cos(phi).*cos(theta);

gamma = atan2(hdot,xFdot);

phitrim = atan(k1*Vtrim^2/g0);

% to find out = 0
out = [pdot; qdot; rdot; udot; vdot; wdot; phidot; thetadot; (psidot-k1*Vtrim);...
    (hdot-Vtrim*sin(gammatrim)); (phi-phitrim); (Va-Vtrim); h; psi];   % gamma - gammatrim


end