function dx = f_CZ150EKF(~,xk,uk,bk,param)

% STATES:
xF = xk(1); yF = xk(2); hF = xk(3);
p = xk(4); q = xk(5); r = xk(6);
ua = xk(7); va = xk(8); wa = xk(9);
phi = xk(10); theta = xk(11); psi = xk(12);

dE = uk(1); dA = uk(2); dR = uk(3); dT = uk(4);
                  
% ROTATION (earth to body)
R_b2i = rotation_Body2Earth([phi;theta;psi]);


% V/AOA/AOS:
Va = sqrt(ua^2+va^2+wa^2);

alpha = atan2(wa,ua);
beta = asin(va/Va);

%Actuator Saturations
eMin = -0.4;
eMax = 0.4;
aMin = -0.4;
aMax = 0.4;
rMin = -0.3;
rMax = 0.3;

if dE > eMax
    dE = eMax;
elseif dE < eMin
    dE = eMin;
end

if dA > aMax
    dA = aMax;
elseif dA < aMin
    dA = aMin;
end

if dR > rMax
    dR = rMax;
elseif dR < rMin
    dR = rMin;
end

if dT<0
    dT = 0;
end

if dT>1
    dT = 1;
end

% Get actuator deflections
[deltaE, deltaA, deltaR] = model_actuator_CZ150(dE,dA,dR);
%Get rotation rate from throttle model
n = model_throttle_CZ150(dT);

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

%Propeller Model from Simmons2023
D = 16/12/m2f; %radius of propeller (16 inch diameter)
J = n*D/Va;
J0 = 2;
Jc = J-J0;

%Nondimensional angular rates
phat = p*bSpan/(2*Va);
qhat = q*cbar/(2*Va);
rhat = r*bSpan/(2*Va);

%Model Coefficients
CX = [alpha^2 Jc Jc^2 deltaE*alpha 1]*param(1:5);
CY = [beta phat rhat deltaA deltaR 1]*param(6:11);
CZ = [alpha qhat deltaE 1]*param(12:15);

CL = [beta phat rhat deltaA 1]*param(16:20);
CM = [alpha alpha^3 qhat deltaE qhat*deltaE 1]*param(21:26);
CN = [beta phat rhat deltaA deltaR 1]*param(27:32);
    

% CALCULATE BODY FORCES AND MOMENTS
qbar = 0.5*rho*Va^2;
Fx = CX*qbar*S;
Fy = CY*qbar*S;
Fz = CZ*qbar*S;
L = CL*qbar*S*bSpan;
M = CM*qbar*S*cbar;
N = CN*qbar*S*bSpan;

ax = Fx/mass;
ay = Fy/mass;
az = Fz/mass;

% EQUATIONS OF MOTION
Pdot = R_b2i*[ua;va;wa]+bk;
xFdot = Pdot(1);
yFdot = Pdot(2);
hdot = -Pdot(3);

pdot = (Iz*L+Ixz*N+(Iy*Iz-Iz^2-Ixz^2)*q*r+Ixz*(Ix-Iy+Iz)*p*q)/(Ix*Iz-Ixz^2);
qdot = (M+(Iz-Ix)*p*r+Ixz*(r^2-p^2))/Iy;
rdot = (Ix*N+Ixz*L+(Ix^2-Ix*Iy+Ixz^2)*p*q-Ixz*(Ix-Iy+Iz)*q*r)/(Ix*Iz-Ixz^2);

uadot = -q*wa + r*va + 1/mass*Fx - g*sin(theta);
vadot = -r*ua + p*wa + 1/mass*Fy + g*cos(theta)*sin(phi);
wadot = -p*va + q*ua + 1/mass*Fz + g*cos(theta)*cos(phi);

phidot = p + (q.*sin(phi)+r.*cos(phi)).*tan(theta);
thetadot = q.*cos(phi) - r.*sin(phi);
psidot = (q.*sin(phi) + r.*cos(phi)).*sec(theta);

dx = [xFdot;yFdot;hdot;pdot;qdot;rdot;uadot;vadot;wadot;phidot;thetadot;psidot];
end