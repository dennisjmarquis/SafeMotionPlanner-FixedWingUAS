%Modified to have pNoise, nargin = 7
function [dx,ax,ay,az] = StandardSimDynamics_CZ150(~,x,command,wind,param,pNoise,psi0)

global H0 windBody

if nargin<7
    psi0 = 0;
end

% STATES:
p = x(1); q = x(2); r = x(3);         %roll/pitch/yaw rate (rad/s)
u = x(4); v = x(5); w = x(6);         %u/v/w velocity components (m/s)
phi = x(7); theta = x(8); psi = x(9); %attitude (rad)
xF = x(10); yF = x(11); h = x(12);    %#ok<NASGU>
dE = x(13);             %elevator dynamics (PWM ms)
dA = x(14);            %aileron dynamics (PWM ms)
dR = x(15);            %rudder dynamics (PWM ms)
dT = x(16);            %throttle dynamics (PWM ms)

dryU = x(17);                         %Dryden states
dryV = [x(18); x(19)];
dryW = [x(20); x(21)];
                        

% INPUTS:
dEc = command(1); %elevator command (PWM ms)
dAc = command(2); %aileron command (PWM ms)
dRc = command(3); %rudder command (PWM ms)
dTc = command(4); %thrust command (PWM ms)
windU = wind(1);   %random variable to drive u disturbance
windV = wind(2);   %random variable to drive v disturbance
windW = wind(3);   %random variable to drive w disturbance
magCase = wind(4); %gust magnitude identifier
windN = wind(5);   %north wind component
windE = wind(6);   %east wind component
windH = wind(7);   %vertical wind component

% ROTATION (earth to body)
REarthToBody = rotation_Body2Earth([phi;theta;psi])';


% DRYDEN CALCULATIONS:
Vt1 = sqrt(u^2+v^2+w^2);
if isnan(Vt1)
    error('Airspeed is NaN');
end
[Au,Av,Aw,Bu,Bv,Bw,Cu,Cv,Cw,Du,Dv,Dw] = DrydenParams(h+H0,Vt1,magCase);
dryUdot = Au*dryU + Bu*windU;
dryVdot = Av*dryV + Bv*windV;
dryWdot = Aw*dryW + Bw*windW;
fts2ms = 0.3048;
uWind = fts2ms*(Cu*dryU + Du*windU);
vWind = fts2ms*(Cv*dryV + Dv*windV);
wWind = fts2ms*(Cw*dryW + Dw*windW);
windBody = REarthToBody*[uWind; vWind; wWind] + REarthToBody*[windN; windE; windH];
dryDyn = [dryUdot; dryVdot; dryWdot];
% disp(norm(windBody))

% V/AOA/AOS:
ua = u - windBody(1);
va = v - windBody(2);
wa = w - windBody(3);
Va = sqrt(ua^2+va^2+wa^2);

alpha = atan2(wa,ua);
beta = asin(va/Va);

%Values are centered PWM ms, based on PWM Ranges: Ele 1100-1900, Ail
%1100-1900, Rud 1200-1800
eMin = -0.4;
eMax = 0.4;
aMin = -0.4;
aMax = 0.4;
rMin = -0.3;
rMax = 0.3;

%Saturations 
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

if dEc > eMax
    dEc = eMax;
elseif dEc < eMin
    dEc = eMin;
end

if dAc > aMax
    dAc = aMax;
elseif dAc < aMin
    dAc = aMin;
end

if dRc > rMax
    dRc = rMax;
elseif dRc < rMin
    dRc = rMin;
end

if dTc<0
    dTc = 0;
end

if dTc>1
    dTc = 1;
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

rho = std_atm(h+H0); % air density at altitude h

%Propeller Model from Simmons2023
D = 16/12/m2f; %radius of propeller (16 inch diameter)
J = n*D/Va;
J0 = 2;
Jc = J-J0;

%Nondimensional angular rates
phat = p*bSpan/(2*Va);
qhat = q*cbar/(2*Va);
rhat = r*bSpan/(2*Va);

%Generate process noise
dCX = pNoise(1);
dCY = pNoise(2);
dCZ = pNoise(3);
dCL = pNoise(4);
dCM = pNoise(5);
dCN = pNoise(6);

%Model Coefficients
CX = [alpha^2 Jc Jc^2 deltaE*alpha 1]*param(1:5) + dCX;
CY = [beta phat rhat deltaA deltaR 1]*param(6:11) + dCY;
CZ = [alpha qhat deltaE 1]*param(12:15) + dCZ;

CL = [beta phat rhat deltaA 1]*param(16:20) + dCL;
CM = [alpha alpha^3 qhat deltaE qhat*deltaE 1]*param(21:26) + dCM;
CN = [beta phat rhat deltaA deltaR 1]*param(27:32) + dCN;
    

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
pdot = (Iz*L+Ixz*N+(Iy*Iz-Iz^2-Ixz^2)*q*r+Ixz*(Ix-Iy+Iz)*p*q)/(Ix*Iz-Ixz^2);
qdot = (M+(Iz-Ix)*p*r+Ixz*(r^2-p^2))/Iy;
rdot = (Ix*N+Ixz*L+(Ix^2-Ix*Iy+Ixz^2)*p*q-Ixz*(Ix-Iy+Iz)*q*r)/(Ix*Iz-Ixz^2);

udot = -q*w + r*v + 1/mass*Fx - g*sin(theta);
vdot = -r*u + p*w + 1/mass*Fy + g*cos(theta)*sin(phi);
wdot = -p*v + q*u + 1/mass*Fz + g*cos(theta)*cos(phi);

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

% ACTUATOR DYNAMICS:
tauEle = 0.055 + 0.016;%Time constants obtained from NSL with first order transfer function 1/(tau*s+1)
tauAil = 0.055 + 0.028;
tauRud = 0.055 + 0.016;
tauThr = 0.082; %Obtained from OE method applied to multisine in 20210831_CZ150-2_Flt1

dEdot = -1/tauEle*dE + 1/tauEle*dEc;
dAdot = -1/tauAil*dA + 1/tauAil*dAc;
dRdot = -1/tauRud*dR + 1/tauRud*dRc;

dTdot = -1/tauThr*dT + 1/tauThr*dTc;

actdot = [dEdot; dAdot; dRdot; dTdot];

dx = [pdot;qdot;rdot;udot;vdot;wdot;phidot;thetadot;psidot;xFdot;yFdot;hdot;actdot;dryDyn];
