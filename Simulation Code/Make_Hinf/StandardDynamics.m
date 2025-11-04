
function dx = StandardDynamics(x,uc,d,param,psi0,flag)

global H0;

if nargin<6 
    psi0 = 0;
end

% STATES:
p = x(1); q = x(2); r = x(3);         %roll/pitch/yaw rate
u = x(4); v = x(5); w = x(6);         %u/v/w velocity components
phi = x(7); theta = x(8); psi = x(9); %attitude
xF = x(10); yF = x(11); h = x(12);    %#ok<NASGU>

R_b2v = rotation_Body2Earth([phi;theta;psi]); % rotation matrix from body to virtual frame
R_v2b = R_b2v';


%% actuator dynamics
if flag == 1
    dE = uc(1); %elevator command (pwm)
    dA = uc(2); %aileron command (pwm)
    dR = uc(3); %rudder command (pwm)
    dT = uc(4); %thrust command (pwm)
       
elseif flag == 2 || flag == 3
    
    dE = x(13);            %elevator dynamics
    dA = x(14);             %aileron dynamics
    dR = x(15);             %rudder dynamics
    dT = x(16);

    dEc = uc(1); %elevator command (pwm)
    dAc = uc(2); %aileron command (pwm)
    dRc = uc(3); %rudder command (pwm)
    dTc = uc(4); %thrust command (pwm)

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
end

%% Dryden dynamics
if flag == 3
    dryU = x(17);                         %Dryden states
    dryV = [x(18); x(19)];
    dryW = [x(20); x(21)];
    
    windU = d(1);   %random variable to drive u disturbance
    windV = d(2);   %random variable to drive v disturbance
    windW = d(3);   %random variable to drive w disturbance
    magCase = 1;       %gust magnitude identifier
    windN = d(end-1);  %north wind component
    windE = d(end);    %east wind component
    windH = 0;   %vertical wind component
    
    % DRYDEN CALCULATIONS:
    Vt1 = sqrt(u^2+v^2+w^2);
    [Au,Av,Aw,Bu,Bv,Bw,Cu,Cv,Cw,Du,Dv,Dw] = DrydenParams(h+H0,Vt1,magCase);
    dryUdot = Au*dryU + Bu*windU;
    dryVdot = Av*dryV + Bv*windV;
    dryWdot = Aw*dryW + Bw*windW;
    fts2ms = 0.3048;
    uWind = fts2ms*(Cu*dryU + Du*windU);
    vWind = fts2ms*(Cv*dryV + Dv*windV);
    wWind = fts2ms*(Cw*dryW + Dw*windW);
    windBody = R_v2b*[uWind; vWind; wWind] + R_v2b*[windN; windE; windH];
    dryDyn = [dryUdot; dryVdot; dryWdot];    
else
    % disturbance
    windBody = d(1:3,1);
end

% V/AOA/AOS:
ua = u - windBody(1);
va = v - windBody(2);
wa = w - windBody(3);
Va = sqrt(ua^2+va^2+wa^2);

alpha = atan2(wa,ua);
beta = asin(va/Va);

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

if flag == 1
    dx = [pdot;qdot;rdot;udot;vdot;wdot;phidot;thetadot;psidot;xFdot;yFdot;hdot]; 
elseif flag == 2
    dx = [pdot;qdot;rdot;udot;vdot;wdot;phidot;thetadot;psidot;xFdot;yFdot;hdot;actdot];
elseif flag == 3
    dx = [pdot;qdot;rdot;udot;vdot;wdot;phidot;thetadot;psidot;xFdot;yFdot;hdot;actdot;dryDyn];
end

