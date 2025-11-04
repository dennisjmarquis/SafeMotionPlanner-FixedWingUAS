function z = RightPerformance(x,uc,d,z_wt, flag)

global H0;
% STATES:
p = x(1); q = x(2); r = x(3);           %roll/pitch/yaw rate
u = x(4); v = x(5); w = x(6);           %u/v/w velocity components
phi = x(7); theta = x(8); psi = x(9);   %euler angles
xF = x(10); yF = x(11); h = x(12);      %position and altitude

R_b2v = rotation_Body2Earth([phi;theta;psi]); % rotation matrix from body to virtual frame
R_v2b = R_b2v';

% INPUTS:
if flag == 1
    dE = uc(1); %elevator command (pwm)
    dA = uc(2); %aileron command (pwm)
    dR = uc(3); %rudder command (pwm)
    dT = uc(4); %thrust command (pwm)
elseif flag == 2  || flag == 3
    dE = x(13);              %elevator dynamics
    dA = x(14);              %aileron dynamics
    dR = x(15);              %rudder dynamics
    dT = x(16);              %throttle dynamics
end

% Dryden dynamics
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
    [~,~,~,~,~,~,Cu,Cv,Cw,Du,Dv,Dw] = DrydenParams(h+H0,Vt1,magCase);
    fts2ms = 0.3048;
    uWind = fts2ms*(Cu*dryU + Du*windU);
    vWind = fts2ms*(Cv*dryV + Dv*windV);
    wWind = fts2ms*(Cw*dryW + Dw*windW);
    windBody = R_v2b*[uWind; vWind; wWind] + R_v2b*[windN; windE; windH];
else
    % disturbance
    windBody = d(1:3,1);
end

uw = windBody(1);
vw = windBody(2);
ww = windBody(3);

% alpha = atan2(w,u);
% Va = sqrt(u^2 + v^2 + w^2);
% beta = asin(v/Va);

alpha = atan2(w-ww, u-uw);
Va = sqrt((u-uw)^2 + (v-vw)^2 + (w-ww)^2);
beta = asin((v-vw)/Va);

% PENALTIES
pPen = z_wt(1);
qPen = z_wt(2);
rPen = z_wt(3);
VPen = z_wt(4); 
phiPen = z_wt(5);
thetaPen = z_wt(6);
psiPen = z_wt(7);
hPen = z_wt(8); 
xFPen = z_wt(9); 
yFPen = z_wt(10); 
dEPen = z_wt(11);
dAPen = z_wt(12);
dRPen = z_wt(13);
dTPen = z_wt(14);



% performance vector
z = [p*pPen; q*qPen; r*rPen; Va*VPen ; phi*phiPen; theta*thetaPen; psi*psiPen;
     h*hPen; xF*xFPen; yF*yFPen; dE*dEPen; dA*dAPen; dR*dRPen; dT*dTPen];


end
