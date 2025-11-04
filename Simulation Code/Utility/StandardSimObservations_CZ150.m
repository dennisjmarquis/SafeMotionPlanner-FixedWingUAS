function y = StandardSimObservations_CZ150(x,wind,noise)
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
global H0;

% INPUTS:
windU = wind(1);   %random variable to drive u disturbance
windV = wind(2);   %random variable to drive v disturbance
windW = wind(3);   %random variable to drive w disturbance
magCase = wind(4); %gust magnitude identifier
windN = wind(5);   %north wind component
windE = wind(6);   %east wind component
windH = wind(7);   %vertical wind component


% ROTATION (earth to body)
REarthToBody = rotation_Body2Earth([phi;theta;psi])';

% DRYDEN OBSERVATION CALCULATIONS
Vt = sqrt(u^2+v^2+w^2);
[~,~,~,~,~,~,Cu,Cv,Cw,Du,Dv,Dw] = DrydenParams(h+H0,Vt,magCase);
fts2ms = 0.3048;
uWind = fts2ms*(Cu*dryU + Du*windU);
vWind = fts2ms*(Cv*dryV + Dv*windV);
wWind = fts2ms*(Cw*dryW + Dw*windW);
windBody = REarthToBody*[uWind; vWind; wWind] + REarthToBody*[windN; windE; windH];

% V/AOA/AOS:
ua = u - windBody(1);
va = v - windBody(2);
wa = w - windBody(3);
Va = sqrt(ua^2+va^2+wa^2);

AOA = atan2(wa,ua);
AOS = asin(va/Va);


% MEASUREMENTS
y = [p; q; r; Va; phi; theta; psi; xF; yF; h] + noise;
