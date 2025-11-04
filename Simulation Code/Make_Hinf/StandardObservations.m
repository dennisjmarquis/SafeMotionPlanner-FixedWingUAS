function y = StandardObservations(x,uc,d,flag)

global H0;

% STATES:
p = x(1); q = x(2); r = x(3);               %roll/pitch/yaw rate
u = x(4); v = x(5); w = x(6);               %u/v/w velocity components
phi = x(7); theta = x(8); psi = x(9);       %euler angles
xF = x(10); yF = x(11); h = x(12);          %positions

R_b2v = rotation_Body2Earth([phi;theta;psi]); % rotation matrix from body to virtual frame
R_v2b = R_b2v';

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

noise = d(4:13,1);

% V/AOA/AOS:
ua = u - windBody(1);
va = v - windBody(2);
wa = w - windBody(3);
Va = sqrt(ua^2+va^2+wa^2);


% MEASUREMENTS
y = [p; q; r; Va; phi; theta; psi; xF; yF; h] + noise ;
