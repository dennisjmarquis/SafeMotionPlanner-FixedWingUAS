function [deltaE, deltaA, deltaR] = model_actuator_CZ150(ePWM,aPWM,rPWM)
%Inputs: aPWM, ePWM, rPWM : actuator PWM in ms
%Outputs: dealtA, deltaE, deltaR : actuator deflection in radian

%CZ150_4 Fit performed by Dennis 2/2024

%PWM ms thresholds
    % Ail : 1.1-1.9 --> -0.4 to 0.4
    % Ele : 1.1-1.9 --> -0.4 to 0.4
    % Rud : 1.2-1.8 --> -0.3 to 0.3
poly_ail = [3.981e-08 -1.082e-08 0.06235 0.2563];
poly_ele = [-8.241e-08 -9.535e-06 0.07207 -0.8866];
poly_rud = [-2.724e-07 7.581e-06 0.1156 -0.7398];

deltaE = polyval(poly_ele,ePWM*1000)*pi/180;
deltaA = -polyval(poly_ail,aPWM*1000)*pi/180;
deltaR = polyval(poly_rud,rPWM*1000)*pi/180;
end