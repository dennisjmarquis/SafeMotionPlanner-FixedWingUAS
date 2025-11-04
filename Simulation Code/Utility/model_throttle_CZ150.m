function n = model_throttle_CZ150(tPWM)
%Inputs: tPWM : centered actuator PWM in ms (0 to 1)
%Outputs: n : rotation rate of motor in rev/s

% Obtained from LS fit applied to multisine in 20240503_CZ150-4_Flt4

poly_thr = [6.56E3 2.91E3];

RPM = polyval(poly_thr,tPWM);
n = RPM/60;


end