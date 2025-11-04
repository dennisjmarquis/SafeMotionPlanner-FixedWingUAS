function [Au,Av,Aw,Bu,Bv,Bw,Cu,Cv,Cw,Du,Dv,Dw] = DrydenParams(h,V,MagCase)
% altitude (m) - convert to feet
hfeet = (h+30)*3.28084;
if hfeet < -30
    error('negative altitude in DrydenParam')
end
% wind speed at 6m (knots -> ft/s)
switch MagCase
    case 1
        w20 = 15*1.68781;
    case 2
        w20 = 30*1.68781;  
    case 3 
        w20 = 45*1.68781;
end

% velocity (m/s) - convert to ft/s
V = V*3.28084;

sigW = 0.1*w20;
sigU = sigW/((0.177+0.000823*hfeet)^0.4);
sigV = sigU;

Lw = hfeet;
Lu = hfeet/((0.177+0.000823*hfeet)^1.2);
Lv = Lu;

numU = sigU*sqrt(2*Lu/(pi*V));
denU = [Lu/V, 1];

numV = sigV*sqrt(Lv/(pi*V))*[sqrt(3)*Lv/V, 1];

%numV = sigV*sqrt(Lv/(pi*V));
denV = [(Lv/V)^2, 2*Lv/V, 1];

 numW = sigW*sqrt(Lw/(pi*V))*[sqrt(3)*Lw/V, 1];

%numW = sigW*sqrt(Lw/(pi*V));
denW = [(Lw/V)^2, 2*Lw/V, 1];

[Au,Bu,Cu,Du] = tf2ss(numU,denU);

[Av,Bv,Cv,Dv] = tf2ss(numV,denV);

[Aw,Bw,Cw,Dw] = tf2ss(numW,denW);
return