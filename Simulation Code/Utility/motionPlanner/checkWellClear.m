function [flag_wc, tcrit,tca] = checkWellClear(Rrel,Vrel)
global Robs Dsep Va
tca = -dot(Rrel,Vrel)/norm(Vrel)^2;
dca = norm(Rrel+tca*Vrel);

if dca < Robs + Dsep && tca > 0
    flag_wc = true;
    tcrit = (norm(Rrel) - Dsep - Robs) / norm(Vrel);
else
    flag_wc = false;
    tcrit = NaN;
end
