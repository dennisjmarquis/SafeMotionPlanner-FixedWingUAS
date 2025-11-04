clear; close all; clc;

load CircularTrim


figure
mesh(gmval,k1val,pTrim);
xlabel('gamma'); ylabel('k1'); zlabel('p_{tr}')

figure
mesh(gmval,k1val,qTrim);
xlabel('gamma'); ylabel('k1'); zlabel('q_{tr}')

figure
mesh(gmval,k1val,rTrim);
xlabel('gamma'); ylabel('k1'); zlabel('r_{tr}')


figure
mesh(gmval,k1val,uuTrim);
xlabel('gamma'); ylabel('k1'); zlabel('u_{tr}')

figure
mesh(gmval,k1val,vTrim);
xlabel('gamma'); ylabel('k1'); zlabel('v_{tr}')

figure
mesh(gmval,k1val,wTrim);
xlabel('gamma'); ylabel('k1'); zlabel('w_{tr}')


figure
mesh(gmval,k1val,phiTrim);
xlabel('gamma'); ylabel('k1'); zlabel('\phi_{tr}')


figure
mesh(gmval,k1val,thetaTrim);
xlabel('gamma'); ylabel('k1'); zlabel('\theta_{tr}')


figure
mesh(gmval,k1val,dETrim);
xlabel('gamma'); ylabel('k1'); zlabel('dE_{tr}')

figure
mesh(gmval,k1val,dATrim);
xlabel('gamma'); ylabel('k1'); zlabel('dA_{tr}')

figure
mesh(gmval,k1val,dRTrim);
xlabel('gamma'); ylabel('k1'); zlabel('dR_{tr}')

figure
mesh(gmval,k1val,dTTrim);
xlabel('gamma'); ylabel('k1'); zlabel('dT_{tr}')






%% trim fit coefficient %%%Need to do this for CZ150 trims, the below are for Telemaster
%  f(x,y) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 (x=k1,y=gamma)

pfit = [0;-1.1476;0;-32.5405;-20.3307;0];
qfit = [0;0;0;726.0481;0;0];
rfit = [0;16.7388;0;0;0;0];
ufit = [20.9391;0.9053;0;-115.2370;2.2471;0.0461]; %u varies between 20.88-20.94, not worth fitting
vfit = [0.9727;-20.1476;0;-99.8007;-49.8086;0];
wfit = [1.2665;0;0;1485.6;0;0];
phfit = [0;38.6736;0;0;0;0];
thfit = [0.0569;0;0.9844;0;0;0]; %Ignored p10 for 0.97 R^2
dEfit = [0.0654;0;0;-128.0343;0;0];
dAfit = [0.0280;-1.0589;0;-7.1126;-4.1713;0];
dRfit = [0.0315;-2.0171;0;0;0;0];
dTfit = [0.5349;0;1.2584;0;0;-2.0742];

save trimFitCoeff pfit qfit rfit ufit vfit wfit phfit thfit dEfit dAfit dRfit dTfit


%% fitting check

fit = [pfit qfit rfit ufit vfit wfit phfit thfit zeros(6,2) dEfit dAfit dRfit dTfit];

for i = [1:8 11:14]
    
    fitCoeff = fit(:,i);
    
    trimVals = Xtrim(i,:);
    k1vals   = Xtrim(15,:);
    gmvals   = Xtrim(16,:);
    
    for j = 1:5000
        ind = randi([1,numel(k1vals)],1);
        
        k1 = k1vals(ind);
        gm = gmvals(ind);
        
        trimEst = sum(fitCoeff.*[1; k1; gm; k1^2; gm*k1; gm^2]);
        
        err(i,j) = trimVals(ind) - trimEst;
        
    end
    
end

mean(err')
max(err')