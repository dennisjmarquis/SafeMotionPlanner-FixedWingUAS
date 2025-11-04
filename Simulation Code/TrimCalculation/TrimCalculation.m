close all; clear; clc;

filepath = '../Utility'; addpath(filepath)

global H0 Vtrim k1 gammatrim

H0 = 630.0; % trim altitude
Vtrim = 21;  % trim velocity

k1max = 0.02;
gammamax = pi/15;

k1val = linspace(-k1max,k1max,11);       % 1/ROC
gmval = linspace(-gammamax,gammamax,31); % FPA


Ni=numel(k1val);
Nj=numel(gmval);
for ii=1:Ni
    for jj=1:Nj
    
    k1 = k1val(ii);
    gammatrim = gmval(jj);

    % states: p,q,r,u,v,w,phi,theta,psi,h,de,da,dr,dt
    x0 =[zeros(3,1);Vtrim;0;0;0;gammatrim+0.02;0;0;0;0;0;0.1]; 


    optfun = @(x) Trim_3D(x);
    [x,xval] = fsolve(optfun,x0,optimoptions('fsolve','Display','final'));
    fval = sum(xval);

    Xtrim(:,(ii-1)*Nj+jj)=[x; k1; gammatrim; fval];
    
    
    p = x(1); q = x(2); r = x(3); u = x(4); v = x(5); w = x(6); phi = x(7);
    theta = x(8); psi = x(9); h = x(10); dE = x(11); dA = x(12); dR = x(13); dT = x(14);
    

    pTrim(ii,jj) = p;
    qTrim(ii,jj) = q;
    rTrim(ii,jj) = r;
    uuTrim(ii,jj) = u;
    vTrim(ii,jj) = v;
    wTrim(ii,jj) = w;
    phiTrim(ii,jj) = phi;
    thetaTrim(ii,jj) = theta;
    psiTrim(ii,jj) = psi;
    hTrim(ii,jj) = h;
    dETrim(ii,jj) = dE;
    dATrim(ii,jj) = dA;
    dRTrim(ii,jj) = dR;
    dTTrim(ii,jj) = dT;
    
    
    if ii+jj == 2
        trimMax = x;
        trimMin = x;     
    else
        trimMax(x>trimMax) = x(x>trimMax);
        trimMin(x<trimMin) = x(x<trimMin);
    end
    
    end
end

plot(Xtrim(17,:))
ylabel('fsolve error')

[trimMin trimMax]

save LevelTrim Xtrim k1max gammamax k1val gmval pTrim qTrim rTrim uuTrim vTrim wTrim phiTrim thetaTrim psiTrim hTrim dETrim dATrim dRTrim dTTrim trimMax trimMin
rmpath(filepath)  