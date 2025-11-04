% clear; clc; close all; 
tic

filepath = 'Utility'; addpath(filepath)
filepath = 'Utility/motionPlanner'; addpath(filepath)
filepath = 'Make_2EKF'; addpath(filepath)

global H0 windBody
global Rrod Robs Dsep Va Vb k1max g1max gammamax deltatTrajec tcritmax Xgoal F wvec flyMode tfin
d2r = pi/180; r2d = 1/d2r; H0 = 630.0;

load param_CZ150.mat
load TrimCalculation/CircularTrim.mat Xtrim k1max gammamax k1val gmval
load TrimCalculation/trimFitCoeff.mat
load primArrayExtreme.mat


%Load Process Noise
load param_CZ150_pNoise.mat

%Load Controller
load Make_Hinf/cont_HINF_sw_UIC.mat
contHinf = contHinfSwitchedUIC;

% Simulation Time
dt = 0.04;
time = 0:dt:20; 
pLen = numel(time); 
N_k1 = 11;
N_gm = 31;

% Disturbance Selection
simParam.gustToggle    = 1;
simParam.noiseToggle   = 1;
simParam.pNoiseToggle  = 1;   
simParam.delayToggle   = 1;
simParam.windToggle    = 1;
simParam.controlToggle = 1;

% Motion Planner parameters -------------------------------------------------
Rrod = 185; %radius of detection (m)
Robs = 5; %radius of obstacles (m)
Dsep = 35; %minimum separation distance (m)
Va = 21; %ownship velocity (m/s)
k1max = 0.02; %max curvature of ownship turns
gam1max = pi/15;
phimax = atand(k1max*Va^2/9.81); %max bank angle of ownship
tcritmax = 6; %time until collision where action should be taken

%Geofence initial conditions
Pgeo=[-1450 1450 1450 -1450;...
      -1450 -1450 1450 1450]; %x and y coordinates of geofence endpoints
% Identify hyperplanes of geofence [px py mx my]
ngeo = size(Pgeo,2);
F = zeros(ngeo,4);
indexGeo = [1:ngeo 1];
for ii=1:ngeo
    pseg = (Pgeo(:,indexGeo(ii+1))+Pgeo(:,indexGeo(ii)))/2;
    mOrthogseg = [Pgeo(1,indexGeo(ii+1))-Pgeo(1,indexGeo(ii)) Pgeo(2,indexGeo(ii+1))-Pgeo(2,indexGeo(ii))];
    mseg = [-mOrthogseg(2) mOrthogseg(1)];
    mseg = mseg/norm(mseg);
    F(ii,1)=pseg(1);
    F(ii,2)=pseg(2);
    F(ii,3)=mseg(1);
    F(ii,4)=mseg(2);
end

%Cost function weighing
wvec = [1/3 1/3 1/3];

Nsim = 1;

for kk = 1:Nsim
tic;

%Obstacle and goal initial conditions
obsIC=[pi/6*rand-pi/12 pi/6*rand-pi/12 pi-pi/12+pi/6*rand pi/6*rand-pi/12 16.2]; %Load Obstacle IC, [theta phi psi gamma Vb]
Nconfig=size(obsIC,1);
Xgoal = [60*21;0;0;0;0]; %x y z psi gamma
numObs = 1; %number of obstacles
simWC = false;
    while simWC==false
   obsIC=[pi/6*rand-pi/12 pi/6*rand-pi/12 pi-pi/12+pi/6*rand pi/6*rand-pi/12 16.2]; %Load Obstacle IC, [theta phi psi gamma Vb]
    randConfig=randi(Nconfig);
    theta_obs = obsIC(randConfig,1);
    phi_obs = obsIC(randConfig,2);
    psi_obs = obsIC(randConfig,3);
    gamma_obs = obsIC(randConfig,4);
    Vb = obsIC(randConfig,5); %obstacle velocity (m/s)
    Xb0 = [(Rrod+Robs+5)*cos(theta_obs)*cos(phi_obs);(Rrod+Robs+5)*sin(theta_obs)*cos(phi_obs);(Rrod+Robs+5)*sin(phi_obs);wrapToPi(psi_obs);gamma_obs]; %x y z psi gamma
    Rrel0 = Xb0(1:3);
    Vrel0=[Vb*cos(psi_obs)*cos(gamma_obs)-Va; Vb*sin(psi_obs)*cos(gamma_obs);Vb*sin(gamma_obs)];
    simWC = checkWellClear(Rrel0,Vrel0);
    end

%Make Process Noise
sf = 3;
pNoise = zeros(6,pLen);
for ii = 1:6
    dC = zeros(1,pLen);
    dC(1) = sf*(2*stdC(ii)*rand-stdC(ii));
    for jj = 2:pLen
        sMax = min(dC(jj-1)+derivCmax(ii)*dt,sf*stdC(ii));
        sMin = max(dC(jj-1)+derivCmin(ii)*dt,-sf*stdC(ii));
        dC(jj) = sMin + (sMax-sMin)*rand;
    end
        pNoise(ii,:)=dC;
end

%Generate Simulation Disturbance
simParam.dn          = randn(13,pLen);
simParam.dn2         = randn(6,pLen); %used for EKF u v w ax ay az
simParam.pNoise      = pNoise;
simParam.magCase     = 3; % gust magntude identifier 1/2/3
windMag = 3; windAz = 360*rand; windEle = 20*rand-10; windUnit = [cosd(windEle)*cosd(windAz);cosd(windEle)*sind(windAz);sind(windEle)];
simParam.constWind   = windMag*windUnit; %NORTH/EAST/DOWN
simParam.delaymin    = 0.038;
simParam.delaymax    = 0.042;

%Generate Level Reference Path
leveltrim   = Xtrim(:,171);
Xa0=[0;0;0;0;0];
deltatTrajec = dt;
tTrajec = time;
Ntrajec = length(tTrajec);
flyMode = 'linear';
tfin=0;
trajec = genTrajecLinear(Xa0,tTrajec,leveltrim);

Xref   = trajec(2,:)';
Yref   = trajec(3,:)';
Zref   = trajec(4,:)';
Varef = trajec(5,:)';
Psiref = trajec(6,:)';
k1ref  = trajec(7,:)';
gamref = trajec(8,:)';
pref = trajec(9,:)';
qref = trajec(10,:)';
rref = trajec(11,:)';
phiref = trajec(12,:)';
thetaref = trajec(13,:)';
dEref = trajec(14,:)';
dAref = trajec(15,:)';
dRref = trajec(16,:)';
dTref = trajec(17,:)';

reference.Xref = Xref; reference.Yref = Yref; reference.Zref = Zref; reference.Varef = Varef; reference.Psiref = Psiref; reference.k1ref = k1ref; reference.gamref = gamref;
reference.pref = pref; reference.qref = qref; reference.rref = rref; reference.phiref = phiref; reference.thetaref = thetaref; reference.dEref = dEref; reference.dAref = dAref; reference.dRref = dRref; reference.dTref = dTref;
reference.k1val = k1val; reference.gmval = gmval;
reference.noPosFlag = 0; %Set to 1 to track attitude only, used for straight and level

%Initial Condition
ind_k1_level = 6;
ind_gam_level = 16;
ind_level = (ind_k1_level-1)*N_gm+ind_gam_level; %171: straight and level
levelTrim = Xtrim(:,ind_level);
uTrim     = levelTrim(11:14);
x0  = [levelTrim(1:9); 0; 0; levelTrim(10,1); uTrim; zeros(5,1)];

%Trim for desired maneuver
ind_k1_trim = 1;
ind_gam_trim =16;
ind_trim = (ind_k1_trim-1)*N_gm+ind_gam_trim;
Trim = Xtrim(:,ind_trim);

%Simulate Hinf UIC
[xOut, uOut, yOut, mu_r, dist, noise, reference,Xbvec] = simHinfSwitchedUIC_MP(time,x0,Xtrim,contHinf,simParam,param,reference,Xb0,primArray,trajec,leveltrim);
posError = vecnorm([xOut(10,:)-reference.Xref(1:pLen)'; xOut(11,:)-reference.Yref(1:pLen)'; xOut(12,:)+reference.Zref(1:pLen)'])';
outHinf.xOut{kk} = xOut;
outHinf.uOut{kk} = uOut;
outHinf.mu_r{kk} = mu_r;
outHinf.emean{kk} = mean(posError);  
disp(['Hinf Mean Position Error : ' num2str(mean(posError))])
outHinf.emax{kk} = max(posError);  
disp(['Hinf Max Position Error : ' num2str(max(posError))])

ttt = toc;
fprintf('Simulation %d took %f seconds\n\n',kk,ttt)

referenceArray{kk} = reference;
simParamArray{kk} = simParam;

windUnits{kk} = windUnit;
windMags{kk} = windMag;

for ii=1:pLen
    outDist(ii)=norm([xOut(10:11,ii);-xOut(12,ii)]-Xbvec(1:3,ii),2);
    refDist(ii)=norm([reference.Xref(ii);reference.Yref(ii);reference.Zref(ii)]-Xbvec(1:3,ii),2);
end
minRef(kk)=min(refDist);
minOut(kk)=min(outDist);
XbInfo.Xbvec{kk}=Xbvec;
XbInfo.Xb0(:,kk)=Xb0;
XbInfo.Vb(kk)=Vb;


end
