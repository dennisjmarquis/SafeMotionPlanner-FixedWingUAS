clear; clc; close all; 
tic

filepath = 'Utility'; addpath(filepath)
filepath = 'Make_2EKF'; addpath(filepath)

global H0 windBody
d2r = pi/180; r2d = 1/d2r; H0 = 630.0;

load param_CZ150.mat
load TrimCalculation/CircularTrim.mat Xtrim k1max gammamax k1val gmval
load TrimCalculation/trimFitCoeff.mat

load param_CZ150_pNoise.mat

load Make_Hinf/cont_HINF.mat
contHinf = contHinfTest;

% Simulation Time
dt = 0.04;
time = 0:dt:10; 
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

Nsim = 1;

for kk = 1:Nsim
tic;

%Make Process Noise
sf = 3;
pNoise = zeros(6,pLen);
for ii = 1:6
    dC = zeros(1,pLen);
    dC(1) = sf*(2*stdC(ii)*rand-stdC(ii)); %generate first value within pm sf*std
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
simParam.magCase     = 1; % gust magntude identifier 1/2/3
windMag = 5; windAz = 360*rand; windEle = 20*rand-10; windUnit = [cosd(windEle)*cosd(windAz);cosd(windEle)*sind(windAz);sind(windEle)];
simParam.constWind   = windMag*windUnit; %NORTH/EAST/DOWN
simParam.delaymin    = 0.038;
simParam.delaymax    = 0.042;

%Generate Reference Path
[Psiref, Xref, Yref, Zref, k1ref, gamref] = generateReferencePath(time,Xtrim,param);


error_norm_mean = 0.3;      % m/s
error_norm_std = 0.1; % m/s
error_norm = error_norm_mean + error_norm_std * randn;
error_direction = [randn(2, 1);0];       % Random direction in 3D space
error_direction = error_direction / norm(error_direction);
wind_error = error_norm * error_direction;
adjusted_wind = simParam.constWind + wind_error;
Xref = Xref + adjusted_wind(1) * time';
Yref = Yref + adjusted_wind(2) * time';
Zref = Zref + adjusted_wind(3) * time';

reference.Psiref = Psiref; reference.Xref = Xref; reference.Yref = Yref; reference.Zref = Zref; reference.k1ref = k1ref; reference.gamref = gamref;
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
ind_k1_trim = 9;
ind_gam_trim =16;
ind_trim = (ind_k1_trim-1)*N_gm+ind_gam_trim;
Trim = Xtrim(:,ind_trim);

%Simulate Hinf
[xOut, uOut, yOut, dy, xc, mu_r, dist, noise, noise2, pathError, yOutEKF] = simHinf(time,x0,Trim,contHinf,simParam,param,reference);
posError = vecnorm([xOut(10,:)-Xref'; xOut(11,:)-Yref'; -xOut(12,:)-Zref'])';
outHinf.xOut{kk} = xOut;
outHinf.uOut{kk} = uOut;
outHinf.mu_r{kk} = mu_r;  %disp(['Hinf Control effort  : ', num2str(mu_r)]);
outHinf.erms{kk} = rms(posError);  
outHinf.emax{kk} = max(posError);  

%Simulate ATEKF (Kim 2009)
nx_EKF = 12; nb_EKF = 3; nm_EKF = 16;
controlhist = uOut;
Pxb0 = ones(nx_EKF,nb_EKF); Px0 = 10*eye(nx_EKF); Pb0 = 10*eye(nb_EKF); V0 = Pxb0/Pb0; b0 = zeros(3,1);
xhist(:,1) = [xOut(10:12,1);xOut(1:9,1)]; %x y h p q r u v w phi theta psi
xbarhist(:,1) = xhist(:,1)-V0*b0;
bhist(:,1) = b0;
Vhist(:,:,1) = V0;
Px_barhist(:,:,1) = Px0-V0*Pb0*V0';
Pb_barhist(:,:,1) = Pb0;
M = 25; %Window size for covariance estimates
etaxM = zeros(nm_EKF,M);
etabM = zeros(nm_EKF,M);
lambdax = 1; lambdab = 1;
yhist = [xOut(10:12,:);yOutEKF];
Q = diag([zeros(1,3) 0.0451 0.0227 0.0053 0.0026 0.0101 0.1190 zeros(1,3)]);
Qb = 9e-6*eye(3);
R = diag([1e-3 1e-3 1e-4 7.6E-5 7.6E-5 7.6E-5 0.01 0.01 0.01 4E-5 4E-5 9E-3 4 0.2 0.5 0.5]); %Values on live EKF

for ii=1:pLen-1
[xhist(:,ii+1), xbarhist(:,ii+1), bhist(:,ii+1),Vhist(:,:,ii+1),Px_barhist(:,:,ii+1),Pb_barhist(:,:,ii+1),Pxhist(:,:,ii),~,etaxM,etabM,lambdax,lambdab] = myATEKF(xhist(:,ii), xbarhist(:,ii), bhist(:,ii), controlhist(:,ii), controlhist(:,ii+1), yhist(:,ii+1), Vhist(:,:,ii), Px_barhist(:,:,ii), Pb_barhist(:,:,ii), etaxM, etabM, lambdax, lambdab, Q, Qb, R, param,dt);
end


ttt = toc;
fprintf('Simulation %d took %f seconds\n',kk,ttt)

referenceArray{kk} = reference;
simParamArray{kk} = simParam;

figure
hold on;
plot3(outHinf.xOut{kk}(11,:),outHinf.xOut{kk}(10,:),outHinf.xOut{kk}(12,:))
plot3(Yref,Xref,-Zref)
ylabel('X, forward (m)')
xlabel('Y, side (m)')
zlabel('H, height (m)')
legend('Hinf','Reference')
axis equal

windInertial=zeros(3,pLen);
for ii=1:pLen
    R_b2i = rotation_Body2Earth(xOut(7:9,ii));
    windInertial(:,ii)=R_b2i*dist(1:3,ii);
end
figure
hold on
plot(time,bhist(1,:),'r:')
plot(time,bhist(2,:),'g:')
plot(time,bhist(3,:),'b:')
plot(time, windInertial(1,:),'r')
plot(time, windInertial(2,:),'g')
plot(time, windInertial(3,:),'b')
xlabel('time')
ylabel('wind estimate')
legend('x','y','z')

estErrs(kk) = mean(vecnorm(windMag*windUnit-bhist(:,floor(pLen/2):end)));

windUnits{kk} = windUnit;
windMags{kk} = windMag;
fprintf('For wind of magnitude %d m/s, estimation error is %f m/s\n',windMag,estErrs(kk))
end

