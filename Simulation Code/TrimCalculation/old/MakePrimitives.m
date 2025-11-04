filepath = '../Utility'; addpath(filepath)


global H0 Vtrim k1 gammatrim
H0 = 630;

load CircularTrim Xtrim
load param_CZ150.mat
load primLevel.mat

d2r = pi/180; r2d = 1/d2r;

kindList=[1 3 6 9 11];
gamindList= [1 8 16 24 31];
tList=[2.5 5];
[K G T] = ndgrid(kindList,gamindList,tList);
indList = [K(:) G(:) T(:)];

Nev = size(indList,1);

for kk = 1:Nev

indK1 = indList(kk,1);
indGam = indList(kk,2);    

% Simulation time ---------------------------------------------------------
timestep = 0.04;
time     = 0:timestep:indList(kk,3); %simulation time
pLen     = numel(time); %number of simultion
jj        = 1;

Nk1 = 11; Ngam = 31;
ind      = (indK1-1)*Ngam+indGam;

k1 = Xtrim(15,ind);
gammatrim = Xtrim(16,ind);
Vtrim = norm(Xtrim(4:6,ind));
disp(['Vtrim, k1, and FPA are ', num2str(Vtrim), ', ', num2str(k1), ' and ', num2str(gammatrim)])


% Initial conditions ------------------------------------------------------

out =  Trim_3D(Xtrim(1:14,ind)); %[pdot;qdot;rdot;udot;vdot;wdot;phidot;thetadot;(psidot-k1*Vtrim);(hdot-Vtrim*sin(gammatrim))...]; 
psidot = k1*Vtrim;
hdot  = Vtrim*sin(gammatrim);

uTrim     = Xtrim(11:14,ind);
uTrimLevel = Xtrim(11:14,171);

xOut(:,jj) = [Xtrim(1:9,ind); 0; 0; Xtrim(10,ind); uTrim; zeros(5,1)];
uOut(:,jj) = uTrim; 


xdot = StandardSimDynamics_CZ150(0,xOut(:,jj),uOut(:,jj),[0;0;0;1;0;0;0],param,zeros(6,1),0);  %[pdot;qdot;rdot;udot;vdot;wdot;phidot;thetadot;psidot;xFdot;yFdot;hdot;actdot;dryDyn]
err = [xdot(1:8)-out(1:8);xdot(9)-psidot-out(9); xdot(12)-hdot-out(10); norm(xdot(10:12))-Vtrim];


% Iteration ---------------------------------------------------------------
% disp('Running nonlinear simulation...')
for ii = 1:(pLen-1)
    jj = jj + 1; 

    odefun    = @(t,x) StandardSimDynamics_CZ150(t,x,uOut(:,jj-1),[0;0;0;1;0;0;0],param,zeros(6,1),0); 
        
    [~,xs]    = ode23(odefun,[0 timestep],xOut(:,jj-1));
    xOut(:,jj) = xs(end,:)';
    
        uOut(:,jj)   = uOut(:,jj-1);      % direct fin angle input
end
    
% results ------------------------------------
% figure
% plot(xOut(10,:),xOut(11,:)); 
% xlabel('X, forward (m)')
% ylabel('Y, side (m)')
% axis equal
% grid on;

% figure
% plot3(xOut(10,:),xOut(11,:),xOut(12,:))
% xlabel('X, forward (m)')
% ylabel('Y, side (m)')
% zlabel('Z (m)')
% axis equal
% grid on;

% figure
% subplot(3,1,1)
% plot(time, xOut(7,:)*180/pi)
% ylabel('\phi (deg)')
% subplot(3,1,2)
% plot(time, xOut(8,:)*180/pi)
% ylabel('\theta (deg)')
% subplot(3,1,3)
% plot(time, xOut(9,:)*180/pi)
% ylabel('\psi (deg)')
% xlabel('time (s)')
% 
% figure
% subplot(3,1,1)
% plot(time, xOut(4,:))
% ylabel('u (m/s)')
% subplot(3,1,2)
% plot(time, xOut(5,:))
% ylabel('v (m/s)')
% subplot(3,1,3)
% plot(time, xOut(6,:))
% ylabel('w (m/s)')
% xlabel('time (s)')
% 
% 
% figure
% subplot(3,1,1)
% plot(time, xOut(1,:)*180/pi)
% ylabel('p (deg/s)')
% subplot(3,1,2)
% plot(time, xOut(2,:)*180/pi)
% ylabel('q (deg/s)')
% subplot(3,1,3)
% plot(time, xOut(3,:)*180/pi)
% ylabel('r (deg/s)')
% xlabel('time (s)')
    primTemp = [time' xOut(10,:)' xOut(11,:)' xOut(12,:)' Vtrim*ones(pLen,1) xOut(9,:)' k1*ones(pLen,1) gammatrim*ones(pLen,1)]; % t x y z psi k gamma %%%%Fix -z***
    
    primEnd = primTemp(end,:);
    primLevelNew = primLevel;
    primLevelNew(:,1)=primLevelNew(:,1)+primEnd(1);
    Rpsi = [cos(primEnd(6)) -sin(primEnd(6));sin(primEnd(6)) cos(primEnd(6))];
    for ll = 1:size(primLevelNew,1)
        primLevelNew(ll,2:3)=(Rpsi*primLevelNew(ll,2:3)')';
    end
    primLevelNew(:,2:4)=primLevelNew(:,2:4)+primEnd(1,2:4);
    primLevelNew(:,6)=primLevelNew(:,6)+primEnd(1,6);
    
    % primArray{kk} = [primTemp;primLevelNew]; % t x y z Va psi k gamma
    primArray{kk} = [primTemp]; % t x y z Va psi k gamma

end

save primArray.mat primArray

rmpath(filepath)  