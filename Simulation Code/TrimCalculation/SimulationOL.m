clear; clc; tic
filepath = '../Utility'; addpath(filepath)

global H0 Vtrim k1 gammatrim
H0 = 630;

load CircularTrim Xtrim
load param_CZ150.mat

d2r = pi/180; r2d = 1/d2r;  


% Simulation time ---------------------------------------------------------
timestep = 0.04;
time     = 0:timestep:20; %simulation time
pLen     = numel(time); %number of simultion
j        = 1;

indK1 = 11;
indGam = 16;
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

xOut(:,j) = [Xtrim(1:9,ind); 0; 0; Xtrim(10,ind); uTrim; zeros(5,1)];
uOut(:,j) = uTrim; 


xdot = StandardSimDynamics_CZ150(0,xOut(:,j),uOut(:,j),[0;0;0;1;0;0;0],param,zeros(6,1),0);  %[pdot;qdot;rdot;udot;vdot;wdot;phidot;thetadot;psidot;xFdot;yFdot;hdot;actdot;dryDyn]
err = [xdot(1:8)-out(1:8);xdot(9)-psidot-out(9); xdot(12)-hdot-out(10); norm(xdot(10:12))-Vtrim]


% Iteration ---------------------------------------------------------------
disp('Running nonlinear simulation...')
for i = 1:(pLen-1)
    j = j + 1; 

    odefun    = @(t,x) StandardSimDynamics_CZ150(t,x,uOut(:,j-1),[0;0;0;1;0;0;0],param,zeros(6,1),0); 
        
    [~,xs]    = ode23(odefun,[0 timestep],xOut(:,j-1));
    xOut(:,j) = xs(end,:)';
               
    uOut(:,j)   = uOut(:,j-1);      % direct fin angle input
end
    
% results ------------------------------------
figure
plot(xOut(10,:),xOut(11,:)); 
xlabel('X, forward (m)')
ylabel('Y, side (m)')
axis equal
grid on;

figure
plot3(xOut(10,:),xOut(11,:),xOut(12,:))
xlabel('X, forward (m)')
ylabel('Y, side (m)')
zlabel('Z (m)')
axis equal
grid on;

figure
subplot(3,1,1)
plot(time, xOut(7,:)*180/pi)
ylabel('\phi (deg)')
subplot(3,1,2)
plot(time, xOut(8,:)*180/pi)
ylabel('\theta (deg)')
subplot(3,1,3)
plot(time, xOut(9,:)*180/pi)
ylabel('\psi (deg)')
xlabel('time (s)')

figure
subplot(3,1,1)
plot(time, xOut(4,:))
ylabel('u (m/s)')
subplot(3,1,2)
plot(time, xOut(5,:))
ylabel('v (m/s)')
subplot(3,1,3)
plot(time, xOut(6,:))
ylabel('w (m/s)')
xlabel('time (s)')


figure
subplot(3,1,1)
plot(time, xOut(1,:)*180/pi)
ylabel('p (deg/s)')
subplot(3,1,2)
plot(time, xOut(2,:)*180/pi)
ylabel('q (deg/s)')
subplot(3,1,3)
plot(time, xOut(3,:)*180/pi)
ylabel('r (deg/s)')
xlabel('time (s)')

rmpath(filepath)  
toc