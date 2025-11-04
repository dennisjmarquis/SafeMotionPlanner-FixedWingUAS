function [xk1_p, xk1_p_bar, bk1_p, Vk1, Pxk1_p_bar, Pbk1_p_bar, Pxk1_p, Pxk1_m, etaxM, etabM,lambdax, lambdab] = myATEKF(xk_p, xk_p_bar, bk_p, inputk, inputk1, zk1, Vk, Pxk_p_bar, Pbk_p_bar, etaxM, etabM, lambdax, lambdab, Qx, Qb, R, param,dt)
%Implemented from Kim 2009
nx = length(xk_p); nb = 3; nm = 16;
M = size(etaxM,2);
Ak = eye(nb);

%Get Fk, Hk1
[Fk,Bk] = discretizedJacobians(@f_CZ150EKF,0,xk_p,inputk,bk_p,param,dt,1);


[~, xs1] = rungeKutta4(@f_CZ150EKF,[0 dt], xk_p, inputk, bk_p, param, 10);
f = xs1(:,end);

Pbk1_m_bar = lambdab*(Ak*Pbk_p_bar*Ak' + Qb); %19b

Uk1_bar = (Fk*Vk + Bk)*inv(Ak); %20c
Uk1 = Uk1_bar*(eye(nb)-lambdab*Qb*inv(Pbk1_m_bar)); %20b
uk = (Uk1_bar - Uk1)*Ak*bk_p; %20e
Qxk_bar = Qx + Uk1*Qb*Uk1_bar'; %20f

Ck = f - Fk*xk_p -Bk*bk_p;
xk1_m_bar = Fk*xk_p_bar + Ck + uk; %18a

bk1_m = Ak*bk_p; %19a
xk1_m = xk1_m_bar + Uk1*bk1_m; %17a

[Hk1,Dk1]=discretizedJacobians(@h_CZ150EKF,0,xk1_m,inputk1,bk1_m,param,dt,2);

h = h_CZ150EKF(xk1_m,inputk1,bk1_m,param);

Pxk1_m_bar = lambdax*(Fk*Pxk_p_bar*Fk' + Qxk_bar); %18b
Kxk1 = Pxk1_m_bar*Hk1'*inv(Hk1*Pxk1_m_bar*Hk1' + R); %18c
Pxk1_p_bar = (eye(nx)-Kxk1*Hk1)*Pxk1_m_bar; %18d

Ek1 = h - Hk1*xk1_m - Dk1*bk1_m;

innov = zk1 - Hk1*xk1_m_bar;
innov(12)=wrapToPi(innov(12));

etaxk1 = innov - Ek1; %18e

etaxM = etaxM(:,2:end); etaxM = [etaxM etaxk1];
xk1_p_bar = xk1_m_bar + Kxk1*etaxk1; %18f
Cxk1 = Hk1*Pxk1_m_bar*Hk1' + R; %18g
Cxk1_bar = 1/(M-1)*(etaxM*etaxM'); %18h
lambdax = max(1,trace(Cxk1_bar)/trace(Cxk1));

Nk1 = Hk1*Uk1+Dk1; %20a

Kbk1 = Pbk1_m_bar*Nk1'*inv(Hk1*Pxk1_m_bar*Hk1' + R + Nk1*Pbk1_m_bar*Nk1'); %19c
Pbk1_p_bar = (eye(nb) - Kbk1*Nk1)*Pbk1_m_bar; %19d
etabk1 = etaxk1 - Nk1*bk1_m; %19e
etabM = etabM(:,2:end); etabM = [etabM etabk1];
bk1_p = bk1_m + Kbk1*etabk1; %19f
Cbk1 = Hk1*Pxk1_m_bar*Hk1' + R + Nk1*Pbk1_m_bar*Nk1'; %19g
Cbk1_bar = 1/(M-1)*(etabM*etabM'); %19h
lambdab = max(1,trace(Cbk1_bar)/trace(Cbk1));

Vk1 = Uk1 - Kxk1*Nk1; %20d

xk1_p = xk1_p_bar + Vk1*bk1_p; %17b
Pxk1_m = Pxk1_m_bar + Uk1*Pbk1_m_bar*Uk1'; %17c
Pxk1_p = Pxk1_p_bar + Vk1*Pbk1_p_bar*Vk1'; %17d
end