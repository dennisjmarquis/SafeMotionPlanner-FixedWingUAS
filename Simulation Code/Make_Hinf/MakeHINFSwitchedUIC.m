clear; close all; clc

filepath = '../Utility'; addpath(filepath)

global H0; H0 = 630;
load param_CZ150.mat; % aerodynamic coefficients
load ../TrimCalculation/CircularTrim

%N-eventually time-invariant controller
N_UIC = 5;

%% TRIM POINT:

k1list = [1 3 6 9 11];     % 1/ROC
gmlist = [1 9 16 23 31]; % FPA


N_k1 = length(k1val);
N_gm = length(gmval);

for ii = 1:length(k1list)
    for jj = 1:length(gmlist)
    index_k1 = k1list(ii); 
    index_gm = gmlist(jj);
    
    index = (index_k1-1)*N_gm+index_gm;
    fprintf('%d\n',index)
    
    x = Xtrim(:,index);
    
    % Model (1: only aircraft dyn, 2: + actuator dyn, 3: + Dryden dyn)
    modelFlag = 2; 
    
    % performance weight
    %     pPen_p = 0.191; pPen_dA = 2.0;
    %     qPen_q = 0.191; qPen_dE = 0.1;
    %     rPen_r = 0.191; rPen_dR = 0.1;
    %     VPen_V = 0.1; VPen_dT = 1;
    %     phiPen = 0.191;
    %     thetaPen = 0.191;
    %     psiPen = 0.191;
    %     hPen = 0.03;
    %     xFPen = 0.03;
    %     yFPen = 0.03;
    %     dEPen = 0.4;
    %     dAPen = 0.1;
    %     dRPen = 1.0;
    %     dTPen = 2.0;
    % Z_wt = [pPen_p;pPen_dA;qPen_q;qPen_dE;rPen_r;rPen_dR;VPen_V;VPen_dT;phiPen;thetaPen;psiPen;hPen;xFPen;yFPen;dEPen;dAPen;dRPen;dTPen];
    % Z_names = {'p_p','p_dA','q_q','q_dE','r_r','r_dR','Va_Va','Va_dT','phi','theta','psi','h','x','y','dE','dA','dR','dT'};
    [Z_wt, Z_names] = PerfWt(index_k1,index_gm,modelFlag);
    

    p = x(1); q = x(2); r = x(3); u = x(4); v = x(5); w = x(6); phi = x(7);
    theta = x(8); psi = x(9); h = x(10); dE = x(11); dA = x(12); dR = x(13);
    dT = x(14); 
    clear x;
    
    
    if modelFlag == 1
        nx = 12; d1 = 3; d2 = 10; nd = d1 + d2; nu = 4; ny = 10; nz = 14;
        xTrim = [p; q; r; u; v; w; phi; theta; psi; 0; 0; h];
        uTrim = [dE; dA; dR; dT]; 
    elseif modelFlag == 2
        nx = 16; d1 = 3; d2 = 10; nd = d1 + d2; nu = 4; ny = 10; nz = 14;
        xTrim = [p; q; r; u; v; w; phi; theta; psi; 0; 0; h; dE; dA; dR; dT];
        uTrim = [dE; dA; dR; dT]; 
    elseif modelFlag == 3
        nx = 21; d1 = 3; d2 = 10; d3 = 2; nd = d1 + d2 + d3; nu = 4; ny = 10; nz = 14;
        xTrim = [p; q; r; u; v; w; phi; theta; psi; 0; 0; h; dE; dA; dR; dT; zeros(5,1)];
        uTrim = [dE; dA; dR; dT]; 
    end
    
     %Uncertain Initial Condition %using weights in r, phi, theta
    Lambda = zeros(nx,5);
    Lambda(3,1) = 0.1; %r uncertainty
    Lambda(7,2) = 0.5; %phi uncertainty
    Lambda(8,3) = 0.5; %theta uncertainty
    Lambda(13,4) = 0.2; %dE uncertainty
    Lambda(16,5) = 0.2; %dT uncertainty

    %% LINEARIZE EQUATIONS OF MOTION:
    
    fprintf('Linearizing equations of motion for k = %f and gamma = %f \n',k1val(index_k1),gmval(index_gm))
    
    %Initialize matrix dimensions
    % (nn = states, d = exogenous inputs, m = control inputs)
    % (q = performance vector, p = output)
    
    A = zeros(nx,nx);    B1 = zeros(nx,nd);   B2 = zeros(nx,nu);
    C1 = zeros(nz,nx);    D11 = zeros(nz,nd);   D12 = zeros(nz,nu);
    C2 = zeros(ny,nx);    D21 = zeros(ny,nd);   D22 = zeros(ny,nu);
    dTrim = zeros(nd,1);
    e = 1e-5;
    dt = 0.04;
    % A matrix
    for j = 1:nx
        dx = zeros(nx,1);
        dx(j) = e;
        dxp = StandardDynamics(xTrim+dx,uTrim,dTrim,param, 0, modelFlag);
        dxm = StandardDynamics(xTrim-dx,uTrim,dTrim,param, 0, modelFlag);
        A(:,j) = (dxp-dxm)/2/e;
    end
    % B1 matrix
    for j = 1:nd
        dd = zeros(nd,1);
        dd(j) = e;
        dxp =  StandardDynamics(xTrim,uTrim,dTrim+dd,param, 0, modelFlag);
        dxm =  StandardDynamics(xTrim,uTrim,dTrim-dd,param, 0, modelFlag);
        B1(:,j) = (dxp-dxm)/2/e;
    end
    % B2 matrix
    for j = 1:3
        du = zeros(nu,1);
        du(j) = e;
        dxp =  StandardDynamics(xTrim,uTrim+du,dTrim,param, 0, modelFlag);
        dxm =  StandardDynamics(xTrim,uTrim-du,dTrim,param, 0, modelFlag);
        B2(:,j) = (dxp-dxm)/2/e;
    end
    for j = 4
        du = zeros(nu,1);
        du(j) = e*900;
        dxp =  StandardDynamics(xTrim,uTrim+du,dTrim,param, 0, modelFlag);
        dxm =  StandardDynamics(xTrim,uTrim-du,dTrim,param, 0, modelFlag);
        B2(:,j) = (dxp-dxm)/2/e/900;
    end
    % C1 matrix
    for j = 1:nx
        dx = zeros(nx,1);
        dx(j) = e;
        dzp = RightPerformance(xTrim+dx,uTrim,dTrim,Z_wt, modelFlag);
        dzm = RightPerformance(xTrim-dx,uTrim,dTrim,Z_wt, modelFlag);
        C1(:,j) = (dzp-dzm)/2/e;
    end
    % D11 matrix
    for j = 1:nd
        dd = zeros(nd,1);
        dd(j) = e;
        dzp = RightPerformance(xTrim,uTrim,dTrim+dd,Z_wt, modelFlag);
        dzm = RightPerformance(xTrim,uTrim,dTrim-dd,Z_wt, modelFlag);
        D11(:,j) = (dzp-dzm)/2/e;
    end
    % D12 matrix
    for j = 1:nu
        du = zeros(nu,1);
        du(j) = e;
        dzp = RightPerformance(xTrim,uTrim+du,dTrim,Z_wt, modelFlag);
        dzm = RightPerformance(xTrim,uTrim-du,dTrim,Z_wt, modelFlag);
        D12(:,j) = (dzp-dzm)/2/e;
    end
    % C2 matrix
    for j = 1:nx
        dx = zeros(nx,1);
        dx(j) = e;
        dyp = StandardObservations(xTrim+dx,uTrim,dTrim, modelFlag);
        dym = StandardObservations(xTrim-dx,uTrim,dTrim, modelFlag);
        C2(:,j) = (dyp-dym)/2/e;
    end
    % D21 matrix
    for j = 1:nd
        dd = zeros(nd,1);
        dd(j) = e;
        dzp = StandardObservations(xTrim,uTrim,dTrim+dd, modelFlag);
        dzm = StandardObservations(xTrim,uTrim,dTrim-dd, modelFlag);
        D21(:,j) = (dzp-dzm)/2/e;
    end
    
    D22 = zeros(ny,nu);
    % D21 = [zeros(10,3) diag([0.01*ones(1,3) 2 0.01*ones(1,3) 2*ones(1,3)])];
    
    %disturbance weight matrix
    sf = 3; %worst-case scale factor
    windsf = 10;
    wwind = windsf*eye(3);
    wnoise = diag([sf*0.01*ones(1,3) sf*2 sf*0.01*ones(1,3) sf*2*ones(1,3)]);
    WD = blkdiag(wwind,wnoise); %13*13 weight matrix
    if modelFlag == 3
        WD = blkdiag(WD,sf*eye(2)); %15*15 weight matrix
    end
    B1 = B1*WD;
    D11 = D11*WD;
    D21 = D21*WD;
    
    % Discretization
    csys = ss(A,[B1 B2],[C1;C2],[D11 D12; D21 D22]);
    dsys = c2d(csys,dt,'zoh');
    
    A = dsys.a;
    B1 = dsys.b(:,1:nd);
    B2 = dsys.b(:,nd+1:end);
    C1 = dsys.c(1:nz,:);
    C2 = dsys.c(nz+1:end,:);
    D11 = dsys.d(1:nz,1:nd);
    D12 = dsys.d(1:nz,nd+1:end);
    D21 = dsys.d(nz+1:end,1:nd);
    D22 = dsys.d(nz+1:end,nd+1:end);
    
    %% CHECK STABILIZABILITY
    unstabilizable = 0;
    evals = eig(A);
    n = length(evals);
    for j = 1:n
        test = rank([evals(j)*eye(n)-A, B2]);
        if test < n
            if abs(evals(j)) > 1
                fprintf('unstabilizable: %2.2f \n',evals(j))
                unstabilizable = 1;
            end
        end
    end
    if unstabilizable == 0
        fprintf('All uncontrollable modes are stable! \n')
    end
    
    
    %% MAKE H-INFINITY UIC CONTROLLER
    gam = [];
    tic
    [~,~,~,gamOpt] = LTIoptsolUIC(A, B1, B2, C1, C2, D11, D12, D21, Lambda, N_UIC, gam);
    
    ttt=toc;
    fprintf('Optimal performance value: %d\n',gamOpt);
    
    %% relaxing and control matrices synthesis
    gam = 1.3*gamOpt;
    [R,S,scal,gam] = LTIoptsolUIC(A,B1,B2,C1,C2,D11,D12,D21, Lambda, N_UIC, gam);
    fprintf('Relaxed performance value: %d\n',gam);
    
    [Ak, Bk, Ck, Dk] = LTIcontrollerUIC(R, S, scal, A, B1, B2, C1, C2, D11, D12, D21, gam);
    
    %% NOMINAL CLOSED LOOP SYSTEM
    flag = 0;
    for kk = 1:N_UIC+1
        [Acl, Bcl, Ccl, Dcl] =  CloseWController(A, B1, B2, C1, C2, D11, D12, D21, Ak{kk}, Bk{kk}, Ck{kk}, Dk{kk});
        if max(abs(eig(Acl))) >= 1
            flag = 1;
        end
    end
    if flag == 0
        fprintf('Closed loop eigenvalues are stable!\n\n')
    else
        fprintf('***WARNING*** Closed loop eigenvalues are not stable\n\n')
    end
    
    contHinfSwitchedUIC.flag{ii,jj} = flag;
    contHinfSwitchedUIC.Ak{ii,jj} = Ak;
    contHinfSwitchedUIC.Bk{ii,jj} = Bk;
    contHinfSwitchedUIC.Ck{ii,jj} = Ck;
    contHinfSwitchedUIC.Dk{ii,jj} = Dk;
    contHinfSwitchedUIC.gam{ii,jj} = gam;
    contHinfSwitchedUIC.gamOpt{ii,jj} = gamOpt;
    contHinfSwitchedUIC.modelFlag{ii,jj} = modelFlag;
    contHinfSwitchedUIC.Z_wt = Z_wt;
    contHinfSwitchedUIC.Z_names = Z_names;
    contHinfSwitchedUIC.k1{ii,jj} = k1val(k1list(ii));
    contHinfSwitchedUIC.gamma{ii,jj} = gmval(gmlist(jj));
    contHinfSwitchedUIC.dt = dt;
    contHinfSwitchedUIC.k1val = k1val(k1list);
    contHinfSwitchedUIC.gmval = gmval(gmlist);

    sysHinfSwitchedUIC.A{ii,jj} = A;
    sysHinfSwitchedUIC.B1{ii,jj} = B1;
    sysHinfSwitchedUIC.B2{ii,jj} = B2;
    sysHinfSwitchedUIC.C1{ii,jj} = C1;
    sysHinfSwitchedUIC.C2{ii,jj} = C2;
    sysHinfSwitchedUIC.D11{ii,jj} = D11;
    sysHinfSwitchedUIC.D12{ii,jj} = D12;
    sysHinfSwitchedUIC.D21{ii,jj} = D21;
    sysHinfSwitchedUIC.D22{ii,jj} = D22;
    

    end
end

contHinf = contHinfSwitchedUIC;
save('cont_RARL_HINF_sw_UIC2.mat','contHinf')
% save('sys_HINF_sw_UIC.mat','sysHinfSwitchedUIC')

rmpath(filepath)

%% supporting functions

function [Z_wt, Z_names] = PerfWt(index_k1,index_gm,modelFlag)

    % performance weight
    % if index_k1 == 6 && index_gm == 16
    %     % level
    %     pPen = 0.1;
    %     qPen = 0.1;
    %     rPen = 0.1;
    %     VPen = 0.06; 
    %     phiPen = 0.5;
    %     thetaPen = 0.5;
    %     psiPen = 0.4;
    %     hPen = 0.12; 
    %     xFPen = 0.12; 
    %     yFPen = 0.12; 
    %     dEPen = 0.4;
    %     dAPen = 0.4;
    %     dRPen = 0.6;
    %     dTPen = 1.0;
    % else
        % circular, weights used in controller 13 from Hinf test
        pPen = 1/(2*pi);
        qPen = 1/(2*pi);
        rPen = 1/(2*pi);
        VPen = 0.1;
        phiPen = 4/pi;
        thetaPen = 4/pi;
        psiPen = 1/pi;
        hPen = 0.2;
        xFPen = 0.2;
        yFPen = 0.2;
        dEPen = 5;
        dAPen = 5;
        dRPen = 5;
        dTPen = 5;
    % end
    Z_wt = [pPen;qPen;rPen;VPen;phiPen;thetaPen;psiPen;hPen;xFPen;yFPen;dEPen;dAPen;dRPen;dTPen];
    Z_names = {'p','q','r','Va','phi','theta','psi','h','x','y','dE','dA','dR','dT'};
   
end