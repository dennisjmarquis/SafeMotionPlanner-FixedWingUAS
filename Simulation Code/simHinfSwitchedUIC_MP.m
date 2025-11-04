function [xOut, uOut, yOut, mu_r, dist, noise, reference_out,Xbvec] = simHinfSwitchedUIC_MP(time,x0,Xtrim,cont,simParam,param,reference,Xb0,primArray,trajec,leveltrim)
    global windBody
    global Rrod Robs Dsep Va Vb k1max g1max gammamax deltatTrajec tcritmax Xgoal F wvec flyMode tfin
    dt = mean(diff(time));

       %Extract Sim Params
    gustToggle = simParam.gustToggle;
    noiseToggle = simParam.noiseToggle;
    pNoiseToggle = simParam.pNoiseToggle;
    delayToggle = simParam.delayToggle;
    windToggle = simParam.windToggle;
    controlToggle = simParam.controlToggle;
    dn = simParam.dn;
    dn2 = simParam.dn2;
    pNoise = pNoiseToggle*simParam.pNoise;
    magCase = simParam.magCase; % gust magnitude identifier 1/2/3
    constWind = simParam.constWind; %NORTH/EAST/DOWN
    delaymin = simParam.delaymin;
    delaymax = simParam.delaymax;
    
    %Extract Controller
    Ak = cont.Ak;
    Bk = cont.Bk;
    Ck = cont.Ck;
    Dk = cont.Dk;

    N_UIC = length(Ak{1,1})-1;
    
    %Extract Reference
    Xref = reference.Xref; Yref = reference.Yref; Zref = reference.Zref; Varef = reference.Varef; Psiref = reference.Psiref; k1ref = reference.k1ref; gamref = reference.gamref;
    pref = reference.pref; qref = reference.qref; rref = reference.rref; phiref = reference.phiref; thetaref = reference.thetaref; dEref = reference.dEref; dAref = reference.dAref; dRref = reference.dRref; dTref = reference.dTref;
    k1val = reference.k1val; gmval = reference.gmval;
    k1valcont = cont.k1val; gmvalcont = cont.gmval;
    
    %Get level trim commands
    N_gm = length(gmval);
    ind_k1_level = 6;
    ind_gam_level = 16;
    ind_level = (ind_k1_level-1)*N_gm+ind_gam_level; %171: straight and level
    levelTrim = Xtrim(:,ind_level);
    uTrim     = levelTrim(11:14);

    % initial conditions -----------------------------------------------
    xOut(:,1) = x0; % p q r u v w phi tht psi x y h xe1 xe2 xa1 ... dryden
    uOut(:,1) = uTrim; % dE dA dR dT
    
    uDelay    = uOut(:,1);
    ubar(:,1) = zeros(4,1);
    
    measNoise   = simParam.noiseToggle*[0.01*ones(1,3) 2 0.01*ones(1,3) 2*ones(1,3) ]; %standard deviation
    
    wind(:,1)   = [gustToggle*dn(1:3,1); magCase; windToggle*constWind]; wind(:,2)=wind;
    dist(:,1)   = [windToggle*constWind; dn(4:end,1)];
    noise(:,1)  = dn(4:end,1).*measNoise';
    yOut(:,1) = StandardSimObservations_CZ150(x0,wind,noise(:,1));
    yTrim     = StandardSimObservations_CZ150(x0,[zeros(3,1);1;zeros(3,1)], zeros(10,1)); % trim measurement
    
    xc       = zeros(size(Ak{1,1}{1},1),1);

    ind_k1_old = 0;
    ind_gam_old = 0;


    %% iteration
    pLen = length(time);
    jj = 1;

    Xbvec = zeros(5,pLen);
    Xbvec(:,1) = Xb0;
    
    for ii = 1:(pLen-1)
        jj = jj + 1;

        Xa=[Xref(ii);Yref(ii);Zref(ii);Psiref(ii);gamref(ii)]; %x y z psi gamma
        Xb=Xbvec(:,ii); %x y z psi k gamma
    
        %Plan trajectory with motion planner
        [trajec,flyMode,flag_wc,tfin] = motionPlanner3D(Xa,Xb,trajec,time(ii),time(end),flyMode,tfin,primArray,leveltrim);
        Xref   = [Xref(1:ii); trajec(2,:)'];
        Yref   = [Yref(1:ii); trajec(3,:)'];
        Zref   = [Zref(1:ii); trajec(4,:)'];
        Varef = [Varef(1:ii); trajec(5,:)'];
        Psiref = [Psiref(1:ii); trajec(6,:)'];
        k1ref  = [k1ref(1:ii); trajec(7,:)'];
        gamref = [gamref(1:ii); trajec(8,:)'];
        pref = [pref(1:ii); trajec(9,:)'];
        qref = [qref(1:ii); trajec(10,:)'];
        rref = [rref(1:ii); trajec(11,:)'];
        phiref = [phiref(1:ii); trajec(12,:)'];
        thetaref = [thetaref(1:ii); trajec(13,:)'];
        dEref = [dEref(1:ii); trajec(14,:)'];
        dAref = [dAref(1:ii); trajec(15,:)'];
        dRref = [dRref(1:ii); trajec(16,:)'];
        dTref = [dTref(1:ii); trajec(17,:)'];
    


        %Switch Controller and Assign Trims
        ind_k1_con = find(abs(k1valcont-k1ref(jj))==min(abs(k1valcont-k1ref(jj))),1);
        ind_gam_con = find(abs(gmvalcont-gamref(jj))==min(abs(gmvalcont-gamref(jj))),1);

        ind_k1 = find(abs(k1val-k1ref(jj))==min(abs(k1val-k1ref(jj))),1);
        ind_gam = find(abs(gmval-gamref(jj))==min(abs(gmval-gamref(jj))),1);

        if ind_k1_con ~= ind_k1_old || ind_gam_con ~= ind_gam_old
            % fprintf('Switching to controller where k = %f and gamma = %f at t = %f\n',k1valcont(ind_k1_con),gmvalcont(ind_gam_con),time(ii))
            ind_UIC = 1;
        end    

        Acon = Ak{ind_k1_con,ind_gam_con}{ind_UIC};
        Bcon = Bk{ind_k1_con,ind_gam_con}{ind_UIC};
        Ccon = Ck{ind_k1_con,ind_gam_con}{ind_UIC};
        Dcon = Dk{ind_k1_con,ind_gam_con}{ind_UIC};
        % fprintf('for time step %f, controller (%d,%d) is being used with index %d\n',ii,ind_k1,ind_gam,ind_UIC)
        ind = (ind_k1-1)*N_gm+ind_gam; %171: straight and level
        trim  = Xtrim(:,ind);
        uTrim     = trim(11:14);
        xOutTrim  = [trim(1:9); 0; 0; trim(10,1); uTrim; zeros(5,1)];
        yTrim     = StandardSimObservations_CZ150(xOutTrim,[zeros(3,1);1;zeros(3,1)], zeros(10,1));

        wind(1:3,jj) = gustToggle*dn(1:3,jj); % dryden model input
    
        tdel = rand(1)*(delaymax-delaymin)+delaymin;      
        if delayToggle==0
            tdel = 0.000001;
        end
        delayvec = [0 tdel];
    
        % state and measurement update -------------------------
        odefundelay = @(t,x) StandardSimDynamics_CZ150(t,x,real(uDelay),wind(:,jj-1),param,pNoise(:,jj-1),0);        
        [~,xd]      = ode23(odefundelay,delayvec,xOut(:,jj-1));
        xDelay      = xd(end,:)';
        dtvec       = [0 dt-tdel];
    
        odefun    = @(t,x) StandardSimDynamics_CZ150(t,x,real(uOut(:,jj-1)),wind(:,jj-1),param,pNoise(:,jj-1),0);      
        [~,xs]    = ode23(odefun,dtvec,xDelay);
        xOut(:,jj) = xs(end,:)';
        
        if abs(xOut(9,jj))>2*pi
            xOut(9,jj) = xOut(9,jj)-2*pi*sign(xOut(9,jj));
        end
        
        dist(:,jj)  = [windBody; dn(4:end,jj)];
        noise(:,jj) = dn(4:end,jj).*measNoise';
        yOut(:,jj)  = StandardSimObservations_CZ150(xOut(:,jj),wind(:,jj),noise(:,jj)); 
    
    
        % error computation ---------------------------------
        if abs(Psiref(jj,1))>2*pi
            Psiref(jj,1) = Psiref(jj,1)-2*pi*sign(Psiref(jj,1));
        end
        Xe(jj,1)     = yOut(8,jj)-Xref(jj,1);
        Ye(jj,1)     = yOut(9,jj)-Yref(jj,1); 
        He(jj,1)     = yOut(10,jj) + Zref(jj,1); % Confirm sign!!!
        psie(jj,1)   = wrapToPi(yOut(7,jj)-Psiref(jj,1));
    
        dyx = Xe(jj,1); dyy = Ye(jj,1);
        Xe(jj,1) = sin(yOut(7,jj))*dyy + cos(yOut(7,jj))*dyx;
        Ye(jj,1) = cos(yOut(7,jj))*dyy - sin(yOut(7,jj))*dyx;
        
        measRef = [pref(jj); qref(jj); rref(jj); Varef(jj); phiref(jj); thetaref(jj)]; 
        dy(1:6,1)  = yOut(1:6,jj) - measRef; % del(p q r V phi theta)
        dy(7:10,1) = [psie(jj,1); Xe(jj,1); Ye(jj,1); He(jj,1)];    
        
            
        % control and update --------------------------------------------------
        uRef = [dEref(jj);dAref(jj);dRref(jj);dTref(jj)];

       if controlToggle == 1  
            ubar(:,jj) = Ccon*xc + Dcon*dy;
        end
    
        uOut(:,jj)   = ubar(:,jj)+uRef;
        xc          = Acon*xc + Bcon*dy;
        wind(:,jj+1) = wind(:,jj);
        uDelay      = uOut(:,jj-1);

        ind_k1_old = ind_k1_con;
        ind_gam_old = ind_gam_con;

        ind_UIC = ind_UIC+1;
        if ind_UIC > N_UIC + 1
            ind_UIC = N_UIC +1;
        end

        %Propagate obstacle
        Xb1=[Xb(1:3)+Vb*dt*[cos(Xb(4))*cos(Xb(5));sin(Xb(4))*cos(Xb(5));sin(Xb(5))];Xb(4);Xb(5)];
        Xbvec(:,jj)=Xb1;
    end
        
     
    %% Results
    mu_r           = (sum(power(ubar,2), 'all'))^0.5;

    reference_out = reference;
    reference_out.Xref = Xref; reference_out.Yref = Yref; reference_out.Zref = Zref; reference_out.Varef = Varef;
    reference_out.Psiref = Psiref; reference_out.k1ref = k1ref; reference_out.gamref = gamref;
    reference_out.pref = pref; reference_out.qref = qref; reference_out.rref = rref; reference_out.phiref = phiref; reference_out.thetaref = thetaref;
    reference_out.dEref = dEref; reference_out.dAref = dAref; reference_out.dRref = dRref; reference_out.dTref = dTref;

end