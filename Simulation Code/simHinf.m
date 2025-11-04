function [xOut, uOut, yOut, dy, xc, mu_r, dist, noise, noise2, pathError,yOutEKF] = simHinf(time,x0,Trim,cont,simParam,param,reference)
    global windBody
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
    Acon = cont.Ak;
    Bcon = cont.Bk;
    Ccon = cont.Ck;
    Dcon = cont.Dk;
    
    %Extract Reference
     Xref = reference.Xref; Yref = reference.Yref; Zref = reference.Zref;
    Psiref = reference.Psiref; k1ref = reference.k1ref; gamref = reference.gamref;
    noPosFlag = reference.noPosFlag;
    
    uTrim     = Trim(11:14);

    % initial conditions -----------------------------------------------
    xOut(:,1) = x0; % p q r u v w phi tht psi x y h dE dA dR dT dryden
    uOut(:,1) = x0(13:16); % dE dA dR dT
    
    uDelay    = uOut(:,1);
    ubar(:,1) = zeros(4,1);
    
    measNoise   = noiseToggle*[0.01*ones(1,3) 2 0.1*ones(1,3) 0.1*ones(1,3) ]; %standard deviation
    measNoise2 = noiseToggle*[0.1*ones(1,3) 0.3*ones(1,3)];

    wind(:,1)   = [gustToggle*dn(1:3,1); magCase; windToggle*constWind]; wind(:,2)=wind;
    dist(:,1)   = [windToggle*constWind; dn(4:end,1)];
    noise(:,1)  = dn(4:end,1).*measNoise';
    noise2(:,1) = dn2(:,1).*measNoise2';
    yOut(:,1) = StandardSimObservations_CZ150(x0,wind,noise(:,1));
    [~,ax,ay,az] = StandardSimDynamics_CZ150(0,xOut(:,1),uOut(:,1),wind(:,1),param,pNoise(:,1),0);
    yOutEKF(:,1) = [yOut(1:3,1);xOut(4:6,1)+noise2(1:3,1);yOut(5:7,1);yOut(4,1);[ax;ay;az]+noise2(4:6,1)];
    xTrim = [Trim(1:8);0;0;0;0; uTrim; zeros(5,1)];
    yTrim     = StandardSimObservations_CZ150(xTrim,[zeros(3,1);1;zeros(3,1)], zeros(10,1)); % trim measurement
    
    Xe(1,1)     = yOut(8,1)-Xref(1,1);
    Ye(1,1)     = yOut(9,1)-Yref(1,1); 
    He(1,1)     = yOut(10,1) - (-Zref(1,1)); % yOut(10) = H 
    psie(1,1)   = wrapToPi(yOut(7,1)-Psiref(1,1));
    
    dyx = Xe(1,1); dyy = Ye(1,1);
    Xe(1,1) = sin(yOut(7,1))*dyy + cos(yOut(7,1))*dyx;
    Ye(1,1) = cos(yOut(7,1))*dyy - sin(yOut(7,1))*dyx;
    dy(1:6,1)  = yOut(1:6,1) - yTrim(1:6); % del(p q r V phi theta)
    dy(7:10,1) = [psie(1,1); Xe(1,1); Ye(1,1); He(1,1)];
    
    xc(:,1)       = zeros(size(Acon,1),1);
    
    %% iteration
    pLen = length(time);
    jj = 1;
    for ii = 1:(pLen-1)
        jj = jj + 1; 
    
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

    
        xOut(9,jj) = wrapToPi(xOut(9,jj));
        
        dist(:,jj)  = [windBody; dn(4:end,jj)];
        noise(:,jj) = dn(4:end,jj).*measNoise';
        noise2(:,jj) = dn2(:,jj).*measNoise2';
        yOut(:,jj)  = StandardSimObservations_CZ150(xOut(:,jj),wind(:,jj),noise(:,jj)); 
        [~,ax,ay,az] = StandardSimDynamics_CZ150(dt-tdel,xOut(:,jj),uOut(:,jj-1),wind(:,jj-1),param,pNoise(:,jj-1),0);
        yOutEKF(:,jj) = [yOut(1:3,jj);xOut(4:6,jj)+noise2(1:3,jj);yOut(5:7,jj);yOut(4,jj);[ax;ay;az]+noise2(4:6,jj)];
    
        % error computation ---------------------------------
        Psiref(jj,1) = wrapToPi(Psiref(jj,1));
        
        Xe(jj,1)     = yOut(8,jj)-Xref(jj,1);
        Ye(jj,1)     = yOut(9,jj)-Yref(jj,1); 
        He(jj,1)     = yOut(10,jj) - (-Zref(jj,1)); % yOut(10) = H 
        psie(jj,1)   = wrapToPi(yOut(7,jj)-Psiref(jj,1));
    
        dyx = Xe(jj,1); dyy = Ye(jj,1);
        Xe(jj,1) = sin(yOut(7,jj))*dyy + cos(yOut(7,jj))*dyx;
        Ye(jj,1) = cos(yOut(7,jj))*dyy - sin(yOut(7,jj))*dyx;
        
        if noPosFlag == 1
            Xe(jj,1) = 0; Ye(jj,1) = 0; He(jj,1) = 0;
        end

        dy(1:6,jj)  = yOut(1:6,jj) - yTrim(1:6); % del(p q r V phi theta)
        dy(7:10,jj) = [psie(jj,1); Xe(jj,1); Ye(jj,1); He(jj,1)];    
        
            
        % control and update --------------------------------------------------
        uOut(:,jj) = uTrim;

       if controlToggle == 1  
            uOut(:,jj) = Ccon*xc(:,jj-1) + Dcon*dy(:,jj) + uOut(:,jj);
       else
           uOut(:,jj) = simParam.uout(jj,:)';
       end
    
        ubar(:,jj)   = uOut(:,jj)-uTrim;
        xc(:,jj)          = Acon*xc(:,jj-1) + Bcon*dy(:,jj);
        wind(:,jj+1) = wind(:,jj);
        uDelay      = uOut(:,jj-1);
    end
        
     
    %% Results
    mu_r           = (sum(power(ubar,2), 'all'))^0.5;
    curve          = [Xref Yref -Zref];
    [~,pathError,~] = distance2curve(curve, xOut(10:12,:)');

end