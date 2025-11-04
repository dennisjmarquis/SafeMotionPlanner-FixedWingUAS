function [pathData,pathK1Data,pathLenData,pathTotLen] = generateRefPath(trim, path, tf, dt)
    
    V_trim  = norm(trim(4:6,1));  % reference velocity       
       
    if ~isfield(path,'R')
        R = 100;
    else
        R = path.R;
    end
    
    if ~isfield(path, 'slip')
        slip = asin(trim(5,1)/V_trim);
    else
        slip = path.slip;
    end
    
    time    = 0:dt:tf;
    
      
    if strcmpi(path.type, 'circle') %============================================
        
        pathTotLen  = 2*pi*R;
        omega       = V_trim/R;  % reference turn rate
        
        for i = 1:numel(time)
            pathK1Data(1,i) = 1/R;                      % ref curvature
            pathData(1,i)   = omega*time(i);            % ref yaw angle
            pathData(2,i)   = R*sin(omega*time(i));     % ref X pos
            pathData(3,i)   = R - R*cos(omega*time(i)); % ref Y pos
            pathData(4,i)   = 0;                        % ref Z pos
            
            if i ~= 1
                pathLenData(1,i) = pathLenData(1,i-1) +...
                    norm(pathData(2:4,i)-pathData(2:4,i-1));
            else
                pathLenData(1,i) = 0;
            end
        end
       
    elseif strcmpi(path.type, 'straight') %============================================
        
        pathTotLen  = inf;
        
        for i = 1:numel(time)
            pathK1Data(1,i) = 0;                        % ref curvature
            pathData(1,i)   = 0;                        % ref yaw angle
            pathData(2,i)   = V_trim*time(i)*cos(slip); % ref X pos
            pathData(3,i)   = V_trim*time(i)*sin(slip); % ref Y pos
            pathData(4,i)   = 0;                        % ref Z pos
            
            if i ~= 1
                pathLenData(1,i) = pathLenData(1,i-1) +...
                    norm(pathData(2:4,i)-pathData(2:4,i-1));
            else
                pathLenData(1,i) = 0;
            end        
        end
        
        
    elseif strcmpi(path.type, 'lemniscate') %============================================
        
        pathTotLen  = 15.725*R;
        kmax_traj   = 1/R;
        omega       = 2*pi*V_trim/pathTotLen;
        
        N    = -3*sin(omega*time)/kmax_traj./(1+cos(omega*time).^2); 
        E    = -3*sin(omega*time).*cos(omega*time)/kmax_traj./(1+cos(omega*time).^2);
        Xref = -(N+E)/sqrt(2); 
        Yref = (N-E)/sqrt(2);
        Zref = zeros(1,size(Xref,2));
        
        Xdot  = gradient(Xref,dt); 
        Xddot = gradient(Xdot,dt); 
        Ydot  = gradient(Yref,dt);
        Yddot = gradient(Ydot,dt);

        Psiref = -atan2(Ydot,Xdot); % Psi goes negative, we take '-' to make psi +ve
        Psiref = wrapTo2Pi(Psiref); % then wrap it to [0,2pi]
        Psiref = -Psiref; % and change sign to wrap it to [-2pi,0]

        k1 = (-Ydot.*Xddot + Xdot.*Yddot)./(Xdot.^2+Ydot.^2).^1.5;
              
        pathK1Data = k1;                        % ref curvature
        pathData   = [Psiref;Xref;Yref;Zref];   % ref yaw angle

        for i = 1:numel(time)
            if i ~= 1
                pathLenData(1,i) = pathLenData(1,i-1) +...
                    norm(pathData(2:4,i)-pathData(2:4,i-1));
            else
                pathLenData(1,i) = 0;
            end  
        end
            
    else %============================================
        error('unknown path')        
    end
end

