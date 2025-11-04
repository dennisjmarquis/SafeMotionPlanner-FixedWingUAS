function [newtrajec,flyMode,flag_wc,tfin] = motionPlanner3D(Xa,Xb,trajec,t,tendsim,flyMode,tfin,primArray,leveltrim)
%%inputs 
 %Xa : current state of ownship
 %Xb : current state of obstacle
 %trajec : current planned trajectory
%%outputs 
 %trajec : time parameterized trajectory matrix [t;x;y;psi]
 
 global Rrod Robs Dsep Va Vb kmax deltatTrajec tcritmax Xgoal F wvec

 Nev = length(primArray);

  %Get relative position and velocity (X is x y z psi gamma)
 Rrel = Xb(1:3)-Xa(1:3);
 Vrel = [Vb*cos(Xb(4))*cos(Xb(5))-Va*cos(Xa(4))*cos(Xa(5)); Vb*sin(Xb(4))*cos(Xb(5))-Va*sin(Xa(4))*cos(Xa(5));Vb*sin(Xb(5))-Va*sin(Xa(5))];
  

 %Check for obstacle in ROD
 if norm(Rrel)<=(Rrod+Robs)
     %Check for well-clear violation
    [flag_wc, tcrit,tca] = checkWellClear(Rrel,Vrel);

     if flag_wc == true && ~strcmp(flyMode,'evade') && tcrit <= tcritmax
         Jmin = inf;
%          fprintf('Avoidance maneuver attempted at t = %f\n',t)
         for ii = 1:Nev
             prim = primArray{ii};
                    Ra_evade = Xa(1:3)+ prim(end,2:4)';
                    psia_evade = Xa(4) + prim(end,6);
                    Xa_evade = [Ra_evade;psia_evade;0];
                    tevade = prim(end,1);
                    kevade = prim(1,7);
                    gamevade = prim(1,8);
                   
                    Rb_evade = Xb(1:3)+[Vb*tevade*cos(Xb(4))*cos(Xb(5));Vb*tevade*sin(Xb(4))*cos(Xb(5));Vb*tevade*sin(Xb(5))];
                    Rrel_evade = Rb_evade-Ra_evade;
                    Vrel_evade = [Vb*cos(Xb(4))*cos(Xb(5))-Va*cos(psia_evade);Vb*sin(Xb(4))*cos(Xb(5))-Va*sin(psia_evade);Vb*sin(Xb(5))];
                    flag_wc_evade = checkWellClear(Rrel_evade,Vrel_evade);
                 if flag_wc_evade == false && tevade < tcrit
                     J = getJMotionPlanner(Xa,Xa_evade,Xgoal,F,wvec); %compute new J
                     if J < Jmin
                        %Generate evasive trajectory
                         Jmin = J;
                         kevade_opt = kevade;
                         gamevade_opt = gamevade;
                         tturn_opt = tevade;
                         Xa_evade_opt = Xa_evade;
                         prim_opt = ii;
                         flyMode = 'evade';
                     end
                 end
         end
         if ~strcmp(flyMode,'evade')
             tTrajec=t:deltatTrajec:tendsim;
             newtrajec=genTrajecLinear(Xa,tTrajec',leveltrim);
             disp('Avoidance failed, could not find valid trajectory')
         else
             newtrajec_prim = primArray{prim_opt}';
             newtrajec_prim(2:4,:)=newtrajec_prim(2:4,:)+Xa(1:3);
             tfin=t+tturn_opt;
             tTrajec=tfin:deltatTrajec:tendsim;
             if isempty(tTrajec)
                 newtrajec = newtrajec_prim;
             else
                newtrajec_lin = genTrajecLinear(Xa_evade_opt,tTrajec',leveltrim);
                newtrajec = [newtrajec_prim newtrajec_lin(:,2:end)];
             end
             fprintf('Primitive %d selected (k = %f, gamma = %f) at t = %f to be performed for %f s\n',prim_opt,kevade_opt,gamevade_opt,t,tturn_opt)
         end
     elseif tca<0 && strcmp(flyMode,'evade') && t>=tfin
         fprintf('Replanning maneuver attempted at t = %f\n',t)
         tTrajec=t:deltatTrajec:tendsim;
         newtrajec=genTrajecLinear(Xa,tTrajec',leveltrim);
         psi_goal = atan2(Xgoal(2)-Xa(2),Xgoal(1)-Xa(1));
         flyMode = 'linear';
     else
         newtrajec = trajec(:,2:end); %do nothing, keep current plan
     end
 else
     flag_wc = false;
     newtrajec = trajec(:,2:end); %do nothing, keep current plan
 end
 
end