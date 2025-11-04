function trajec=genTrajecLinear(Xa,tTrajec,leveltrim)
global Va
%trajec stores [t;x;y;z;Va;psi;k1;gamma1] (same as ROS)
%Xa is [x y z psi gamma]
trajec = zeros(7,length(tTrajec));

%trim is p q r u v w phi theta psi h dE dA dR dT
%trajec is t x y z Va psi k gamma p q r phi theta dE dA dR dT

%Assume FPA = 0
trajec(1,:) = tTrajec;
trajec(2,:) = Xa(1)+Va*cos(Xa(4)).*(tTrajec-tTrajec(1));
trajec(3,:) = Xa(2)+Va*sin(Xa(4)).*(tTrajec-tTrajec(1));
trajec(4,:) = Xa(3)+Va*0.*(tTrajec-tTrajec(1));
trajec(5,:) = Va;
trajec(6,:) = Xa(4);
trajec(7,:) = 0;
trajec(8,:) = 0;
trajec(9,:) = leveltrim(1);
trajec(10,:) = leveltrim(2);
trajec(11,:) = leveltrim(3); 
trajec(12,:) = leveltrim(7);
trajec(13,:) = leveltrim(8);
trajec(14,:) = leveltrim(11);
trajec(15,:) = leveltrim(12);
trajec(16,:) = leveltrim(13);
trajec(17,:) = leveltrim(14);
end