function J = getJMotionPlanner(X0,X1,Xg,F,wvec)

%%Inputs
%X0 - current configuration of aircraft [x y z psi gamma]
%X1 - proposed configuration of aircraft after maneuver [x y z psi gamma]
%Xg - goal configuration of aircraft [x y z psi gamma]
%F - matrix of geofence where each row is a geofence segment [px py mx my]
%wvec - vector of weights for each cost term

%%Outputs
%J - cost of maneuver

%J1, distance traveled during maneuver
J1 = norm(X1(1:3)-X0(1:3));

%J2, distance from goal
J2 = norm(Xg(1:3)-X1(1:3));

%J3, distance from geofence
nseg = size(F,1);
q = X1(1:2);
J3 = 0;
for ii = 1:nseg
    p = F(ii,1:2)';
    m = F(ii,3:4)';
    geoDist = dot(q-p,m);
    J3 = J3 + log(geoDist); %Option 2: using the geometric mean
end
J3=exp(-1/nseg*J3);%convert back to distance (geometric mean)
J = wvec(1)*J1+wvec(2)*J2+wvec(3)*J3;
end