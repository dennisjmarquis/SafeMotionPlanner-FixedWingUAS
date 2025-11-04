function R_b2i= rotation_Body2Earth( anglesout)
phi=anglesout(1);
theta=anglesout(2);
psi=anglesout(3);
R_b2i=zeros(3,3);
R_b2i(1,1)= cos(theta)*cos(psi);
R_b2i(2,1)= cos(theta)*sin(psi);
R_b2i(3,1)= -sin(theta);

R_b2i(1,2)= sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi);
R_b2i(2,2)= sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi);
R_b2i(3,2)= sin(phi)*cos(theta);

R_b2i(1,3)= cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi);
R_b2i(2,3)= cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi);
R_b2i(3,3)= cos(phi)*cos(theta);

end
