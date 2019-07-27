% source:
% [1] Diebel2006 - Representing Attitude: Euler Angles, Unit Quaternions, and Rotation Vectors

function q = eulerAnglesToQuaterion(gamma)

phi = gamma(1,:);
theta = gamma(2,:);
psi = gamma(3,:);

q0 = cos(phi/2).*cos(theta/2).*cos(psi/2)+sin(phi/2).*sin(theta/2).*sin(psi/2);
q1 = -cos(phi/2).*sin(theta/2).*sin(psi/2)+cos(theta/2).*cos(psi/2).*sin(phi/2);
q2 = cos(phi/2).*cos(psi/2).*sin(theta/2)+sin(phi/2).*cos(theta/2).*sin(psi/2);
q3 = cos(phi/2).*cos(theta/2).*sin(psi/2)-sin(phi/2).*cos(psi/2).*sin(theta/2);

q = [q0;
     q1;
     q2;
     q3];
end