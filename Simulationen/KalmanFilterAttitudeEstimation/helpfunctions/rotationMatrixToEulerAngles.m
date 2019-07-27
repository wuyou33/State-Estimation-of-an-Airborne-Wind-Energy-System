% see EulerAngleCalculation.wxm or 
% Diebel2006 eq. 72
function gamma = rotationMatrixToEulerAngles(M)

    phi= atan2(M(2,3), M(3,3));
    theta= -asin(M(1,3));
    psi= atan2(M(1,2), M(1,1));

    gamma = [ phi;
              theta;
              psi ];
end