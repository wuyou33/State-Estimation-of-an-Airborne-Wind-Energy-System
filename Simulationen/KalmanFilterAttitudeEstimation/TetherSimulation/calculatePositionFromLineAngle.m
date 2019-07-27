% Calculates the Position of the kite, when the line Angles and the length
% of the tether is present
% sources: 
%   -[1]:  center2015 - Control Line Equations
% TODO: vectorize this function
function [x, y, z, a, C1, C2] = calculatePositionFromLineAngle(theta_base, theta_kite, phi, l)
 
C1 = asinh(tan(theta_base));

XK = (2*l*exp(asinh(sin(theta_kite)/cos(theta_kite))+C1)*asinh(sin(theta_kite)/cos(theta_kite))-2*C1*l*exp(asinh(sin(theta_kite)/cos(theta_kite))+C1))/(exp(2*asinh(sin(theta_kite)/cos(theta_kite))+C1)-exp(asinh(sin(theta_kite)/cos(theta_kite))+2*C1)+exp(asinh(sin(theta_kite)/cos(theta_kite)))-exp(C1));

a = XK/(asinh(tan(theta_kite))-C1);

C2 = -a*cosh(C1);

z = C2 + a * cosh(XK/a + C1);

x = XK * cos(phi);
y = XK * sin(phi);

end