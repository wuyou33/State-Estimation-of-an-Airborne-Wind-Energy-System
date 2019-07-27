% Inverse calculation of "calculatePositionFromLineAngle.m"
% Calculates the line angles, when the position of the kite is known
% This function is used to simulate the tether sag
% sources: 
%   -[1]:  center2015 - Control Line Equations
% Arguments:
%               - l: tether length
%               -x,y,z: position of the kite
function [theta_base, theta_kite, phi, C1, b] = calculateAnglesFromPosition(xinput, yinput, zinput, l)
    xpos = sqrt(xinput.^2 + yinput.^2); % Ground distance
    ypos = zinput;
  
    fun = @(x)functionHandle(x, xpos, ypos, l);
    x0 = [0.0612, 0.612];
    % don't know how to vectorice
    %x0temp = ones(length(xpos), 2);
    %x0temp = x0temp .* [0.0612, 0.612];
    %x = arrayfun(@(x0) fsolve(fun, x0), x0temp);
    x = fsolve(fun, x0);
    
    C1 = x(1);
    b = x(2);
    theta_base = atan(sinh(C1));
    theta_kite = atan(sinh(b));
    phi = atan(yinput / xinput);
    
end

function F = functionHandle(x, x1, y1, l)
    C1 = x(1);
    b = x(2);
    F(1) = -y1+(cosh(b)*x1)/(b-C1)-(cosh(C1)*x1)/(b-C1);
    F(2) = l-(exp(-b)*(exp(2*b+C1)*x1+(1-exp(2*C1))*exp(b)*x1-exp(C1)*x1))/(2*exp(C1)*b-2*C1*exp(C1));
end