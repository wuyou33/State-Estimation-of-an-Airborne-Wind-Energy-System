%% 
% Plot the catenary curve, when two line angle sensors and the length of
% the tether is available
% for the reverse step see LineAnglePositionEstimationNASA.wxm
% sources: 
%   -[1]:  center2015 - Control Line Equations

% For the derivation see LineAnglePositionEstimationNASA.wxm
%%

clear; clc; close all;

%% Settings begin
 % set to one, if the functions 
 % calculatePositioinFromLineAngle and calculateAnglesFromPostion should be
 % tested
testCalculationFunctions = 0;
plotLineCalculatedFromAngles = 1;
plotLineCalculatedFromPosition = 1;

if plotLineCalculatedFromAngles
    l = 180; % [m] total length of the tether

    theta_b = 10 * pi / 180; % [rad] angle measured by the lineAngleSensor on the base
    theta_k = 50 * pi / 180; % [rad] angle measured by the line angle sensor on the kite
    phi = 0;
    %% Settings end

    % this function calculates only the position of the last point
    [x_test, y_test, z_test, a, C1, C2] = calculatePositionFromLineAngle(theta_b, theta_k, phi, l); % only 2D here

    C1 = asinh(tan(theta_b));

    XK = sqrt(x_test^2 + y_test^2); %(2*l*exp(asinh(sin(theta_k)/cos(theta_k))+C1)*asinh(sin(theta_k)/cos(theta_k))-2*C1*l*exp(asinh(sin(theta_k)/cos(theta_k))+C1))/(exp(2*asinh(sin(theta_k)/cos(theta_k))+C1)-exp(asinh(sin(theta_k)/cos(theta_k))+2*C1)+exp(asinh(sin(theta_k)/cos(theta_k)))-exp(C1));
    x = linspace(0, XK, 1000);


    C2 = -a*cosh(C1);
    y = C2 + a * cosh(x/a + C1);

    figure('Name', 'LineCalculatedFromAngles');
    plot(x, y);
end

if plotLineCalculatedFromPosition
    l = 154.5461;
    xpos = 92.191062927246100;
    ypos = 1.158024673461914e+02;
    zpos = 24.576808929443360;
    
    XK = sqrt(xpos^2+ypos^2);
    x = linspace(0, XK, 1000);
    
    [theta_base, theta_kite, phi_sag, C1, b] = calculateAnglesFromPosition(xpos, ypos, zpos, l);
    
    a = XK/(b-C1);
    C2 = -a*cosh(C1);
    
    y = C2 + a*cosh(x/a + C1);
    
    figure('Name', 'LineCalculatedFromPosition');
    plot(x, y);
    
    clear xpos ypos zpos l C2 y x a theta_base theta_kite phi_sag C1 b
end

if testCalculationFunctions
    theta_b_test = 17*pi/180; % [rad]
    theta_k_test = 80*pi/180; % [rad]
    phi_test = 10 * pi/180; % [rad]
    l_test = 150; % tether lenght [m]
    
    [xpos_test, ypos_test, zpos_test] = calculatePositionFromLineAngle(theta_b_test, theta_k_test, phi_test, l_test);
    
    [theta_b_sol, theta_k_sol, phi_sol] = calculateAnglesFromPosition(xpos_test, ypos_test, zpos_test, l_test);
    
    if abs(theta_b_sol - theta_b_test) < 0.01 && abs(theta_k_sol - theta_k_test) < 0.01 && abs(phi_sol - phi_test) < 0.01
        display('Test successfully completed');
    else
        warning('Test not successfully completed');
        theta_b_sol
        theta_b_test
        theta_k_sol
        theta_k_test
        phi_sol
        theta_k
    end
    clear theta_b_test theta_k_test phi_test l_test xpos_test ypos_test zpos_teset theta_b_sol theta_k_sol phi_sol;
end

