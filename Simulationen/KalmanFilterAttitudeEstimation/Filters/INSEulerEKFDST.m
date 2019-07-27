%% INS AHRS euler angles [1]
% INS with euler angles, but respects the different sample times. The
% INSEulerEKF doesn't this
% Version: 1.0
% For the derivation of the matrices see "Systemmatrix.wxm"
classdef INSEulerEKFDST < handle

    properties
        EKF;
        sampleTime;
        delta_i;
        delta_d;
        g_mps2;
    end
    
    methods
        
        function obj = INSEulerEKFDST(sampleTime, delta_i, delta_d, g_mps2, x0, P0, sigma_2_process, sigma_2_measurement) 
            obj.sampleTime = sampleTime;
            obj.delta_i = delta_i;
            obj.delta_d = delta_d;
            obj.g_mps2 = g_mps2;
            obj.EKF = EKF(sampleTime, sigma_2_process, sigma_2_measurement, P0, x0);
        end
        
        % estimate states
        function performEKF(obj, u, y, available)
            if isnan(obj.EKF.x(1))
                display('isnan');
            end
            u(4:6) = u(4:6) * obj.g_mps2;
            [x_dot, G_lin] = obj.x_dot_and_G_lin(obj.EKF.x, u);
            obj.EKF.predictorStep(x_dot, G_lin);
            
            %for i=1:4 % gps, magnetometer, baro, line
                [y_meas, y_hat, H_lin] = obj.y_hat_and_H_lin(y, obj.EKF.x, u, i, available);
                obj.EKF.correctorStep(y_meas, y_hat, H_lin);
            %end
            
            obj.EKF.x(1:3) = wrapToPi(obj.EKF.x(1:3));
            if isnan(obj.EKF.x(1))
                display('isnan');
            end
        end
        
        % return states
        function [gamma_i, pos, vel, b_w] = states(obj)
            gamma_i = obj.EKF.x(1:3);
            pos = obj.EKF.x(4:6);
            vel = obj.EKF.x(7:9);
            b_w = obj.EKF.x(10:12);
        end

        % x_dot and G were calculated in maxima:
        % "KalmanfilterLagebestimmungEKF_BARTONSystemmatrizen/Systemmatrix.wxmx"
        
        % calculate the derivative of the states
        function [x_dot, G_lin] = x_dot_and_G_lin(obj, x, u)
            phi = x(1);
            theta = x(2);
            psi = x(3);
            px = x(4);
            py = x(5);
            pz = x(6);
            vx = x(7);
            vy = x(8);
            vz = x(9);
            b_w_phi = x(10);
            b_w_theta = x(11);
            b_w_psi = x(12);
            
            w_phi = u (1);
            w_theta = u (2);
            w_psi = u (3);
            a_bx = u (4);
            a_by = u (5);
            a_bz = u (6);
            
            g_mps2 = obj.g_mps2;
            
            x_dot = [(sin(phi)*sin(theta)*(w_theta-b_w_theta))/cos(theta)+(cos(phi)*sin(theta)*(w_psi-b_w_psi))/cos(theta)+w_phi-b_w_phi;
                    cos(phi)*(w_theta-b_w_theta)-sin(phi)*(w_psi-b_w_psi);
                    (sin(phi)*(w_theta-b_w_theta))/cos(theta)+(cos(phi)*(w_psi-b_w_psi))/cos(theta);
                    vx;
                    vy;
                    vz;
                    a_by*(sin(phi)*cos(psi)*sin(theta)-cos(phi)*sin(psi))+a_bz*(cos(phi)*cos(psi)*sin(theta)+sin(phi)*sin(psi))+a_bx*cos(psi)*cos(theta);
                    a_by*(sin(phi)*sin(psi)*sin(theta)+cos(phi)*cos(psi))+a_bz*(cos(phi)*sin(psi)*sin(theta)-sin(phi)*cos(psi))+a_bx*sin(psi)*cos(theta);
                    -a_bx*sin(theta)+a_by*sin(phi)*cos(theta)+a_bz*cos(phi)*cos(theta)-g_mps2;
                    0;
                    0;
                    0];
                
            G_lin = [(cos(phi)*sin(theta)*(w_theta-b_w_theta))/cos(theta)-(sin(phi)*sin(theta)*(w_psi-b_w_psi))/cos(theta), (sin(phi)*sin(theta)^2*(w_theta-b_w_theta))/cos(theta)^2+sin(phi)*(w_theta-b_w_theta)+(cos(phi)*sin(theta)^2*(w_psi-b_w_psi))/cos(theta)^2+cos(phi)*(w_psi-b_w_psi), 0, 0, 0, 0, 0, 0, 0, -1, -(sin(phi)*sin(theta))/cos(theta), -(cos(phi)*sin(theta))/cos(theta);
                    -sin(phi)*(w_theta-b_w_theta)-cos(phi)*(w_psi-b_w_psi), 0, 0, 0, 0, 0, 0, 0, 0, 0, -cos(phi), sin(phi);
                    (cos(phi)*(w_theta-b_w_theta))/cos(theta)-(sin(phi)*(w_psi-b_w_psi))/cos(theta), (sin(phi)*sin(theta)*(w_theta-b_w_theta))/cos(theta)^2+(cos(phi)*sin(theta)*(w_psi-b_w_psi))/cos(theta)^2, 0, 0, 0, 0, 0, 0, 0, 0, -sin(phi)/cos(theta), -cos(phi)/cos(theta);
                    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
                    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;
                    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0;
                    a_bz*(cos(phi)*sin(psi)-sin(phi)*cos(psi)*sin(theta))+a_by*(cos(phi)*cos(psi)*sin(theta)+sin(phi)*sin(psi)), -a_bx*cos(psi)*sin(theta)+a_by*sin(phi)*cos(psi)*cos(theta)+a_bz*cos(phi)*cos(psi)*cos(theta), a_by*(-sin(phi)*sin(psi)*sin(theta)-cos(phi)*cos(psi))+a_bz*(sin(phi)*cos(psi)-cos(phi)*sin(psi)*sin(theta))-a_bx*sin(psi)*cos(theta), 0, 0, 0, 0, 0, 0, 0, 0, 0;
                    a_bz*(-sin(phi)*sin(psi)*sin(theta)-cos(phi)*cos(psi))+a_by*(cos(phi)*sin(psi)*sin(theta)-sin(phi)*cos(psi)), -a_bx*sin(psi)*sin(theta)+a_by*sin(phi)*sin(psi)*cos(theta)+a_bz*cos(phi)*sin(psi)*cos(theta), a_by*(sin(phi)*cos(psi)*sin(theta)-cos(phi)*sin(psi))+a_bz*(cos(phi)*cos(psi)*sin(theta)+sin(phi)*sin(psi))+a_bx*cos(psi)*cos(theta), 0, 0, 0, 0, 0, 0, 0, 0, 0;
                    a_by*cos(phi)*cos(theta)-a_bz*sin(phi)*cos(theta), -a_by*sin(phi)*sin(theta)-a_bz*cos(phi)*sin(theta)-a_bx*cos(theta), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

        end
        
        function [y_meas, y_hat, H_lin] = y_hat_and_H_lin(obj,y_meas, x, u, nr, available)
            phi = x(1);
            theta = x(2);
            psi = x(3);
            px = x(4);
            py = x(5);
            pz = x(6);
            vx = x(7);
            vy = x(8);
            vz = x(9);
            b_w_phi = x(10);
            b_w_theta = x(11);
            b_w_psi = x(12);
            
            delta_i = obj.delta_i;
            delta_d = obj.delta_d;
            
            y_hat = [sin((pi*delta_i)/180)*sin(theta)-sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*sin(psi)*cos(theta)+cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*cos(psi)*cos(theta);
                    -sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(sin(phi)*sin(psi)*sin(theta)+cos(phi)*cos(psi))+cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(sin(phi)*cos(psi)*sin(theta)-cos(phi)*sin(psi))-sin((pi*delta_i)/180)*sin(phi)*cos(theta);
                    -sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(cos(phi)*sin(psi)*sin(theta)-sin(phi)*cos(psi))+cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(cos(phi)*cos(psi)*sin(theta)+sin(phi)*sin(psi))-sin((pi*delta_i)/180)*cos(phi)*cos(theta);
                    pz;
                    px;
                    py;
                    pz;
                    px;
                    py;
                    pz;
                    vx;
                    vy;
                    vz];
                
            y_hat(~available) = 0;
            y_meas(~available) = 0; 
                    
           H_lin = 	[0, sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*sin(psi)*sin(theta)-cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*cos(psi)*sin(theta)+sin((pi*delta_i)/180)*cos(theta), -cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*sin(psi)*cos(theta)-sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*cos(psi)*cos(theta), 0, 0, 0, 0, 0, 0, 0, 0, 0;
                    -sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(cos(phi)*sin(psi)*sin(theta)-sin(phi)*cos(psi))+cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(cos(phi)*cos(psi)*sin(theta)+sin(phi)*sin(psi))-sin((pi*delta_i)/180)*cos(phi)*cos(theta), sin((pi*delta_i)/180)*sin(phi)*sin(theta)-sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*sin(phi)*sin(psi)*cos(theta)+cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*sin(phi)*cos(psi)*cos(theta), cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(-sin(phi)*sin(psi)*sin(theta)-cos(phi)*cos(psi))-sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(sin(phi)*cos(psi)*sin(theta)-cos(phi)*sin(psi)), 0, 0, 0, 0, 0, 0, 0, 0, 0;
                    -sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(-sin(phi)*sin(psi)*sin(theta)-cos(phi)*cos(psi))+cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(cos(phi)*sin(psi)-sin(phi)*cos(psi)*sin(theta))+sin((pi*delta_i)/180)*sin(phi)*cos(theta), sin((pi*delta_i)/180)*cos(phi)*sin(theta)-sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*cos(phi)*sin(psi)*cos(theta)+cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*cos(phi)*cos(psi)*cos(theta), cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(sin(phi)*cos(psi)-cos(phi)*sin(psi)*sin(theta))-sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(cos(phi)*cos(psi)*sin(theta)+sin(phi)*sin(psi)), 0, 0, 0, 0, 0, 0, 0, 0, 0;
                    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
                    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
                    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
                    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
                    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
                    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
                    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
                    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
                    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;
                    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0];
            
       end
    end

end