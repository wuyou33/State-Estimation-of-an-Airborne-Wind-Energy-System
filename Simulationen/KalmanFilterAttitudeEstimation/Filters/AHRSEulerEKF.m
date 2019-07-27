classdef AHRSEulerEKF < handle
    %% EKF AHRS euler angles [1]
    % AHRS in with euler angles
    % Input u: angular velocity in body coords
    % Measurement y: acceleration in body coords, magnetometer measurements
    % y: [acc_x, acc_y, acc_z, m_x, m_y, m_z]';
    % Version: 1.0
    properties
        EKF;
        sampleTime;
        R;
        Q;
        delta_i;
        delta_d;
        g_mps2;
    end
    
    methods
        
        function obj = AHRSEulerEKF(sampleTime, delta_i, delta_d, g_mps2, x0, P0, sigma_2_w, sigma_2_a, sigma_2_b_w, sigma_2_yaw) 
            obj.Q  = diag([sigma_2_w',sigma_2_b_w']); % contious-time process noise covariance matrix
            obj.R = diag([sigma_2_a',sigma_2_yaw']); % continous-time measurement noise covariance matrix
            obj.sampleTime = sampleTime;
            obj.EKF = EKF(sampleTime, obj.Q, obj.R, P0, x0);
            obj.delta_i = delta_i;
            obj.delta_d = delta_d;
            obj.g_mps2 = g_mps2;
        end
        
        % estimate states
        function performEKF(obj, u, y)
            y(1:3) = y(1:3) * obj.g_mps2;
            [x_dot, G_lin] = obj.x_dot_and_G_lin(obj.EKF.x, u);
            obj.EKF.predictorStep(x_dot, G_lin);
            [y_hat, H_lin] = obj.y_hat_and_H_lin(obj.EKF.x, u);
            obj.EKF.correctorStep(y, y_hat, H_lin);
            
            obj.EKF.x(1:3) = wrapToPi(obj.EKF.x(1:3));
        end
        
        % return states
        function [gamma_i, b_w] = states(obj)
            gamma_i = obj.EKF.x(1:3);
            b_w = obj.EKF.x(4:6);
        end

        % x_dot and G were calculated in maxima:
        % "KalmanfilterLagebestimmungEKF_BARTONSystemmatrizen/Systemmatrix.wxmx"
        
        % calculate the derivative of the states
        function [x_dot, G_lin] = x_dot_and_G_lin(obj, x,w)
            
            w_phi = w(1);
            w_theta = w(2);
            w_psi = w(3);

            phi = x(1);
            theta = x(2);
            psi = x(3);
            b_w_phi = x(4);
            b_w_theta = x(5);
            b_w_psi = x(6);
            x_dot = [(sin(phi)*sin(theta)*w_theta+cos(phi)*sin(theta)*w_psi+cos(theta)*w_phi+(-b_w_theta*sin(phi)-b_w_psi*cos(phi))*sin(theta)-b_w_phi*cos(theta))/cos(theta);
                    cos(phi)*w_theta-sin(phi)*w_psi+b_w_psi*sin(phi)-b_w_theta*cos(phi);
                    (sin(phi)*w_theta+cos(phi)*w_psi-b_w_theta*sin(phi)-b_w_psi*cos(phi))/cos(theta);
                    0;
                    0;
                    0];
                
            G_lin = [(cos(phi)*sin(theta)*w_theta-sin(phi)*sin(theta)*w_psi+(b_w_psi*sin(phi)-b_w_theta*cos(phi))*sin(theta))/cos(theta),	(sin(theta)*(sin(phi)*sin(theta)*w_theta+cos(phi)*sin(theta)*w_psi+cos(theta)*w_phi+(-b_w_theta*sin(phi)-b_w_psi*cos(phi))*sin(theta)-b_w_phi*cos(theta)))/cos(theta)^2+(sin(phi)*cos(theta)*w_theta+cos(phi)*cos(theta)*w_psi-sin(theta)*w_phi+b_w_phi*sin(theta)+(-b_w_theta*sin(phi)-b_w_psi*cos(phi))*cos(theta))/cos(theta),	0,	-1,	-(sin(phi)*sin(theta))/cos(theta),	-(cos(phi)*sin(theta))/cos(theta);
                    -sin(phi)*w_theta-cos(phi)*w_psi+b_w_theta*sin(phi)+b_w_psi*cos(phi),	0,	0,	0,	-cos(phi),	sin(phi);
                    (cos(phi)*w_theta-sin(phi)*w_psi+b_w_psi*sin(phi)-b_w_theta*cos(phi))/cos(theta),	(sin(theta)*(sin(phi)*w_theta+cos(phi)*w_psi-b_w_theta*sin(phi)-b_w_psi*cos(phi)))/cos(theta)^2,	0,	0,	-sin(phi)/cos(theta),	-cos(phi)/cos(theta);
                    0,	0,	0,	0,	0,	0;
                    0,	0,	0,	0,	0,	0;
                    0,	0,	0,	0,	0,	0];
        end
        
        function [y_hat, H_lin] = y_hat_and_H_lin(obj, x, u)
            
            delta_i = obj.delta_i;
            delta_d = obj.delta_d;
            g_mps2 = obj.g_mps2;
            phi = x(1);
            theta = x(2);
            psi = x(3);
            b_w_phi = x(4); % bias
            b_w_theta = x(5);
            b_w_psi = x(6);

            w_phi= u(1);
            w_theta = u(2);
            w_psi = u(3);

            % in earth frame
            vx_ins = u(4);
            vy_ins = u(5);
            vz_ins = u(6);

            % from "Fundamentals of Small Unmanned Aircraft Flight"
            v_b_ref = sqrt(vx_ins^2+vy_ins^2);

            % velocities in body coords
            % a multikopter can also flight sidewart
            v_x = v_b_ref;
            v_y = 0;
            v_z = 0;

            % yhat and H_lin calculated in
            % "KalmanfilterLagebestimmungEKF_BARTONSystemmatrizen/Systemmatrix.wxmx"

            y_hat =	[v_z*(w_theta-b_w_theta)+v_y*(b_w_psi-w_psi)-g_mps2*sin(theta);
                    v_x*(w_psi-b_w_psi)+v_z*(w_phi-b_w_phi)+g_mps2*sin(phi)*cos(theta);
                    v_x*(b_w_theta-w_theta)+v_y*(w_phi-b_w_phi)+g_mps2*cos(phi)*cos(theta);
                    sin((pi*delta_i)/180)*sin(theta)-sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*sin(psi)*cos(theta)+cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*cos(psi)*cos(theta);
                    -sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(sin(phi)*sin(psi)*sin(theta)+cos(phi)*cos(psi))+cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(sin(phi)*cos(psi)*sin(theta)-cos(phi)*sin(psi))-sin((pi*delta_i)/180)*sin(phi)*cos(theta);
                    -sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(cos(phi)*sin(psi)*sin(theta)-sin(phi)*cos(psi))+cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(cos(phi)*cos(psi)*sin(theta)+sin(phi)*sin(psi))-sin((pi*delta_i)/180)*cos(phi)*cos(theta)];

            H_lin = [0,	-g_mps2*cos(theta),	0,	0,	-v_z,	v_y;
                    g_mps2*cos(phi)*cos(theta),	-g_mps2*sin(phi)*sin(theta),	0,	-v_z,	0,	-v_x;
                    -g_mps2*sin(phi)*cos(theta),	-g_mps2*cos(phi)*sin(theta),	0,	-v_y,	v_x,	0;
                    0,	sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*sin(psi)*sin(theta)-cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*cos(psi)*sin(theta)+sin((pi*delta_i)/180)*cos(theta),	-cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*sin(psi)*cos(theta)-sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*cos(psi)*cos(theta),	0,	0,	0;
                    -sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(cos(phi)*sin(psi)*sin(theta)-sin(phi)*cos(psi))+cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(cos(phi)*cos(psi)*sin(theta)+sin(phi)*sin(psi))-sin((pi*delta_i)/180)*cos(phi)*cos(theta),	sin((pi*delta_i)/180)*sin(phi)*sin(theta)-sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*sin(phi)*sin(psi)*cos(theta)+cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*sin(phi)*cos(psi)*cos(theta),	cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(-sin(phi)*sin(psi)*sin(theta)-cos(phi)*cos(psi))-sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(sin(phi)*cos(psi)*sin(theta)-cos(phi)*sin(psi)),	0,	0,	0;
                    -sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(-sin(phi)*sin(psi)*sin(theta)-cos(phi)*cos(psi))+cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(cos(phi)*sin(psi)-sin(phi)*cos(psi)*sin(theta))+sin((pi*delta_i)/180)*sin(phi)*cos(theta),	sin((pi*delta_i)/180)*cos(phi)*sin(theta)-sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*cos(phi)*sin(psi)*cos(theta)+cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*cos(phi)*cos(psi)*cos(theta),	cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(sin(phi)*cos(psi)-cos(phi)*sin(psi)*sin(theta))-sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(cos(phi)*cos(psi)*sin(theta)+sin(phi)*sin(psi)),	0,	0,	0];

        end
    end

end