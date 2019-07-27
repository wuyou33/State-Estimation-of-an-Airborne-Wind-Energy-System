classdef AHRSEulerEKFHIComp < matlab.System
    %% EKF AHRS euler angles [1]
    % AHRS in with euler angles and magnet hard iron compensation
    % This filter includes the magnetic field in earth frame and the
    % magnetic field in body frame to compensate the hard iron distortions
    % Another benefit of this filter is, no inclination and declination is
    % used.
    % Derivation: "AHRS EKF (with hard iron calibration)"
    % States: 3 x euler angles, 3 x gyro bias, 3 x magnetic value in NED
    % frame, 3 x magnetic value in body frame
    % Input u: angular velocity in body coords
    % Measurement y: acceleration in body coords, magnetometer measurements
    % y: [acc_x, acc_y, acc_z, m_x, m_y, m_z]';
    % Version: 1.0
    
    % These variables must be initialized in the simulink model
    properties
        sampleTime;
        delta_i;
        delta_d;
        g_mps2;
        x0;
        P0;
        sigma_2_w;
        sigma_2_a;
        sigma_2_b_w;
        sigma_2_yaw;
    end
    
    properties (Access = protected)
        EKF;
        R;
        Q;
    end
    
    methods (Access = protected)
        
        function setupImpl(obj)
            obj.Q  = diag([obj.sigma_2_w',obj.sigma_2_b_w', obj.sigma_2_yaw', obj.sigma_2_yaw']); % contious-time process noise covariance matrix ; Assumption: the derivatives of the magentic values have the same variance than the magnetometer values it self (not correct)
            obj.R = diag([obj.sigma_2_a',obj.sigma_2_yaw']); % continous-time measurement noise covariance matrix
            obj.EKF = EKF(obj.sampleTime, obj.Q, obj.R, obj.P0, obj.x0);
        end
        
        function [attitude, b_w, m_e, m_b] = stepImpl(obj, u, y)
            % Implement algorithm
            obj.performEKF(u, y);
            states = obj.EKF.x;
            
            attitude = states(1:3);
            b_w = states(4:6);
            m_e = states(7:9);
            m_b = states(10:12);
        end
   end
        
   methods
        % Constructor must be empty for matlab.System. In Matlab call
        % initFilter after the object was created. In simulink setupImpl()
        % will be called
        function obj = AHRSEulerEKFHIComp()  
        end
        
        function initFilter(obj, TA, delta_i, delta_d, g_mps2, x0, P0, sigma_2_w_b, sigma_2_a_b, sigma_2_b_w, sigma_2_mag)
            obj.sampleTime = TA;
            obj.delta_i = delta_i;
            obj.delta_d = delta_d;
            obj.g_mps2 = g_mps2;
            obj.sigma_2_w = sigma_2_w_b;
            obj.sigma_2_a = sigma_2_a_b;
            obj.sigma_2_b_w = sigma_2_b_w;
            obj.sigma_2_yaw = sigma_2_mag;
            obj.x0 = x0;
            obj.P0 = P0;
            obj.setupImpl();
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
            x_dot = [(sin(phi)*sin(theta)*(w_theta-b_w_theta)+cos(phi)*sin(theta)*(w_psi-b_w_psi)+cos(theta)*(w_phi-b_w_phi))/cos(theta);
                    cos(phi)*(w_theta-b_w_theta)-sin(phi)*(w_psi-b_w_psi);
                    (sin(phi)*(w_theta-b_w_theta)+cos(phi)*(w_psi-b_w_psi))/cos(theta);
                    0;
                    0;
                    0;
                    0;
                    0;
                    0;
                    0;
                    0;
                    0];
                
            G_lin = [(cos(phi)*sin(theta)*(w_theta-b_w_theta)-sin(phi)*sin(theta)*(w_psi-b_w_psi))/cos(theta), (sin(theta)*(sin(phi)*sin(theta)*(w_theta-b_w_theta)+cos(phi)*sin(theta)*(w_psi-b_w_psi)+cos(theta)*(w_phi-b_w_phi)))/cos(theta)^2+(sin(phi)*cos(theta)*(w_theta-b_w_theta)+cos(phi)*cos(theta)*(w_psi-b_w_psi)-sin(theta)*(w_phi-b_w_phi))/cos(theta), 0, -1, -(sin(phi)*sin(theta))/cos(theta), -(cos(phi)*sin(theta))/cos(theta), 0, 0, 0, 0, 0, 0;
                    -sin(phi)*(w_theta-b_w_theta)-cos(phi)*(w_psi-b_w_psi), 0, 0, 0, -cos(phi), sin(phi), 0, 0, 0, 0, 0, 0;
                    (cos(phi)*(w_theta-b_w_theta)-sin(phi)*(w_psi-b_w_psi))/cos(theta), (sin(theta)*(sin(phi)*(w_theta-b_w_theta)+cos(phi)*(w_psi-b_w_psi)))/cos(theta)^2, 0, 0, -sin(phi)/cos(theta), -cos(phi)/cos(theta), 0, 0, 0, 0, 0, 0;
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
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
            mag_e_N = x(7);
            mag_e_E = x(8);
            mag_e_D = x(9);
            mag_b_x = x(10);
            mag_b_y = x(11);
            mag_b_z = x(12);

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
                    v_x*(w_psi-b_w_psi)+v_z*(b_w_phi-w_phi)+g_mps2*sin(phi)*cos(theta);
                    v_x*(b_w_theta-w_theta)+v_y*(w_phi-b_w_phi)+g_mps2*cos(phi)*cos(theta);
                    -mag_e_D*sin(theta)+mag_e_E*sin(psi)*cos(theta)+mag_e_N*cos(psi)*cos(theta)+mag_b_x;
                    mag_e_E*(sin(phi)*sin(psi)*sin(theta)+cos(phi)*cos(psi))+mag_e_N*(sin(phi)*cos(psi)*sin(theta)-cos(phi)*sin(psi))+mag_e_D*sin(phi)*cos(theta)+mag_b_y;
                    mag_e_E*(cos(phi)*sin(psi)*sin(theta)-sin(phi)*cos(psi))+mag_e_N*(cos(phi)*cos(psi)*sin(theta)+sin(phi)*sin(psi))+mag_e_D*cos(phi)*cos(theta)+mag_b_z];

            H_lin = [0, -g_mps2*cos(theta), 0, 0, -v_z, v_y, 0, 0, 0, 0, 0, 0;
                    g_mps2*cos(phi)*cos(theta), -g_mps2*sin(phi)*sin(theta), 0, v_z, 0, -v_x, 0, 0, 0, 0, 0, 0;
                    -g_mps2*sin(phi)*cos(theta), -g_mps2*cos(phi)*sin(theta), 0, -v_y, v_x, 0, 0, 0, 0, 0, 0, 0;
                    0, -mag_e_E*sin(psi)*sin(theta)-mag_e_N*cos(psi)*sin(theta)-mag_e_D*cos(theta), mag_e_E*cos(psi)*cos(theta)-mag_e_N*sin(psi)*cos(theta), 0, 0, 0, cos(psi)*cos(theta), sin(psi)*cos(theta), -sin(theta), 1, 0, 0;
                    mag_e_E*(cos(phi)*sin(psi)*sin(theta)-sin(phi)*cos(psi))+mag_e_N*(cos(phi)*cos(psi)*sin(theta)+sin(phi)*sin(psi))+mag_e_D*cos(phi)*cos(theta), -mag_e_D*sin(phi)*sin(theta)+mag_e_E*sin(phi)*sin(psi)*cos(theta)+mag_e_N*sin(phi)*cos(psi)*cos(theta), mag_e_N*(-sin(phi)*sin(psi)*sin(theta)-cos(phi)*cos(psi))+mag_e_E*(sin(phi)*cos(psi)*sin(theta)-cos(phi)*sin(psi)), 0, 0, 0, sin(phi)*cos(psi)*sin(theta)-cos(phi)*sin(psi), sin(phi)*sin(psi)*sin(theta)+cos(phi)*cos(psi), sin(phi)*cos(theta), 0, 1, 0;
                    mag_e_E*(-sin(phi)*sin(psi)*sin(theta)-cos(phi)*cos(psi))+mag_e_N*(cos(phi)*sin(psi)-sin(phi)*cos(psi)*sin(theta))-mag_e_D*sin(phi)*cos(theta), -mag_e_D*cos(phi)*sin(theta)+mag_e_E*cos(phi)*sin(psi)*cos(theta)+mag_e_N*cos(phi)*cos(psi)*cos(theta), mag_e_N*(sin(phi)*cos(psi)-cos(phi)*sin(psi)*sin(theta))+mag_e_E*(cos(phi)*cos(psi)*sin(theta)+sin(phi)*sin(psi)), 0, 0, 0, cos(phi)*cos(psi)*sin(theta)+sin(phi)*sin(psi), cos(phi)*sin(psi)*sin(theta)-sin(phi)*cos(psi), cos(phi)*cos(theta), 0, 0, 1];

        end
    end

end