classdef AHRSQuaternionEKF < matlab.System %handle
    %% EKF AHRS Quaternion [1]
    % AHRS in with quaternions
    % Input u: angular velocity in body coords
    % Measurement y: acceleration in body coords, magnetometer measurements
    % y: [acc_x, acc_y, acc_z, m_x, m_y, m_z]';
    % Version: 1.0
    
    % These variables must be initialized in the simulink model
    properties
       sigma_2_quaternion; 
       sigma_2_b_w;
       sigma_2_a;
       sigma_2_yaw;
       TA;
       delta_i;
       delta_d;
       g_mps2;
       P0;
       x0;
    end
    
    % Public, non-tunable properties
    properties(Nontunable)

    end

    properties(DiscreteState)

    end

    properties(Access = private) % These variables must be initialised. Here or in the setupImpl function
        EKF;
        sampleTime;
        R;
        Q;
        y_hat;
    end
    
    methods (Access = protected)
        
        %% Common functions
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
            obj.Q  = diag([obj.sigma_2_quaternion',obj.sigma_2_b_w']); % contious-time process noise covariance matrix
            obj.R = diag([obj.sigma_2_a',obj.sigma_2_yaw']); % continous-time measurement noise covariance matrix
            obj.sampleTime = obj.TA;
            obj.EKF = EKF(obj.TA, obj.Q, obj.R, obj.P0, obj.x0);
        end
        
        function states = stepImpl(obj,u, y)
            % Implement algorithm.
            obj.performEKF(u, y);
            [q_hat, b_w] = obj.states();
            q = Quaternion(q_hat);
            gamma_i_ahrsQuaternion = rad2deg(q.quatToEuler());
            b_w_ahrsQuaternion = b_w;
            states = [gamma_i_ahrsQuaternion;
                        b_w_ahrsQuaternion];
        end
    end
    
    methods
        function obj = AHRSQuaternionEKF(varargin) %sampleTime, x0, P0, delta_i, delta_d, g_mps2, sigma_2_quaternion, sigma_2_a, sigma_2_b_w, sigma_2_yaw)               
            % Support name-value pair arguments when constructing object
        end
        
        function initFilter(obj, varargin)
            sampleTime = varargin{1};
            x0 = varargin{2};
            P0 = varargin{3};
            delta_i = varargin{4};
            delta_d = varargin{5};
            g_mps2 = varargin{6};
            sigma_2_quaternion = varargin{7};
            sigma_2_a = varargin{8};
            sigma_2_b_w = varargin{9};
            sigma_2_yaw = varargin{10};
            
            obj.Q  = diag([sigma_2_quaternion',sigma_2_b_w']); % contious-time process noise covariance matrix
            obj.R = diag([sigma_2_a',sigma_2_yaw']); % continous-time measurement noise covariance matrix
            obj.sampleTime = sampleTime;
            obj.EKF = EKF(sampleTime, obj.Q, obj.R, P0, x0);
            obj.delta_i = delta_i;
            obj.delta_d = delta_d;
            obj.g_mps2 = g_mps2;
        end

        % estimate states
        function performEKF(obj, u, y)
            y(1:3) = y(1:3) * obj.g_mps2; % [g] --> [m/s^2]
            [x_dot, G_lin] = obj.x_dot_and_G_lin(obj.EKF.x, u);
            obj.EKF.predictorStep(x_dot, G_lin);
            [y_hat, H_lin] = obj.y_hat_and_H_lin(obj.EKF.x, u); % oder doch nur x?
            obj.y_hat = y_hat; % just for debugging
            obj.EKF.correctorStep(y, y_hat, H_lin);
            
            % due to the calculations the norm will be removed
            % the quaternions must have unit norm!
            obj.EKF.x(1:4) = obj.EKF.x(1:4)/norm(obj.EKF.x(1:4)); 
        end
        
        % return states
        function [q, b_w] = states(obj)
            q = obj.EKF.x(1:4);
            b_w = obj.EKF.x(5:7);
        end
        
        function y_hat = estimatedOutput(obj)
            y_hat = obj.y_hat;
        end

        % x_dot and G were calculated in maxima:
        % "KalmanfilterLagebestimmungEKF_BARTONSystemmatrizen/Systemmatrix.wxmx"
        
        % calculate the derivative of the states
        function [x_dot, G_lin] = x_dot_and_G_lin(obj, x,w)
            
            w_phi_m = w(1);
            w_theta_m = w(2);
            w_psi_m = w(3);

            q_0 = x(1);
            q_1 = x(2);
            q_2 = x(3);
            q_3 = x(4);
            b_w_phi = x(5);
            b_w_theta = x(6);
            b_w_psi = x(7);

            % x_dot and A were calculated in maxima:
            % "KalmanfilterLagebestimmungEKF_BARTONSystemmatrizen/EKFQuaternion.wxmx"
            x_dot = [(-q_2*(w_theta_m-b_w_theta)-q_3*(w_psi_m-b_w_psi)-q_1*(w_phi_m-b_w_phi))/2;
                    (-q_3*(w_theta_m-b_w_theta)+q_2*(w_psi_m-b_w_psi)+q_0*(w_phi_m-b_w_phi))/2;
                    (q_0*(w_theta_m-b_w_theta)-q_1*(w_psi_m-b_w_psi)+q_3*(w_phi_m-b_w_phi))/2;
                    (q_1*(w_theta_m-b_w_theta)+q_0*(w_psi_m-b_w_psi)-q_2*(w_phi_m-b_w_phi))/2;
                    0;
                    0;
                    0];

            G_lin = [0,	(b_w_phi-w_phi_m)/2,	(b_w_theta-w_theta_m)/2,	(b_w_psi-w_psi_m)/2,	q_1/2,	q_2/2,	q_3/2;
                    (w_phi_m-b_w_phi)/2,	0,	(w_psi_m-b_w_psi)/2,	(b_w_theta-w_theta_m)/2,	-q_0/2,	q_3/2,	-q_2/2;
                    (w_theta_m-b_w_theta)/2,	(b_w_psi-w_psi_m)/2,	0,	(w_phi_m-b_w_phi)/2,	-q_3/2,	-q_0/2,	q_1/2;
                    (w_psi_m-b_w_psi)/2,	(w_theta_m-b_w_theta)/2,	(b_w_phi-w_phi_m)/2,	0,	q_2/2,	-q_1/2,	-q_0/2;
                    0,	0,	0,	0,	0,	0,	0;
                    0,	0,	0,	0,	0,	0,	0;
                    0,	0,	0,	0,	0,	0,	0];
        end
        
        function [y_hat, H_lin] = y_hat_and_H_lin(obj, x, u)
            
            delta_i = obj.delta_i;
            delta_d = obj.delta_d;
            g_mps2 = obj.g_mps2;

            % x: state vector 7x1
            q_0 = x(1);
            q_1 = x(2);
            q_2 = x(3);
            q_3 = x(4);
            b_w_phi = x(5); % bias
            b_w_theta = x(6);
            b_w_psi = x(7);

            w_phi_m= u(1);
            w_theta_m = u(2);
            w_psi_m = u(3);

            % in earth frame
            vx_ins = u(4);
            vy_ins = u(5);
            vz_ins = u(6);

            % from "Fundamentals of Small Unmanned Aircraft Flight"
            v_b_ref = sqrt(vx_ins^2+vy_ins^2);

            % velocities in body coords
            % aircraft flights only forward,
            % the velocity for a multikopter must be defined different
            v_x = v_b_ref;
            v_y = 0;
            v_z = 0;

            % yhat and H_lin calculated in
            % "KalmanfilterLagebestimmungEKF_BARTONSystemmatrizen/Systemmatrix.wxmx"

            y_hat = [v_z*(w_theta_m-b_w_theta)+v_y*(b_w_psi-w_psi_m)+2*g_mps2*(q_1*q_3-q_0*q_2);
                    v_x*(w_psi_m-b_w_psi)+v_z*(b_w_phi-w_phi_m)+2*g_mps2*(q_2*q_3+q_0*q_1);
                    v_x*(b_w_theta-w_theta_m)+v_y*(w_phi_m-b_w_phi)+g_mps2*(1-2*(q_2^2+q_1^2));
                    cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(1-2*(q_3^2+q_2^2))-2*sin((pi*delta_i)/180)*(q_1*q_3-q_0*q_2)-2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q_0*q_3+q_1*q_2);
                    -sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(1-2*(q_3^2+q_1^2))-2*sin((pi*delta_i)/180)*(q_2*q_3+q_0*q_1)+2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q_1*q_2-q_0*q_3);
                    -2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q_2*q_3-q_0*q_1)+2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q_1*q_3+q_0*q_2)-sin((pi*delta_i)/180)*(1-2*(q_2^2+q_1^2))];
            H_lin = [-2*g_mps2*q_2,	2*g_mps2*q_3,	-2*g_mps2*q_0,	2*g_mps2*q_1,	0,	-v_z,	v_y;
                    2*g_mps2*q_1,	2*g_mps2*q_0,	2*g_mps2*q_3,	2*g_mps2*q_2,	v_z,	0,	-v_x;
                    0,	-4*g_mps2*q_1,	-4*g_mps2*q_2,	0,	-v_y,	v_x,	0;
                    2*sin((pi*delta_i)/180)*q_2-2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*q_3,	-2*sin((pi*delta_i)/180)*q_3-2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*q_2,	-4*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*q_2-2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*q_1+2*sin((pi*delta_i)/180)*q_0,	-4*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*q_3-2*sin((pi*delta_i)/180)*q_1-2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*q_0,	0,	0,	0;
                    -2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*q_3-2*sin((pi*delta_i)/180)*q_1,	2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*q_2+4*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*q_1-2*sin((pi*delta_i)/180)*q_0,	2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*q_1-2*sin((pi*delta_i)/180)*q_3,	4*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*q_3-2*sin((pi*delta_i)/180)*q_2-2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*q_0,	0,	0,	0;
                    2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*q_2+2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*q_1,	2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*q_3+4*sin((pi*delta_i)/180)*q_1+2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*q_0,	-2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*q_3+4*sin((pi*delta_i)/180)*q_2+2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*q_0,	2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*q_1-2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*q_2,	0,	0,	0];
        end
    end

end