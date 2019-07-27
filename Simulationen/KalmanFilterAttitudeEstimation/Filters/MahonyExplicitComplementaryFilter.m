% source: 
% [1]: Mahony et al. - Nonlinear Complementary Filters on the Special Orthogonal Group
% [2]: https://github.com/gaochq/IMU_Attitude_Estimator/blob/master/src/Mahony_Attitude.cpp
% [3]: Mattia Giurato - https://github.com/Murmele/Attitude_Estimation
% [4]: Maley2013 - Multiplicative Quaternion Extended Kalman Filtering for Nonspinning Guided Projectiles
% Version: 1.0

% TODO: 
% - Why the performance of psi is so bad, when delta_i and delta_d are not zero?
% - it seems, that the filter with centripetalcompensation is anymore
% everywhere stable. Have to check it. When I optimized again, everything
% was fine


classdef MahonyExplicitComplementaryFilter < handle %matlab.System
 
    % Public, tunable properties
    properties
        TA;
        KP = 10.049404236028266; % filter parameter
        KI =54.934135404052185; % filter parameter
        Kacc = 23.403584221427753; % [3]. Weight accelerometer
        Kmag = 43.458448798901200; % [3]. Weigth magnetometer
        % with these parameters the step response was anymore stable
%         KP = 17.359571688139190; % filter parameter
%         KI = 11.056899210324557; % filter parameter
%         Kacc = 0.847817744656889; % [3]. Weight accelerometer
%         Kmag = 1.669354841447321; % [3]. Weigth magnetometer
        q; % attitude
        b_w;
        delta_i;
        delta_d;
        g_mps2;
        gyro = [0;0;0];
    end
    
    % Public, non-tunable properties
    properties(Nontunable)

    end

    properties(DiscreteState)

    end
    
    % Pre-computed constants
    properties(Access = private)

    end
    
    %% Public methods
    methods (Access = public)
        % Constructor
        function obj = MahonyExplicitComplementaryFilter(TA, delta_i, delta_d, g_mps2, q0, b_w0)
            
            % Support name-value pair arguments when constructing object
            %setProperties(obj,nargin,varargin{:})
            obj.TA = TA;
            obj.q = q0;
            obj.b_w = b_w0;
            obj.delta_i = delta_i;
            obj.delta_d = delta_d;
            obj.g_mps2 = g_mps2;
        end
        
        function update(obj, accelerometer, magnetometer, gyro, position, velocity)
            % normalize
            v_acc = accelerometer/norm(accelerometer);
            v_mag = magnetometer/norm(magnetometer);

            % [1], eq. 29 f, hat{v}_i = \hat{R}^T
            m_i = Rmag_to_e([1;0;0], obj.delta_d, obj.delta_i);
            v_mag_hat = obj.Re_to_b_Quaternion(obj.q, m_i);
            a_i = [0; 0; obj.g_mps2];
            v_vel_b = [sqrt(velocity(1)^2+velocity(2)^2); 0;0];
            % must be multiplied by -1 I think, because the accelerometer
            % measures the centrifugal acceleration??
            v_acc_centripetal = obj.crossproductMatrix(gyro - obj.b_w)*v_vel_b; % ignored in the original paper.
            v_acc_hat = obj.Re_to_b_Quaternion(obj.q, a_i) + v_acc_centripetal; % last part is the centripetal acceleration
            v_acc_hat = v_acc_hat/obj.g_mps2; % [m/s^2] --> [g]
            v_acc_hat = v_acc_hat/norm(v_acc_hat);

            % [1], eq. 34
            acc_corr = obj.Kacc/2 * (v_acc*v_acc_hat' - v_acc_hat*v_acc');
            mag_corr = obj.Kmag/2 * (v_mag*v_mag_hat' - v_mag_hat*v_mag');
            ome_mes_x = (acc_corr + mag_corr); % cross product matrix
            % [1], eq. 48a
            ome_mes = - [ome_mes_x(3,2) ; % vector from ome_mes_x
                         ome_mes_x(1,3) ;
                         ome_mes_x(2,1)];

            % [1], eq. 48c
            obj.b_w = obj.b_w - obj.KI*ome_mes * obj.TA; % euler integrator

            % [1], eq 48b
            euler = 1;
            if euler
                q = Quaternion(obj.q);
                q_dot = 1/2 * q.mult([0; gyro - obj.b_w + obj.KP*ome_mes]);
                obj.q = obj.q + q_dot * obj.TA; % euler integrator % wahrscheinlich nicht richtig integriert!
                obj.q = obj.q / norm(obj.q);
            else
                % [4], eq. 14
                gyro_mean = (obj.gyro+ gyro - obj.b_w + obj.KP*ome_mes)./2;
                
                % [4], eq. 12
                Omega = [0, -gyro_mean';
                         gyro_mean,      - obj.crossproductMatrix(gyro_mean)];
                % [4], eq. 13     
                obj.q = expm(1/2 * Omega * obj.TA) * obj.q;
                % needed to calculate integral of q_dot
                obj.gyro = gyro - obj.b_w + obj.KP*ome_mes;
            end
                
        end
        
        function gamma = eulerAngles(obj)
            q = Quaternion(obj.q);
            gamma = q.quatToEuler();
        end
        
        function bias = bias(obj)
            bias = obj.b_w;
        end
        
        function m = inertiaMagneticField(obj)
            delta_d = obj.delta_d;
            delta_i = obj.delta_i;
            
            m = [cos((pi*delta_d)/180)*cos((pi*delta_i)/180);
                -sin((pi*delta_d)/180)*cos((pi*delta_i)/180);
                -sin((pi*delta_i)/180)];
        end
    end
    
    methods(Access = protected)
        %% Common functions
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
        end

        function stepImpl(obj, accelerometer, magnetometer, gyro)
            obj.update(accelerometer, magnetometer, gyro);
        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end
    end
    
    methods (Static)
        function v = Re_to_b_Quaternion(q, v)
            q_0 = q(1);
            q_1 = q(2);
            q_2 = q(3);
            q_3 = q(4);
            v = [1-2*(q_3^2+q_2^2),	2*(q_0*q_3+q_1*q_2),	2*(q_1*q_3-q_0*q_2);
                2*(q_1*q_2-q_0*q_3),	1-2*(q_3^2+q_1^2),	2*(q_2*q_3+q_0*q_1);
                2*(q_1*q_3+q_0*q_2),	2*(q_2*q_3-q_0*q_1),	1-2*(q_2^2+q_1^2)] * v;
        end
        
        % calculates the crossproduct matrix from a vector
        function M = crossproductMatrix(v)
            M = [0,-v(3),v(2); 
                 v(3),0,-v(1);
                 -v(2),v(1),0];
        end
    end
    
    methods (Access = private)
        
    end
    
    
    
    
end
