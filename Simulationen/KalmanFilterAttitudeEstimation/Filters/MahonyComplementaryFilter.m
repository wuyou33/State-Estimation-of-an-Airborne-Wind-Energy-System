% source: [1] Mahony et al. - Complementary filter design on the special orthogonal group SO(3)
% working directly on SO(3) (estimating rotation matrix)

% not able to determine gyro bias!
% no magnetic inclination/declination respected
% Version: 0.1

classdef MahonyComplementaryFilter < matlab.System
 
    % Public, tunable properties
    properties
        R; % rotation matrix
        k_est; % filter gain
        TA;
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
        function obj = MahonyComplementaryFilter(TA)
            
            % Support name-value pair arguments when constructing object
            %setProperties(obj,nargin,varargin{:})
            
            obj.TA = TA;
            obj.R = [1,0,0;
                     0,1,0;
                     0,0,1]; % no initial rotation
            % The higher this value, the higher ist cut off frequencey of
            % the low pass filter (acc and mag have more weight also in
            % higher frequencies)
            obj.k_est = 0.15;  % must not to big, because otherwise there is a too large error (don't know why this has an influence on it. For me it should'nt)     
        end
        
        function update(obj, accelerometer, magnetometer, gyro)
            % calculate Ry see [EulerAngleCalculation.wxm]
            
                % acceleration is falsified, because no centrifugalforce is
                % included
                % use accelerometer to calculate roll (phi) and pitch
                % (theta)
                a_n = sqrt(accelerometer(1)^2+accelerometer(2)^2+accelerometer(3)^2);
                ax = accelerometer(1)/a_n;
                ay = accelerometer(2)/a_n;
                az = accelerometer(3)/a_n;
                
                phi= atan2(ay,az);
                theta= asin(-ax/1);
                % use magnetometer only for yaw angle
                m_n = sqrt(magnetometer(1)^2+magnetometer(2)^2+magnetometer(3)^2);
                mx = magnetometer(1)/m_n;
                my = magnetometer(2)/m_n;
                mz = magnetometer(3)/m_n;
                psi =  atan2(sin(phi)*mz - cos(phi)*my, cos(theta)*mx + sin(phi)*sin(theta)*my + cos(phi)*sin(theta)*mz);
                
                % declination and inclination not included!!!
                
                % earth to body:
                Ry = [cos(psi)*cos(theta),	sin(psi)*cos(theta),	-sin(theta);
                      sin(phi)*cos(psi)*sin(theta)-cos(phi)*sin(psi),	sin(phi)*sin(psi)*sin(theta)+cos(phi)*cos(psi),	sin(phi)*cos(theta);
                      cos(phi)*cos(psi)*sin(theta)+sin(phi)*sin(psi),	cos(phi)*sin(psi)*sin(theta)-sin(phi)*cos(psi),	cos(phi)*cos(theta)];
                  
                Ry = Ry'; % must be body to earth!

            Omega_cross = crossproductMatrix(gyro);
            
            % [1] Blockdiagram 3/4
            blockdiagram34 = 0;
            passiveComplementaryFilter = 1;
            if blockdiagram34
                if passiveComplementaryFilter
                    % $$ (R\Omega)_\times $$
                    ROmega_cross = obj.R*Omega_cross*obj.R'; % Passive complementary filter
                else
                    ROmega_cross = Ry*Omega_cross*Ry'; % alternative. Direct complementary filter
                end
                % $$ \tilde{R} $$ in source
                R_error = obj.R'*Ry;
                % $$ \pi_a (\tilde{R}) $$
                piAR_error = 1/2*(R_error - R_error');
                % eq: 10 in [1] (in the blockdiagram there is no obj.R', maybe
                % they forgott it?)
                R_dot = (ROmega_cross + obj.k_est*obj.R*piAR_error*obj.R')*obj.R; % should work too, but don't understand why not
                obj.R = obj.R + R_dot * obj.TA; % euler integrator
            else
                % Optimized. [1] Blockdiagram 5
                % eq: 11 in [1]
                R_error = obj.R'*Ry;
                piAR_error = 1/2*(R_error - R_error');
                R_dot = obj.R*(Omega_cross + obj.k_est*piAR_error);
                obj.R = obj.R + R_dot * obj.TA; % euler integrator
            end
            
            % Es ist aufgefallen, dass der Parameter obj.k_est
            % unterschiedlich viel Auswirkung auf die Winkelfehler hat,
            % weshalb man auch versuchen könnte, obj.k_est mehrdimensional
            % zu machen und dann die Werte bestimmen um jeden Winkel
            % einzeln zu optimieren?
                
        end
        
        function gamma = eulerAngles(obj)
            gamma = rotationMatrixToEulerAngles(obj.R');
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
    
    methods (Access = private)
        
    end
    
    
    
    
end
