% discrete Extended kalman filter
% Version: 1.0
classdef EKF < handle % discrete Kalman filter
    
    %% Public properties
    properties (Access = public)
       sampleTime;
       Q;
       R_k;
       P;
       x;
       x_dot_h;
       G_h;
       y_h;
       H_h; 
       K; % for inspecting the Kalman Gain
    end
    
    %% Public methods
    methods (Access = public)
        % Constructor
        % varargin: 
        %             - varargin{1}: sampleTime [s] to calculate discrete measurement
        %             covariancematrix R_k from continuous time measurement
        %             matrix
        %             - varargin{2} Q: continuous time state covariance matrix
        %             - varargin{3} R: continuous time measurement
        %             covariance matrix
        %             - varargin{4} P: initial kalman covariance matrix
        %             - varargin{5} x_k: initial states
        function obj = EKF(varargin)
            ninputarguments = 5;
            if nargin < ninputarguments
               error(['Not enough input arguments: ', nargin, ' of ', ninputarguments]); 
            end
            
            obj.sampleTime = varargin{1};
            obj.Q = varargin{2};
            obj.R_k = varargin{3}; % / obj.sampleTime;
            obj.P = varargin{4};
            obj.x = varargin{5};
            
            % checking input arguments
            tempSize = size(obj.P);
            if (tempSize(1) ~= length(obj.x) || tempSize(2) ~= length(obj.x))
                error(['EKF.m: Kalman Covariance matrix has not the size of (',...
                    num2str(length(obj.x)), ',', num2str(length(obj.x)), '). It has the size',...
                    '(', num2str(tempSize(1)), ',', num2str(tempSize(2)), ')']);
            end
            
            tempSize = size(obj.Q);
            if (tempSize(1) ~= length(obj.x) || tempSize(2) ~= length(obj.x))
                error(['State covariance matrix Q has not the size of (',...
                    num2str(length(obj.x)), ',', num2str(length(obj.x)), '). It has the size',...
                    '(', num2str(tempSize(1)), ',', num2str(tempSize(2)), ')']);
            end
        end
        
        function predictorStep(obj, x_dot, G_lin)
            nStates = length(obj.x);
            %% calculating the discrete-time transition matrix PHI_k and the discrete-time error covarianze process matrix
            % using the VanLoan Method: Introduction to Random Signals and Applied Kalman Filtering
            % can be moved outside of this function, when G_lin and Q are constant
            AA = [-G_lin, obj.Q; zeros(nStates,nStates), G_lin']*obj.sampleTime;
            BB = expm(AA);
            PHI_k = BB(nStates+1:2*nStates, nStates+1: 2*nStates)';
            Q_k = PHI_k*BB(1:nStates,nStates+1:2*nStates);
            
            % calculating x_k_k1:
            % if not using euler to calculate x_k_k1 then B_n*u_n must be discretized
            % too
            x_k_k1 = obj.x + x_dot * obj.sampleTime; % see Barton perform_ekf
            obj.x = x_k_k1;
            
            % calculating P_k_k1:
            P_k_k1 = PHI_k * obj.P * PHI_k' + Q_k;
            obj.P = P_k_k1;
        end
        
        % measurement step
        % must be called after predictorStep
        function correctorStep(obj, y, y_hat, H_lin)
            %%% correction step %%%
            % Kalman Gain
            P_y_y1 = (obj.R_k + H_lin * obj.P * H_lin');
            K_k = (obj.P * H_lin') / P_y_y1; % inv(A)*b is slower than A\b
            obj.K = K_k;
            obj.P = obj.P - K_k*H_lin*obj.P;

            obj.x = obj.x + K_k*(y - y_hat);
        end
    end
end