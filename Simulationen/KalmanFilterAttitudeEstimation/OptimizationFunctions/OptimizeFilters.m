% source:
% - [1]: Mattia Giurato - https://github.com/Murmele/Attitude_Estimation


%Weight tuning
% higher means more important
w_roll = 1;
w_pitch = 1;
w_yaw = 0.5;

if filterOptimization.mahonyExplicitComplementaryFilter

    % Tuning guess
    Kacc = 12.2305; 
    Kmag = 3.6870;
    Kp = 20.95;
    Ki = 12.3743;

    x0 = [Kacc Kmag Kp Ki]';
    
    % Estimator tuning
    % U are the parameters which should be determined
    fun = @(U) mahonyExplicitComplementaryFilter_optimization(U,measurements.meas.w_b,...
                                                              measurements.meas.a_b,...
                                                              measurements.meas.m_b,...
                                                              truthdata.gamma_inertia,...
                                                              measurements.meas.p_gps,...
                                                              measurements.meas.v_gps,...
                                                              w_roll, w_pitch, w_yaw,...
                                                              delta_i, delta_d, g_mps2, TA);
    options = optimoptions('fsolve','algorithm','levenberg-marquardt');
    x = fsolve(fun,x0,options);

    %Optimal parameters
    optimizedParameters.mahonyExplicitComplementaryFilter.Kacc = x(1);
    optimizedParameters.mahonyExplicitComplementaryFilter.Kmag = x(2);
    optimizedParameters.mahonyExplicitComplementaryFilter.Kp = x(3);
    optimizedParameters.mahonyExplicitComplementaryFilter.Ki = x(4);
end

if filterOptimization.MEKFDifferentSampleTimeVarianceChange
    
    % Tuning guess
    x0 = [sigma2.a_b' sigma2.mag'];   
    
    % Estimator tuning
    % U are the parameters which should be determined
    fun = @(U) AHRSMEKFDSTVarianceChange_optimization(...
                                          U,measurements, ...
                                          truthdata.gamma_inertia,...
                                          w_roll, w_pitch, w_yaw,...
                                          delta_i, delta_d, g_mps2, TA,...
                                          sigma2);
   options = optimoptions('fsolve','algorithm','levenberg-marquardt');
   x = fsolve(fun,x0,options);
   
   %Optimal parameters
    optimizedParameters.AHRSMEKFDSTVarianceChange_optimization.sigma_2_a_b_not_available = x(1:3);
    optimizedParameters.AHRSMEKFDSTVarianceChange_optimization.sigma_2_yaw_not_available = x(4:6);
    
end