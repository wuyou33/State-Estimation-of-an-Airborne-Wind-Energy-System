function RMS = AHRSMEKFDSTVarianceChange_optimization(U, measurements,...
                                                              attitude,...
                                                              w_roll, w_pitch, w_yaw, ...
                                                              delta_i, delta_d, ...
                                                              g_mps2, TA, ...
                                                              sigma2)
    %% Initial rotation
    a_n = vecnorm(measurements.meas.a_b);
    a_n = a_n(1);
    % es darf nicht die norm gemacht werden, sonst wird theta immer sehr
    % klein sein!!!
    ax = measurements.meas.a_b(1,1); % ./a_n;
    ay = measurements.meas.a_b(2,1); % ./a_n;
    az = measurements.meas.a_b(3,1); % ./a_n;

    phi_init= atan2(ay,az);
    theta_init= asin(-ax/g_mps2); % hier!
    m_b_init = Rmag_to_e(measurements.meas.m_b(:,1), delta_d, delta_i);
    
    m_n = vecnorm(m_b_init);
    m_n = m_n(1);
    mx = m_b_init(1,1)./m_n;
    my = m_b_init(2,1)./m_n;
    mz = m_b_init(3,1)./m_n;
    psi_init =  atan((sin(phi_init).*mz - cos(phi_init).*my)/ (cos(theta_init).*mx + sin(phi_init).*sin(theta_init).*my + cos(phi_init).*sin(theta_init).*mz));
    
                                                          
    %% Init filter                                                      
    q0 = eulerAnglesToQuaternion([phi_init; theta_init; psi_init]);
    b_w0 = [0;0;0];
    P0 = [0.3, 0, 0, 0, 0, 0;
          0, 1.5, 0, 0, 0, 0;
          0, 0, 0.3, 0, 0, 0;
          0, 0, 0, 5, 0, 0; % bias to a high variance. In real it will never so high. But above I use as bias 3rad/s
          0, 0, 0, 0, 5, 0;
          0, 0, 0, 0, 0, 5];
    ahrsMEKFDifferentSampleTimeVarianceChange = AHRSMEKFDifferentSampleTimeVarianceChange(TA, delta_i, delta_d, g_mps2, q0, b_w0, P0, sigma2.w_b, sigma2.a_b, sigma2.b_w, sigma2.mag);

    clear q0 b_w0 P0;
    
    ahrsMEKFDifferentSampleTimeVarianceChange.sigma_2_a_b_not_available = U(1:3);
    ahrsMEKFDifferentSampleTimeVarianceChange.sigma_2_yaw_not_available = U(4:6);

    b_w_ahrsMEKFDifferentSampleTimeVarianceChange = zeros(3, length(attitude));                                                    
    att_est = zeros(3, length(attitude));
    for i = 1: length(attitude)
        u = [measurements.meas.w_b(:,i); measurements.meas.v_gps(:,i)];
        y = [measurements.meas.a_b(:,i); measurements.meas.m_b(:,i)];
        ahrsMEKFDifferentSampleTimeVarianceChange.performEKF(u, y, measurements.available.acc(i), measurements.available.mag(i));
        q_hat = ahrsMEKFDifferentSampleTimeVarianceChange.attitude();
        q = Quaternion(q_hat);
        att_est(:, i) = rad2deg(q.quatToEuler());
        b_w_ahrsMEKFDifferentSampleTimeVarianceChange(:, i) = ahrsMEKFDifferentSampleTimeVarianceChange.gyroBias(); 
    end

    error = attitude - att_est;
    RMS = (w_roll*rms(error(1,:)) + w_pitch*rms(error(2,:)) + w_yaw*rms(error(3,:)))/(w_roll + w_pitch + w_yaw);
end