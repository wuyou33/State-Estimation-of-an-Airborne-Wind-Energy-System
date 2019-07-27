function RMS = mahonyExplicitComplementaryFilter_optimization(U, w_b_meas,...
                                                              a_b_meas,...
                                                              m_b_meas,...
                                                              attitude,...
                                                              position, ...
                                                              velocity, ...
                                                              w_roll, w_pitch, w_yaw, delta_i, delta_d, g_mps2, TA)
    % MahonyExplicitComplementaryFilter
    q0 = [1;0;0;0];
    b_w0 = [0;0;0];
    mahonyExplicitComplementaryFilter = MahonyExplicitComplementaryFilter(TA, delta_i, delta_d, g_mps2, q0, b_w0);
    
    mahonyExplicitComplementaryFilter.Kacc = U(1);
    mahonyExplicitComplementaryFilter.Kmag = U(2);
    mahonyExplicitComplementaryFilter.KP = U(3);
    mahonyExplicitComplementaryFilter.KI = U(4);

    b_w_mahonyExplicitComplementaryFilter = zeros(3, length(attitude));                                                    
    att_est = zeros(3, length(attitude));
    for i = 1: length(attitude)                                                 
        % MahonyExplicitComplementaryFilterQuaterion
        mahonyExplicitComplementaryFilter.update(a_b_meas(:,i), m_b_meas(:,i), w_b_meas(:,i), position(:,i), velocity(:,i));
        att_est(:, i) = rad2deg(mahonyExplicitComplementaryFilter.eulerAngles());
        b_w_mahonyExplicitComplementaryFilter(:, i) = mahonyExplicitComplementaryFilter.bias(); 
    end


    error = attitude - att_est;
    RMS = (w_roll*rms(error(1,:)) + w_pitch*rms(error(2,:)) + w_yaw*rms(error(3,:)))/(w_roll + w_pitch + w_yaw);
end