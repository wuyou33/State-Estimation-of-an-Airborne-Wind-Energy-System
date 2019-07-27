if ~simulation
    display('No simulation started');
    return;
end

%% Evaluating the results
if plotting
    %close all;
    display('Data evaluation');
    
    indexMaxError = find(time >= maxErrorTime);
    
    %% figure 2: attitude
    [phi_log, theta_log, psi_log] = calculateEulerAnglesFromAccelerometerAndMagnetometer(measurements.meas.a_b, measurements.meas.m_b, delta_d, delta_i);
    try
    [phi_log2, theta_log2, psi_log2] = calculateEulerAnglesFromAccelerometerAndMagnetometer(measurements.meas.a_b, m_b_meas_woutOffsetCompensation, delta_d, delta_i);
    catch
    end
    
    if plotDifferenceAngles
        magnetometerComensationAngles = figure('Name', 'Magnetometer compensation effects');
        figure(magnetometerComensationAngles);
        subplot(2,1,1);
        title('Calculated \psi from different magnetometer values [°]');
        xlabel('Time [s]');
        hold on;
        try
        plot(time, measurements.meas.poti_deg, 'DisplayName','\psi_{potentiometer} (real)'); % directly from the sensordata without filters
        catch
        end
        plot(time, psi_log(1,:), 'DisplayName','\psi_{calc} with hard and soft iron compensation'); % directly from the sensordata without filters
        try
            plot(time, psi_log2(1,:), 'DisplayName','\psi_{calc} without compensation'); % directly from the sensordata without filters
        catch
        end
        
        try
            [phi_log3, theta_log3, psi_log3] = calculateEulerAnglesFromAccelerometerAndMagnetometer(measurements.meas.a_b, measurements.meas.m_b_HI_comp, delta_d, delta_i);
            plot(time, psi_log3(1,:), 'DisplayName','\psi_{calc} with hard iron only compensation'); % directly from the sensordata without filters
            
        catch
        end
        try
        rmsd_HI = RMSD(measurements.meas.poti_deg, psi_log3); % with hard iron compensation only
        rmsd_HISI = RMSD(measurements.meas.poti_deg, psi_log); % with hard and soft iron compensation
        rmsd_Psi = RMSD(measurements.meas.poti_deg, psi_log2); % without offset compensation
        
        rmsd_plot = {['With hard iron compensation only: rmsd = ', num2str(rmsd_HI), ', max = ', num2str(max(abs(measurements.meas.poti_deg - psi_log3)))];
                     ['With hard and soft iron compensation: rmsd = ', num2str(rmsd_HISI), ', max = ', num2str(max(abs(measurements.meas.poti_deg - psi_log)))];
                     ['Without offset compensation: rmsd = ', num2str(rmsd_Psi), ', max = ', num2str(max(abs(measurements.meas.poti_deg - psi_log2)))]};
        
        text(30,80,rmsd_plot, 'clipping','on');
        catch
        end
        
        hold off;
        
        subplot(2,1,2);
        title('Error compared to potentiometer angle [°]');
        hold on;
        try
            plot(time, psi_log(1,:) - measurements.meas.poti_deg, 'DisplayName','\psi_{calc} with hard and soft iron compensation'); % directly from the sensordata without filters
        catch
        end
        try
            plot(time, psi_log2(1,:) - measurements.meas.poti_deg, 'DisplayName','\psi_{calc} without compensation'); % directly from the sensordata without filters
        catch
        end 
        try
            [phi_log3, theta_log3, psi_log3] = calculateEulerAnglesFromAccelerometerAndMagnetometer(measurements.meas.a_b, measurements.meas.m_b_HI_comp, delta_d, delta_i);
            plot(time, psi_log3(1,:) - measurements.meas.poti_deg, 'DisplayName','\psi_{calc} with hard iron only compensation'); % directly from the sensordata without filters
            
        catch
        end
        plot(time, zeros(length(time),1), '-k', 'HandleVisibility','off');
        hold off;
        
        clear rmsd_HI rmsd_HISI rmsd_Psi rmsd_plot;
        
    end
    
    % Error between real attitude and estimated attitude
    skip = 0;
    if strcmp(angleErrorPlot,'eulerEKF') && filterSimulation.eulerEKF
        errorPlot = gamma_i_ahrsEuler;
    elseif strcmp(angleErrorPlot,'eulerEKFHIComp') && filterSimulation.eulerEKFHIComp
        errorPlot = gamma_i_ahrsEulerHIComp;
    elseif strcmp(angleErrorPlot,'eulerEKFWithoutGyroBias') && filterSimulation.eulerEKFWoutGyroBias
        errorPlot = gamma_i_ahrsEulerWoutGyroBias;
    elseif strcmp(angleErrorPlot,'eulerEKFDifferentSampleTime') && filterSimulation.eulerEKFDifferentSampleTime
        errorPlot = gamma_i_ahrsEulerDifferentSampleTime;
    elseif strcmp(angleErrorPlot, 'quaternionEKF') && filterSimulation.quaternionEKF
        errorPlot = gamma_i_ahrsQuaternion;
    elseif strcmp(angleErrorPlot, 'quaternionEKFDifferentSampleTime') && filterSimulation.quaternionEKFDifferentSampleTime
        errorPlot = gamma_i_ahrsQuaternionDifferentSampleTime;
    elseif strcmp(angleErrorPlot, 'MEKF') && filterSimulation.MEKF
        errorPlot = gamma_i_ahrsMEKF;
    elseif strcmp(angleErrorPlot, 'MEKFDifferentSampleTime') && filterSimulation.MEKFDifferentSampleTime
        errorPlot = gamma_i_ahrsMEKFDifferentSampleTime;
    elseif strcmp(angleErrorPlot, 'MEKFDifferentSampleTimeVarianceChange') && filterSimulation.MEKFDifferentSampleTimeVarianceChange
        errorPlot = gamma_i_ahrsMEKFDifferentSampleTimeVarianceChange;
    elseif strcmp(angleErrorPlot, 'MEKFDifferentSampleTimeMagnetometer') && filterSimulation.MEKFDifferentSampleTimeMagnetometer
        errorPlot = gamma_i_ahrsMEKFDifferentSampleTimeMagnetometer;
    elseif strcmp(angleErrorPlot, 'MEKFDifferentSampleTimeMagnetometer2') && filterSimulation.MEKFDifferentSampleTimeMagnetometer2
        errorPlot = gamma_i_ahrsMEKFDifferentSampleTimeMagnetometer2;
    elseif strcmp(angleErrorPlot, 'MEKFwithoutCentripetalCompensation') && filterSimulation.MEKFwithoutCentripetalCompensation
        errorPlot = gamma_i_ahrsMEKFwithoutCentripetalCompensation;
    elseif strcmp(angleErrorPlot, 'mahonyComplementaryFilter') && filterSimulation.mahonyComplementaryFilter
        errorPlot = gamma_i_mahonyComplementaryFilter;
    elseif strcmp(angleErrorPlot, 'mahonyExplicitComplementaryFilter') && filterSimulation.mahonyExplicitComplementaryFilter
        errorPlot = gamma_i_mahonyExplicitComplementaryFilter;
    elseif strcmp(angleErrorPlot, 'INSEulerEKF') && filterSimulation.INSEulerEKF
        errorPlot = gamma_i_insEulerEKF;
    elseif strcmp(angleErrorPlot, 'INSEulerEKFDST') && filterSimulation.INSEulerEKFDST
        errorPlot = gamma_i_insEulerEKFDST;
    elseif strcmp(angleErrorPlot, 'INSMEKF') && filterSimulation.INSMEKF
        errorPlot = gamma_i_insMEKF;
    elseif strcmp(angleErrorPlot, 'INSMEKFDST') && filterSimulation.INSMEKFDST
        errorPlot = gamma_i_INSMEKFDST;
    elseif strcmp(angleErrorPlot, 'INSMEKFDSTVarianceChange') && filterSimulation.INSMEKFDSTVarianceChange
        errorPlot = gamma_i_insMEKF;
    elseif strcmp(angleErrorPlot, 'INSMEKFDST_tetherSag') && filterSimulation.INSMEKFDST_tetherSag
        errorPlot = gamma_i_INSMEKFDST_tetherSag;    
    else
        warning(['ErrorPlot "' angleErrorPlot '"is not available']);
        skip = 1;
    end
    
    f = figure('Name','Attitude','NumberTitle','off');
    figure(f); % must be set, otherwise the plot is not in a new window
    movegui(f,'northwest');
    if skip
        subplot(1,3,1);
    else
        subplot(2,3,1);
    end
    title('\phi [\si{\degree}]');
    xlabel('time [s]');
    hold('on');
    plot(time, truthdata.gamma_inertia(1,:), 'DisplayName','real');
    try
        plot(time, truthdata.multikopterEstimatioGammaInertia(1, :), 'DisplayName', 'Multikopter estimation')
    catch
    end
    %plot(time, phi_log(1,:), 'DisplayName','phi_{calc} with hard iron compensation'); % directly from the sensordata without filters
    try
        plot(time, phi_log2(1,:), 'DisplayName','phi_{calc} without hard iron compensation'); % directly from the sensordata without filters
    catch
    end
    if filterSimulation.eulerEKF
      plot(time, gamma_i_ahrsEuler(1,:),'-', 'DisplayName','ahrsEuler');
    end
    if filterSimulation.eulerEKFHIComp
      plot(time, gamma_i_ahrsEulerHIComp(1,:),'-', 'DisplayName','ahrsEulerHIComp');
    end
    if filterSimulation.eulerEKFWoutGyroBias
      plot(time, gamma_i_ahrsEulerWoutGyroBias(1,:),'-', 'DisplayName','ahrsEulerWoutGyroBias');
    end
    if filterSimulation.eulerEKFDifferentSampleTime
      plot(time, gamma_i_ahrsEulerDifferentSampleTime(1,:),'-', 'DisplayName','ahrsEulerDST');
    end
    if filterSimulation.quaternionEKF
      plot(time, gamma_i_ahrsQuaternion(1,:),'-', 'DisplayName','ahrsQuaternion');
    end
    if filterSimulation.quaternionEKFDifferentSampleTime
      plot(time, gamma_i_ahrsQuaternionDifferentSampleTime(1,:),'-', 'DisplayName','ahrsQuaternionDST');
    end
    if filterSimulation.MEKF
      plot(time, gamma_i_ahrsMEKF(1,:),'-', 'DisplayName','ahrsMEKF');
    end
    if filterSimulation.MEKFDifferentSampleTime
      plot(time, gamma_i_ahrsMEKFDifferentSampleTime(1,:),'-', 'DisplayName','ahrsMEKFDST');
    end
    if filterSimulation.MEKFDifferentSampleTimeVarianceChange
      plot(time, gamma_i_ahrsMEKFDifferentSampleTimeVarianceChange(1,:),'-', 'DisplayName','ahrsMEKFDSTVarChng');
    end
    if filterSimulation.MEKFDifferentSampleTimeMagnetometer
      plot(time, gamma_i_ahrsMEKFDifferentSampleTimeMagnetometer(1,:),'-', 'DisplayName','ahrsMEKFDST_Mag');
    end
    if filterSimulation.MEKFDifferentSampleTimeMagnetometer2
      plot(time, gamma_i_ahrsMEKFDifferentSampleTimeMagnetometer2(1,:),'-', 'DisplayName','ahrsMEKFDST_Mag2');
    end
    if filterSimulation.MEKFwithoutCentripetalCompensation
      plot(time, gamma_i_ahrsMEKFwithoutCentripetalCompensation(1,:),'-', 'DisplayName','ahrsMEKFw/outCC');
    end
    if filterSimulation.mahonyComplementaryFilter
      plot(time, gamma_i_mahonyComplementaryFilter(1,:),'-', 'DisplayName','mahonyComplementaryFilter');
    end
    if filterSimulation.mahonyExplicitComplementaryFilter
      plot(time, gamma_i_mahonyExplicitComplementaryFilter(1,:),'-', 'DisplayName',['mahony explicit' newline 'complementary filter']);
    end
    if filterSimulation.INSEulerEKF
      plot(time, gamma_i_insEulerEKF(1,:),'-', 'DisplayName','INSEulerEKF');
    end
    if filterSimulation.INSEulerEKFDST
      plot(time, gamma_i_insEulerEKFDST(1,:),'-', 'DisplayName','INSEulerEKFDST');
    end
    if filterSimulation.INSMEKF
      plot(time, gamma_i_insMEKF(1,:),'-', 'DisplayName','INSMEKF');
    end
    if filterSimulation.INSMEKFDST
      plot(time, gamma_i_INSMEKFDST(1,:),'-', 'DisplayName','INSMEKFDST');
    end
    if filterSimulation.INSMEKFDSTVarianceChange
      plot(time, gamma_i_INSMEKFDSTVarianceChange(1,:),'-', 'DisplayName',['INSMEKF', newline, 'Interpolate']);
    end
    if filterSimulation.INSMEKFDST_tetherSag
      plot(time, gamma_i_INSMEKFDST_tetherSag(1,:),'-', 'DisplayName','INSMEKFDST_{tetherSag}');
    end
    hold('off');
    %legend;
    if skip
        subplot(1,3,2);
    else
        subplot(2,3,2);
    end
    title('\theta [\si{\degree}]');
    xlabel('time [s]');
    hold('on');
    plot(time, truthdata.gamma_inertia(2,:), 'DisplayName','real');
    try
        plot(time, truthdata.multikopterEstimatioGammaInertia(2, :), 'DisplayName', 'Multikopter estimation')
    catch
    end
    %plot(time, theta_log(1,:), 'DisplayName','calc  with hard iron compensation'); % directly from the sensordata without filters
    try
        plot(time, theta_log2(1,:), 'DisplayName','calc without hard iron compensation'); % directly from the sensordata without filters    
    catch
    end
    if filterSimulation.eulerEKF
        plot(time, gamma_i_ahrsEuler(2,:),'-', 'DisplayName','ahrsEuler');
    end
    if filterSimulation.eulerEKFHIComp
        plot(time, gamma_i_ahrsEulerHIComp(2,:),'-', 'DisplayName','ahrsEulerHIComp');
    end
    if filterSimulation.eulerEKFWoutGyroBias
        plot(time, gamma_i_ahrsEulerWoutGyroBias(2,:),'-', 'DisplayName','ahrsEulerWoutGyroBias');
    end
    if filterSimulation.eulerEKFDifferentSampleTime
        plot(time, gamma_i_ahrsEulerDifferentSampleTime(2,:),'-', 'DisplayName','ahrsEulerDST');
    end
    if filterSimulation.quaternionEKF
        plot(time, gamma_i_ahrsQuaternion(2,:),'-', 'DisplayName','ahrsQuaternion');
    end
    if filterSimulation.quaternionEKFDifferentSampleTime
        plot(time, gamma_i_ahrsQuaternionDifferentSampleTime(2,:),'-', 'DisplayName','ahrsQuaternionDST');
    end
    if filterSimulation.MEKF
        plot(time, gamma_i_ahrsMEKF(2,:),'-', 'DisplayName','ahrsMEKF');
    end
    if filterSimulation.MEKFDifferentSampleTime
        plot(time, gamma_i_ahrsMEKFDifferentSampleTime(2,:),'-', 'DisplayName','ahrsMEKFDST');
    end
    if filterSimulation.MEKFDifferentSampleTimeVarianceChange
        plot(time, gamma_i_ahrsMEKFDifferentSampleTimeVarianceChange(2,:),'-', 'DisplayName','ahrsMEKFDSTVarChng');
    end
    if filterSimulation.MEKFDifferentSampleTimeMagnetometer
        plot(time, gamma_i_ahrsMEKFDifferentSampleTimeMagnetometer(2,:),'-', 'DisplayName','ahrsMEKFDST_Mag');
    end
    if filterSimulation.MEKFDifferentSampleTimeMagnetometer2
        plot(time, gamma_i_ahrsMEKFDifferentSampleTimeMagnetometer2(2,:),'-', 'DisplayName','ahrsMEKFDST_Mag2');
    end
    if filterSimulation.MEKFwithoutCentripetalCompensation
        plot(time, gamma_i_ahrsMEKFwithoutCentripetalCompensation(2,:),'-', 'DisplayName','ahrsMEKFw/outCC');
    end
    if filterSimulation.mahonyComplementaryFilter
        plot(time, gamma_i_mahonyComplementaryFilter(2,:),'-', 'DisplayName','mahonyComplementaryFilter');
    end
    if filterSimulation.mahonyExplicitComplementaryFilter
      plot(time, gamma_i_mahonyExplicitComplementaryFilter(2,:),'-', 'DisplayName',['mahony explicit' newline 'complementary filter']);
    end
    if filterSimulation.INSEulerEKF
      plot(time, gamma_i_insEulerEKF(2,:),'-', 'DisplayName','INSEulerEKF');
    end
    if filterSimulation.INSEulerEKFDST
      plot(time, gamma_i_insEulerEKFDST(2,:),'-', 'DisplayName','INSEulerEKFDST');
    end
    if filterSimulation.INSMEKF
      plot(time, gamma_i_insMEKF(2,:),'-', 'DisplayName','INSMEKF');
    end
    if filterSimulation.INSMEKFDST
      plot(time, gamma_i_INSMEKFDST(2,:),'-', 'DisplayName','INSMEKFDST');
    end
    if filterSimulation.INSMEKFDSTVarianceChange
      plot(time, gamma_i_INSMEKFDSTVarianceChange(2,:),'-', 'DisplayName',['INSMEKF', newline, 'Interpolate']);
    end
    if filterSimulation.INSMEKFDST_tetherSag
      plot(time, gamma_i_INSMEKFDST_tetherSag(2,:),'-', 'DisplayName','INSMEKFDST_tetherSag');
    end
    hold('off');
    %legend;
    
    if skip
        subplot(1,3,3);
    else
        subplot(2,3,3);
    end
    title('\psi [\si{\degree}]');
    xlabel('time [s]');
    hold('on');
    plot(time, truthdata.gamma_inertia(3,:), 'DisplayName','real');
    try
        plot(time, truthdata.multikopterEstimatioGammaInertia(3, :), 'DisplayName', 'Multikopter estimation')
    catch
    end
    %plot(time, psi_log(1,:), 'DisplayName','psi_calc  with hard iron compensation'); % directly from the sensordata without filters
    try
        plot(time, psi_log2(1,:), 'DisplayName','calc without hard iron compensation'); % directly from the sensordata without filters
    catch
    end
    if filterSimulation.eulerEKF
        plot(time, gamma_i_ahrsEuler(3,:),'-', 'DisplayName','ahrsEuler');
    end
    if filterSimulation.eulerEKFHIComp
        plot(time, gamma_i_ahrsEulerHIComp(3,:),'-', 'DisplayName','ahrsEulerHIComp');
    end
    if filterSimulation.eulerEKFWoutGyroBias
        plot(time, gamma_i_ahrsEulerWoutGyroBias(3,:),'-', 'DisplayName','ahrsEulerWoutGyroBias');
    end
    if filterSimulation.eulerEKFDifferentSampleTime
        plot(time, gamma_i_ahrsEulerDifferentSampleTime(3,:),'-', 'DisplayName','ahrsEulerDST');
    end
    if filterSimulation.quaternionEKF
        plot(time, gamma_i_ahrsQuaternion(3,:),'-', 'DisplayName','ahrsQuaternion');
    end
    if filterSimulation.quaternionEKFDifferentSampleTime
        plot(time, gamma_i_ahrsQuaternionDifferentSampleTime(3,:),'-', 'DisplayName','ahrsQuaternionDST');
    end
    if filterSimulation.MEKF
        plot(time, gamma_i_ahrsMEKF(3,:),'-', 'DisplayName','ahrsMEKF');
    end
    if filterSimulation.MEKFDifferentSampleTime
        plot(time, gamma_i_ahrsMEKFDifferentSampleTime(3,:),'-', 'DisplayName','ahrsMEKFDST');
    end
    if filterSimulation.MEKFDifferentSampleTimeVarianceChange
        plot(time, gamma_i_ahrsMEKFDifferentSampleTimeVarianceChange(3,:),'-', 'DisplayName','ahrsMEKFDSTVarChng');
    end
    if filterSimulation.MEKFDifferentSampleTimeMagnetometer
        plot(time, gamma_i_ahrsMEKFDifferentSampleTimeMagnetometer(3,:),'-', 'DisplayName','ahrsMEKFDST_Mag');
    end
    if filterSimulation.MEKFDifferentSampleTimeMagnetometer2
        plot(time, gamma_i_ahrsMEKFDifferentSampleTimeMagnetometer2(3,:),'-', 'DisplayName','ahrsMEKFDST_Mag2');
    end
    if filterSimulation.MEKFwithoutCentripetalCompensation
        plot(time, gamma_i_ahrsMEKFwithoutCentripetalCompensation(3,:),'-', 'DisplayName','ahrsMEKFw/outCC');
    end
    if filterSimulation.mahonyComplementaryFilter
        plot(time, gamma_i_mahonyComplementaryFilter(3,:),'-', 'DisplayName','mahonyComplementaryFilter');
    end
    if filterSimulation.mahonyExplicitComplementaryFilter
      plot(time, gamma_i_mahonyExplicitComplementaryFilter(3,:),'-', 'DisplayName',['mahony explicit' newline 'complementary filter']);
    end
    if filterSimulation.INSEulerEKF
      plot(time, gamma_i_insEulerEKF(3,:),'-', 'DisplayName','insEulerEKF');
    end
    if filterSimulation.INSEulerEKFDST
      plot(time, gamma_i_insEulerEKFDST(3,:),'-', 'DisplayName','insEulerEKFDST');
    end
    if filterSimulation.INSMEKF
      plot(time, gamma_i_insMEKF(3,:),'-', 'DisplayName','insMEKF');
    end
    if filterSimulation.INSMEKFDST
      plot(time, gamma_i_INSMEKFDST(3,:),'-', 'DisplayName','INSMEKFDST');
    end
    if filterSimulation.INSMEKFDSTVarianceChange
      plot(time, gamma_i_INSMEKFDSTVarianceChange(3,:),'-', 'DisplayName',['INSMEKF', newline, 'Interpolate']);
    end
    if filterSimulation.INSMEKFDST_tetherSag
      plot(time, gamma_i_INSMEKFDST_tetherSag(3,:),'-', 'DisplayName','INSMEKFDST_{tetherSag}');
    end
    hold('off')
    
    if ~skip
        subplot(2,3,[4,5,6]);
        xlabel('time [s]');
        hold('on');
        title([angleErrorPlot, ': Errors [\si{\degree}]']);
        plot(time, calculateAngleDifferenceDeg(truthdata.gamma_inertia(1,:), errorPlot(1,:)), 'DisplayName',['\phi euler']);
        plot(time, calculateAngleDifferenceDeg(truthdata.gamma_inertia(2,:),errorPlot(2,:)), 'DisplayName',['\theta euler']);
        plot(time, calculateAngleDifferenceDeg(truthdata.gamma_inertia(3,:),errorPlot(3,:)), 'DisplayName',['\psi euler']);
        hold('off');
    end
    
    %% figure 3 Bias estimation
    skip = 0;
    if strcmp(gyroBiasPlot, 'eulerEKF') && filterSimulation.eulerEKF
        biasPlot = b_w_ahrsEuler;
    elseif strcmp(gyroBiasPlot, 'eulerEKFHIComp') && filterSimulation.eulerEKFHIComp
        biasPlot = b_w_ahrsEulerHIComp;
    elseif strcmp(gyroBiasPlot, 'eulerEKFWoutGyroBias') && filterSimulation.eulerEKFWoutGyroBias
        biasPlot = b_w_ahrsEulerWoutGyroBias;
    elseif strcmp(gyroBiasPlot, 'eulerEKFDifferentSampleTime') && filterSimulation.eulerEKFDifferentSampleTime
        biasPlot = b_w_ahrsEulerDifferentSampleTime;
    elseif strcmp(gyroBiasPlot, 'quaternionEKF') && filterSimulation.quaternionEKF
        biasPlot = b_w_ahrsQuaternion;
    elseif strcmp(gyroBiasPlot, 'quaternionEKFDifferentSampleTime') && filterSimulation.quaternionEKFDifferentSampleTime
        biasPlot = b_w_ahrsQuaternionDifferentSampleTime;
    elseif strcmp(gyroBiasPlot, 'MEKF') && filterSimulation.MEKF
        biasPlot = b_w_ahrsMEKF;
    elseif strcmp(gyroBiasPlot, 'MEKFDifferentSampleTime') && filterSimulation.MEKFDifferentSampleTime
        biasPlot = b_w_ahrsMEKFDifferentSampleTime;
    elseif strcmp(gyroBiasPlot, 'MEKFDifferentSampleTimeVarianceChange') && filterSimulation.MEKFDifferentSampleTimeVarianceChange
        biasPlot = b_w_ahrsMEKFDifferentSampleTimeVarianceChange;
    elseif strcmp(gyroBiasPlot, 'MEKFDifferentSampleTimeMagnetometer') && filterSimulation.MEKFDifferentSampleTimeMagnetometer
        biasPlot = b_w_ahrsMEKFDifferentSampleTimeMagnetometer;
    elseif strcmp(gyroBiasPlot, 'MEKFDifferentSampleTimeMagnetometer2') && filterSimulation.MEKFDifferentSampleTimeMagnetometer2
        biasPlot = b_w_ahrsMEKFDifferentSampleTimeMagnetometer2;
    elseif strcmp(gyroBiasPlot, 'MEKFwithoutCentripetalCompensation') && filterSimulation.MEKFwithoutCentripetalCompensation
        biasPlot = b_w_ahrsMEKFwithoutCentripetalCompensation;
    elseif strcmp(gyroBiasPlot, 'mahonyComplementaryFilter') && filterSimulation.mahonyComplementaryFilter
        skip = 1;
        warning('Bias estimation for mahonyComplementaryFilter Not implemented');
    elseif strcmp(gyroBiasPlot, 'mahonyExplicitComplementaryFilter') && filterSimulation.mahonyExplicitComplementaryFilter
        biasPlot = b_w_mahonyExplicitComplementaryFilter;
    elseif strcmp(gyroBiasPlot, 'INSEulerEKF') && filterSimulation.INSEulerEKF
        biasPlot = b_w_insEulerEKF;
    elseif strcmp(gyroBiasPlot, 'INSEulerEKFDST') && filterSimulation.INSEulerEKFDST
        biasPlot = b_w_insEulerEKFDST;
    elseif strcmp(gyroBiasPlot, 'INSMEKF') && filterSimulation.INSMEKF
        biasPlot = b_w_insMEKF;
    elseif strcmp(gyroBiasPlot, 'INSMEKFDST') && filterSimulation.INSMEKFDST
        biasPlot = b_w_INSMEKFDST;
    elseif strcmp(angleErrorPlot, 'INSMEKFDSTVarianceChange') && filterSimulation.INSMEKFDSTVarianceChange
        biasPlot = b_w_INSMEKFDSTVarianceChange;
    elseif strcmp(gyroBiasPlot, 'INSMEKFDST_tetherSag') && filterSimulation.INSMEKFDST_tetherSag
        biasPlot = b_w_INSMEKFDST_tetherSag;
    else
        skip = 1;
        warning(['BiasPlot "' gyroBiasPlot '"is not available or the simulation was not enabled (See entries in filterSimulation)']);
    end
    if ~skip
        f = figure('Name','Bias estimation','NumberTitle','off');
        figure(f); % otherwise it did not work correctly
        movegui(f,'north');
        hold('on');
        title([gyroBiasPlot, ': Estimated bias [rad/s]']);
        xlabel('time [s]');
        plot(time, b_w(1)*ones(length(time),1),':+r', 'DisplayName', ['Bias phi']);
        plot(time, b_w(2)*ones(length(time),1),':og', 'DisplayName', ['Bias theta']);
        plot(time, b_w(3)*ones(length(time),1),':xb', 'DisplayName', ['Bias psi']);
        
        plot(time, biasPlot(1, :),'--m', 'LineWidth', 2, 'DisplayName', 'Est: Bias phi'); % Estimated
        plot(time, biasPlot(2, :),'--y', 'LineWidth', 2, 'DisplayName', 'Est: Bias theta');
        plot(time, biasPlot(3, :),'--c', 'LineWidth', 2, 'DisplayName', 'Est: Bias psi');
    end
    hold('off'); 
    
    clear biasPlot;
    
    %% figure 4
    % Position and velocity
    skip = 1;
    if strcmp(posVelPlot, 'INSEulerEKF') && filterSimulation.INSEulerEKF
        skip = 0;
        posPlot = p_insEulerEKF;
        velPlot = v_insEulerEKF;
    elseif strcmp(posVelPlot, 'INSEulerEKFDST') && filterSimulation.INSEulerEKFDST
        skip = 0;
        posPlot = p_insEulerEKF;
        velPlot = v_insEulerEKF;
    elseif strcmp(posVelPlot, 'INSMEKF') && filterSimulation.INSMEKF
        skip = 0;
        posPlot = p_insMEKF;
        velPlot = v_insMEKF;
    elseif strcmp(posVelPlot, 'INSMEKFDST') && filterSimulation.INSMEKFDST
        skip = 0;
        posPlot = p_INSMEKFDST;
        velPlot = v_INSMEKFDST;
    elseif strcmp(posVelPlot, 'INSMEKFDSTVarianceChange') && filterSimulation.INSMEKFDSTVarianceChange
        skip = 0;
        posPlot = p_INSMEKFDSTVarianceChange;
        velPlot = v_INSMEKFDSTVarianceChange;
    elseif strcmp(posVelPlot, 'INSMEKFDST_tetherSag') && filterSimulation.INSMEKFDST_tetherSag
        skip = 0;
        posPlot = p_INSMEKFDST_tetherSag;
        velPlot = v_INSMEKFDST_tetherSag;
    else
        skip = 1;
        warning(['posVelPlot "' posVelPlot '"is not available or the simulation was not enabled (See entries in filterSimulation)']);   
    end
    
    if exist('truthdata.p_i')
        f = figure('Name','Position[m] and velocity[m/s]','NumberTitle','off');
        figure(f);
        movegui(f,'northeast');
        subplot(3,1,1);
        title('Position 3D');
        xlabel('x [m]');
        ylabel('y [m]');
        zlabel('z [m]');
        hold on;
        plot3(truthdata.p_i(1,:), truthdata.p_i(2,:), truthdata.p_i(3,:), 'DisplayName', 'Real position');
        if ~skip
            title(['Position 3D: ' posVelPlot]);
            plot3(posPlot(1,:), posPlot(2,:), posPlot(3,:), 'DisplayName', 'Est. position');
        end
        hold off;
        subplot(3,1,2);
        title('Position [m]');
        hold('on');
        xlabel('time [s]');
        plot(time, truthdata.p_i(1,:), 'DisplayName','x');
        plot(time, truthdata.p_i(2,:), 'DisplayName','y');
        plot(time, truthdata.p_i(3,:), 'DisplayName','z');
        if ~skip
            title(['Position [m]: ' posVelPlot]);
            plot(time, posPlot(1,:), 'DisplayName','est. x');
            plot(time, posPlot(2,:), 'DisplayName','est. y');
            plot(time, posPlot(3,:), 'DisplayName','est. z');
        end
        hold('off');
        subplot(3,1,3);
        title('Inertia velocity [m/s]');
        xlabel('time [s]');
        hold('on');
        plot(time, truthdata.v_i(1, :), 'DisplayName', 'x');
        plot(time, truthdata.v_i(2, :), 'DisplayName', 'y');
        plot(time, truthdata.v_i(3, :), 'DisplayName', 'z');
        if ~skip
            title(['Inertia velocity [m/s]: ' posVelPlot]);
            plot(time, velPlot(1,:), 'DisplayName','est. x');
            plot(time, velPlot(2,:), 'DisplayName','est. y');
            plot(time, velPlot(3,:), 'DisplayName','est. z');
        end
        hold('off');
    end
end

%%
if evaluation
    %% Table to show square Error
    colNames = ["rmsd phi [°]";
                "rmsd theta [°]";
                "rmsd psi [°]";
                "max phi error [°]";
                "max theta error [°]";
                "max psi error [°]";
                "rmsd pos x [m]";
                "rmsd pos y [m]";
                "rmsd pos z [m]";
                "rmsd pos [m]";
                "rmsd vel x [m/s]";
                "rmsd vel y [m/s]";
                "rmsd vel z [m/s]";
                "rmsd vel [m/s]";];
    rmsd = [];
    rowNames = {};
    max_error = [];
    rmsd_pos = [];
    rmsd_vel = [];
    rowcount = 0;
    if filterSimulation.eulerEKF
        rmsd_euler = [RMSD_angle(truthdata.gamma_inertia(1,:), gamma_i_ahrsEuler(1,:)), ...
            RMSD_angle(truthdata.gamma_inertia(2,:), gamma_i_ahrsEuler(2,:)), ...
            RMSD_angle(truthdata.gamma_inertia(3,:), gamma_i_ahrsEuler(3,:))];
        rowNames{end+1} = 'ahrs_euler';
        rowcount = rowcount +1;
        rmsd(rowcount,:) = rmsd_euler;
        temp = abs(calculateAngleDifferenceDeg(truthdata.gamma_inertia,gamma_i_ahrsEuler));
        temp = max(temp(:,indexMaxError) , [], 2);
        max_error(rowcount, :) = temp;
        rmsd_pos(rowcount, :) = [-1, -1, -1, -1];
        rmsd_vel(rowcount, :) = [-1, -1, -1, -1];
    end
    if filterSimulation.eulerEKFHIComp
        rmsd_euler = [RMSD_angle(truthdata.gamma_inertia(1,:), gamma_i_ahrsEulerHIComp(1,:)), ...
            RMSD_angle(truthdata.gamma_inertia(2,:), gamma_i_ahrsEulerHIComp(2,:)), ...
            RMSD_angle(truthdata.gamma_inertia(3,:), gamma_i_ahrsEulerHIComp(3,:))];
        rowNames{end+1} = 'ahrs_eulerHIComp';
        rowcount = rowcount +1;
        rmsd(rowcount,:) = rmsd_euler;
        temp = abs(calculateAngleDifferenceDeg(truthdata.gamma_inertia,gamma_i_ahrsEulerHIComp));
        temp = max(temp(:,indexMaxError) , [], 2);
        max_error(rowcount, :) = temp;
        rmsd_pos(rowcount, :) = [-1, -1, -1, -1];
        rmsd_vel(rowcount, :) = [-1, -1, -1, -1];
    end
    if filterSimulation.eulerEKFWoutGyroBias
        rmsd_eulerWoutGyroBias = [RMSD_angle(truthdata.gamma_inertia(1,:), gamma_i_ahrsEulerWoutGyroBias(1,:)), ...
            RMSD_angle(truthdata.gamma_inertia(2,:), gamma_i_ahrsEulerWoutGyroBias(2,:)), ...
            RMSD_angle(truthdata.gamma_inertia(3,:), gamma_i_ahrsEulerWoutGyroBias(3,:))];
        rowNames{end+1} = 'ahrs_eulerWoutGyroBias';
        rowcount = rowcount +1;
        rmsd(rowcount,:) = rmsd_eulerWoutGyroBias;
        temp = abs(calculateAngleDifferenceDeg(truthdata.gamma_inertia,gamma_i_ahrsEulerWoutGyroBias));
        temp = max(temp(:,indexMaxError) , [], 2);
        max_error(rowcount, :) = temp;
        rmsd_pos(rowcount, :) = [-1, -1, -1, -1];
        rmsd_vel(rowcount, :) = [-1, -1, -1, -1];
    end
    if filterSimulation.eulerEKFDifferentSampleTime
        rmsd_eulerDST = [RMSD_angle(truthdata.gamma_inertia(1,:), gamma_i_ahrsEulerDifferentSampleTime(1,:)), ...
            RMSD_angle(truthdata.gamma_inertia(2,:), gamma_i_ahrsEulerDifferentSampleTime(2,:)), ...
            RMSD_angle(truthdata.gamma_inertia(3,:), gamma_i_ahrsEulerDifferentSampleTime(3,:))];
        rowNames{end+1} = 'ahrs_eulerDST';
        rowcount = rowcount +1;
        rmsd(rowcount,:) = rmsd_eulerDST;
        temp = abs(calculateAngleDifferenceDeg(truthdata.gamma_inertia,gamma_i_ahrsEulerDifferentSampleTime));
        temp = max(temp(:,indexMaxError) , [], 2);
        max_error(rowcount, :) = temp;
        rmsd_pos(rowcount, :) = [-1, -1, -1, -1];
        rmsd_vel(rowcount, :) = [-1, -1, -1, -1];
    end
    if filterSimulation.quaternionEKF
        rmsd_quaternionEKF = [RMSD_angle(truthdata.gamma_inertia(1,:), gamma_i_ahrsQuaternion(1,:)), ...
            RMSD_angle(truthdata.gamma_inertia(2,:), gamma_i_ahrsQuaternion(2,:)), ...
            RMSD_angle(truthdata.gamma_inertia(3,:), gamma_i_ahrsQuaternion(3,:))];
        rowNames{end+1} = 'ahrs_quaternion';
        rowcount = rowcount +1;
        rmsd(rowcount,:) = rmsd_quaternionEKF;
        temp = abs(calculateAngleDifferenceDeg(truthdata.gamma_inertia,gamma_i_ahrsQuaternion));
        temp = max(temp(:,indexMaxError) , [], 2);
        max_error(rowcount, :) = temp;
        rmsd_pos(rowcount, :) = [-1, -1, -1, -1];
        rmsd_vel(rowcount, :) = [-1, -1, -1, -1];
    end
    if filterSimulation.quaternionEKFDifferentSampleTime
        rmsd_quaternionEKFDST = [RMSD_angle(truthdata.gamma_inertia(1,:), gamma_i_ahrsQuaternionDifferentSampleTime(1,:)), ...
            RMSD_angle(truthdata.gamma_inertia(2,:), gamma_i_ahrsQuaternionDifferentSampleTime(2,:)), ...
            RMSD_angle(truthdata.gamma_inertia(3,:), gamma_i_ahrsQuaternionDifferentSampleTime(3,:))];
        rowNames{end+1} = 'ahrs_quaternionDST';
        rowcount = rowcount +1;
        rmsd(rowcount,:) = rmsd_quaternionEKFDST;
        temp = abs(calculateAngleDifferenceDeg(truthdata.gamma_inertia,gamma_i_ahrsQuaternionDifferentSampleTime));
        temp = max(temp(:,indexMaxError) , [], 2);
        max_error(rowcount, :) = temp;
        rmsd_pos(rowcount, :) = [-1, -1, -1, -1];
        rmsd_vel(rowcount, :) = [-1, -1, -1, -1];
    end
    if filterSimulation.MEKF
        rmsd_MEKF = [RMSD_angle(truthdata.gamma_inertia(1,:), gamma_i_ahrsMEKF(1,:)), ...
                               RMSD_angle(truthdata.gamma_inertia(2,:), gamma_i_ahrsMEKF(2,:)), ...
                               RMSD_angle(truthdata.gamma_inertia(3,:), gamma_i_ahrsMEKF(3,:))];
        rowNames{end+1} = 'ahrs_MEKF';
        rowcount = rowcount +1;
        rmsd(rowcount,:) = rmsd_MEKF;
        temp = abs(calculateAngleDifferenceDeg(truthdata.gamma_inertia,gamma_i_ahrsMEKF));
        temp = max(temp(:,indexMaxError) , [], 2);
        max_error(rowcount, :) = temp;
        rmsd_pos(rowcount, :) = [-1, -1, -1, -1];
        rmsd_vel(rowcount, :) = [-1, -1, -1, -1];
    end
    if filterSimulation.MEKFDifferentSampleTime
        rmsd_MEKFDST = [RMSD_angle(truthdata.gamma_inertia(1,:), gamma_i_ahrsMEKFDifferentSampleTime(1,:)), ...
                               RMSD_angle(truthdata.gamma_inertia(2,:), gamma_i_ahrsMEKFDifferentSampleTime(2,:)), ...
                               RMSD_angle(truthdata.gamma_inertia(3,:), gamma_i_ahrsMEKFDifferentSampleTime(3,:))];
        rowNames{end+1} = 'ahrs_MEKFDST';
        rowcount = rowcount +1;
        rmsd(rowcount,:) = rmsd_MEKFDST;
        temp = abs(calculateAngleDifferenceDeg(truthdata.gamma_inertia,gamma_i_ahrsMEKFDifferentSampleTime));
        temp = max(temp(:,indexMaxError) , [], 2);
        max_error(rowcount, :) = temp;
        rmsd_pos(rowcount, :) = [-1, -1, -1, -1];
        rmsd_vel(rowcount, :) = [-1, -1, -1, -1];
    end
    if filterSimulation.MEKFDifferentSampleTimeVarianceChange
        rmsd_MEKFDSTVarChng = [RMSD_angle(truthdata.gamma_inertia(1,:), gamma_i_ahrsMEKFDifferentSampleTimeVarianceChange(1,:)), ...
                               RMSD_angle(truthdata.gamma_inertia(2,:), gamma_i_ahrsMEKFDifferentSampleTimeVarianceChange(2,:)), ...
                               RMSD_angle(truthdata.gamma_inertia(3,:), gamma_i_ahrsMEKFDifferentSampleTimeVarianceChange(3,:))];
        rowNames{end+1} = 'ahrs_MEKFDSTVarChng';
        rowcount = rowcount +1;
        rmsd(rowcount,:) = rmsd_MEKFDSTVarChng;
        temp = abs(calculateAngleDifferenceDeg(truthdata.gamma_inertia,gamma_i_ahrsMEKFDifferentSampleTimeVarianceChange));
        temp = max(temp(:,indexMaxError) , [], 2);
        max_error(rowcount, :) = temp;
        rmsd_pos(rowcount, :) = [-1, -1, -1, -1];
        rmsd_vel(rowcount, :) = [-1, -1, -1, -1];
    end
    if filterSimulation.MEKFDifferentSampleTimeMagnetometer
        rmsd_MEKFDST_Mag = [RMSD_angle(truthdata.gamma_inertia(1,:), gamma_i_ahrsMEKFDifferentSampleTimeMagnetometer(1,:)), ...
                               RMSD_angle(truthdata.gamma_inertia(2,:), gamma_i_ahrsMEKFDifferentSampleTimeMagnetometer(2,:)), ...
                               RMSD_angle(truthdata.gamma_inertia(3,:), gamma_i_ahrsMEKFDifferentSampleTimeMagnetometer(3,:))];
        rowNames{end+1} = 'ahrs_MEKFDST_Mag';
        rowcount = rowcount +1;
        rmsd(rowcount,:) = rmsd_MEKFDST_Mag;
        temp = abs(calculateAngleDifferenceDeg(truthdata.gamma_inertia,gamma_i_ahrsMEKFDifferentSampleTimeMagnetometer));
        temp = max(temp(:,indexMaxError) , [], 2);
        max_error(rowcount, :) = temp;
        rmsd_pos(rowcount, :) = [-1, -1, -1, -1];
        rmsd_vel(rowcount, :) = [-1, -1, -1, -1];
    end
    if filterSimulation.MEKFDifferentSampleTimeMagnetometer2
        rmsd_MEKFDST_Mag2 = [RMSD_angle(truthdata.gamma_inertia(1,:), gamma_i_ahrsMEKFDifferentSampleTimeMagnetometer2(1,:)), ...
                               RMSD_angle(truthdata.gamma_inertia(2,:), gamma_i_ahrsMEKFDifferentSampleTimeMagnetometer2(2,:)), ...
                               RMSD_angle(truthdata.gamma_inertia(3,:), gamma_i_ahrsMEKFDifferentSampleTimeMagnetometer2(3,:))];
        rowNames{end+1} = 'ahrs_MEKFDST_Mag2';
        rowcount = rowcount +1;
        rmsd(rowcount,:) = rmsd_MEKFDST_Mag2;
        temp = abs(calculateAngleDifferenceDeg(truthdata.gamma_inertia,gamma_i_ahrsMEKFDifferentSampleTimeMagnetometer2));
        temp = max(temp(:,indexMaxError) , [], 2);
        max_error(rowcount, :) = temp;
        rmsd_pos(rowcount, :) = [-1, -1, -1, -1];
        rmsd_vel(rowcount, :) = [-1, -1, -1, -1];
    end
    if filterSimulation.MEKFwithoutCentripetalCompensation
        rmsd_MEKFwithoutCentripetalCompensation = [RMSD_angle(truthdata.gamma_inertia(1,:), gamma_i_ahrsMEKFwithoutCentripetalCompensation(1,:)), ...
                               RMSD_angle(truthdata.gamma_inertia(2,:), gamma_i_ahrsMEKFwithoutCentripetalCompensation(2,:)), ...
                               RMSD_angle(truthdata.gamma_inertia(3,:), gamma_i_ahrsMEKFwithoutCentripetalCompensation(3,:))];
        rowNames{end+1} = 'ahrs_MEKFw/outCC';
        rowcount = rowcount +1;
        rmsd(rowcount,:) = rmsd_MEKFwithoutCentripetalCompensation;
        temp = abs(calculateAngleDifferenceDeg(truthdata.gamma_inertia,gamma_i_ahrsMEKFwithoutCentripetalCompensation));
        temp = max(temp(:,indexMaxError) , [], 2);
        max_error(rowcount, :) = temp;
        rmsd_pos(rowcount, :) = [-1, -1, -1, -1];
        rmsd_vel(rowcount, :) = [-1, -1, -1, -1];
    end
    
    if filterSimulation.mahonyComplementaryFilter
        rmsd_mahonyEuler = [RMSD_angle(truthdata.gamma_inertia(1,:), gamma_i_mahonyComplementaryFilter(1,:)), ...
                               RMSD_angle(truthdata.gamma_inertia(2,:), gamma_i_mahonyComplementaryFilter(2,:)), ...
                               RMSD_angle(truthdata.gamma_inertia(3,:), gamma_i_mahonyComplementaryFilter(3,:))];
        rowNames{end+1} = 'mahonyComplementaryFilterEuler';
        rowcount = rowcount +1;
        rmsd(rowcount,:) = rmsd_mahonyEuler;
        temp = abs(calculateAngleDifferenceDeg(truthdata.gamma_inertia,gamma_i_mahonyComplementaryFilter));
        temp = max(temp(:,indexMaxError) , [], 2);
        max_error(rowcount, :) = temp;
        rmsd_pos(rowcount, :) = [-1, -1, -1, -1];
        rmsd_vel(rowcount, :) = [-1, -1, -1, -1];
    end
    
    if filterSimulation.mahonyExplicitComplementaryFilter
        rmsd_mahonyExplicitComplementaryFilter = [RMSD_angle(truthdata.gamma_inertia(1,:), gamma_i_mahonyExplicitComplementaryFilter(1,:)), ...
                               RMSD_angle(truthdata.gamma_inertia(2,:), gamma_i_mahonyExplicitComplementaryFilter(2,:)), ...
                               RMSD_angle(truthdata.gamma_inertia(3,:), gamma_i_mahonyExplicitComplementaryFilter(3,:))];
        rowNames{end+1} = ['mahony explicit' newline 'complementary filter'];
        rowcount = rowcount +1;
        rmsd(rowcount,:) = rmsd_mahonyExplicitComplementaryFilter;
        temp = abs(calculateAngleDifferenceDeg(truthdata.gamma_inertia,gamma_i_mahonyExplicitComplementaryFilter));
        temp = max(temp(:,indexMaxError) , [], 2);
        max_error(rowcount, :) = temp;
        rmsd_pos(rowcount, :) = [-1, -1, -1, -1];
        rmsd_vel(rowcount, :) = [-1, -1, -1, -1];
    end

    if filterSimulation.INSEulerEKF
        rmsd_gamma_i_insEulerEKF = [RMSD_angle(truthdata.gamma_inertia(1,:), gamma_i_insEulerEKF(1,:)), ...
                               RMSD_angle(truthdata.gamma_inertia(2,:), gamma_i_insEulerEKF(2,:)), ...
                               RMSD_angle(truthdata.gamma_inertia(3,:), gamma_i_insEulerEKF(3,:))];
        rowNames{end+1} = 'INSEulerEKF';
        rowcount = rowcount +1;
        rmsd(rowcount,:) = rmsd_gamma_i_insEulerEKF;
        temp = abs(calculateAngleDifferenceDeg(truthdata.gamma_inertia,gamma_i_insEulerEKF));
        temp = max(temp(:,indexMaxError) , [], 2);
        max_error(rowcount, :) = temp;
        rmsd_p = [RMSD(truthdata.p_i(1,:), p_insEulerEKF(1,:)), ...
                  RMSD(truthdata.p_i(2,:), p_insEulerEKF(2,:)), ...
                  RMSD(truthdata.p_i(3,:), p_insEulerEKF(3,:))];
        rmsd_p = [rmsd_p, sqrt(rmsd_p(1)^2+rmsd_p(2)^2 + rmsd_p(3)^2)];
        rmsd_v = [RMSD(truthdata.v_i(1,:), v_insEulerEKF(1,:)), ...
                  RMSD(truthdata.v_i(2,:), v_insEulerEKF(2,:)), ...
                  RMSD(truthdata.v_i(3,:), v_insEulerEKF(3,:))];
        rmsd_v = [rmsd_v, sqrt(rmsd_v(1)^2+rmsd_v(2)^2 + rmsd_v(3)^2)];
        rmsd_pos(rowcount, :) = rmsd_p;
        rmsd_vel(rowcount, :) = rmsd_v;
    end
    
    if filterSimulation.INSEulerEKFDST
        rmsd_gamma_i_insEulerEKFDST = [RMSD_angle(truthdata.gamma_inertia(1,:), gamma_i_insEulerEKFDST(1,:)), ...
                               RMSD_angle(truthdata.gamma_inertia(2,:), gamma_i_insEulerEKFDST(2,:)), ...
                               RMSD_angle(truthdata.gamma_inertia(3,:), gamma_i_insEulerEKFDST(3,:))];
        rowNames{end+1} = 'INSEulerEKFDST';
        rowcount = rowcount +1;
        rmsd(rowcount,:) = rmsd_gamma_i_insEulerEKFDST;
        temp = abs(calculateAngleDifferenceDeg(truthdata.gamma_inertia,gamma_i_insEulerEKFDST));
        temp = max(temp(:,indexMaxError) , [], 2);
        max_error(rowcount, :) = temp;
        rmsd_p = [RMSD(truthdata.p_i(1,:), p_insEulerEKFDST(1,:)), ...
                  RMSD(truthdata.p_i(2,:), p_insEulerEKFDST(2,:)), ...
                  RMSD(truthdata.p_i(3,:), p_insEulerEKFDST(3,:))];
        rmsd_p = [rmsd_p, sqrt(rmsd_p(1)^2+rmsd_p(2)^2 + rmsd_p(3)^2)];
        rmsd_v = [RMSD(truthdata.v_i(1,:), v_insEulerEKFDST(1,:)), ...
                  RMSD(truthdata.v_i(2,:), v_insEulerEKFDST(2,:)), ...
                  RMSD(truthdata.v_i(3,:), v_insEulerEKFDST(3,:))];
        rmsd_v = [rmsd_v, sqrt(rmsd_v(1)^2+rmsd_v(2)^2 + rmsd_v(3)^2)];
        rmsd_pos(rowcount, :) = rmsd_p;
        rmsd_vel(rowcount, :) = rmsd_v;
    end
    
    if filterSimulation.INSMEKF
        rmsd_gamma_i_insMEKF = [RMSD_angle(truthdata.gamma_inertia(1,:), gamma_i_insMEKF(1,:)), ...
                               RMSD_angle(truthdata.gamma_inertia(2,:), gamma_i_insMEKF(2,:)), ...
                               RMSD_angle(truthdata.gamma_inertia(3,:), gamma_i_insMEKF(3,:))];
        rowNames{end+1} = 'INSMEKF';
        rowcount = rowcount +1;
        rmsd(rowcount,:) = rmsd_gamma_i_insMEKF;
        temp = abs(calculateAngleDifferenceDeg(truthdata.gamma_inertia,gamma_i_insMEKF));
        temp = max(temp(:,indexMaxError) , [], 2);
        max_error(rowcount, :) = temp;
        rmsd_p = [RMSD(truthdata.p_i(1,:), p_insMEKF(1,:)), ...
                  RMSD(truthdata.p_i(2,:), p_insMEKF(2,:)), ...
                  RMSD(truthdata.p_i(3,:), p_insMEKF(3,:))];
        rmsd_p = [rmsd_p, sqrt(rmsd_p(1)^2+rmsd_p(2)^2 + rmsd_p(3)^2)];
        rmsd_v = [RMSD(truthdata.v_i(1,:), v_insMEKF(1,:)), ...
                  RMSD(truthdata.v_i(2,:), v_insMEKF(2,:)), ...
                  RMSD(truthdata.v_i(3,:), v_insMEKF(3,:))];
        rmsd_v = [rmsd_v, sqrt(rmsd_v(1)^2+rmsd_v(2)^2 + rmsd_v(3)^2)];
        rmsd_pos(rowcount, :) = rmsd_p;
        rmsd_vel(rowcount, :) = rmsd_v;
        
    end
    
    if filterSimulation.INSMEKFDST
        rmsd_gamma_i_INSMEKFDST = [RMSD_angle(truthdata.gamma_inertia(1,:), gamma_i_INSMEKFDST(1,:)), ...
                               RMSD_angle(truthdata.gamma_inertia(2,:), gamma_i_INSMEKFDST(2,:)), ...
                               RMSD_angle(truthdata.gamma_inertia(3,:), gamma_i_INSMEKFDST(3,:))];
        rowNames{end+1} = 'INSMEKFDST';
        rowcount = rowcount +1;
        rmsd(rowcount,:) = rmsd_gamma_i_INSMEKFDST;
        temp = abs(calculateAngleDifferenceDeg(truthdata.gamma_inertia,gamma_i_INSMEKFDST));
        temp = max(temp(:,indexMaxError) , [], 2);
        max_error(rowcount, :) = temp;
        rmsd_p = [RMSD(truthdata.p_i(1,:), p_INSMEKFDST(1,:)), ...
                  RMSD(truthdata.p_i(2,:), p_INSMEKFDST(2,:)), ...
                  RMSD(truthdata.p_i(3,:), p_INSMEKFDST(3,:))];
        rmsd_p = [rmsd_p, sqrt(rmsd_p(1)^2+rmsd_p(2)^2 + rmsd_p(3)^2)];
        rmsd_v = [RMSD(truthdata.v_i(1,:), v_INSMEKFDST(1,:)), ...
                  RMSD(truthdata.v_i(2,:), v_INSMEKFDST(2,:)), ...
                  RMSD(truthdata.v_i(3,:), v_INSMEKFDST(3,:))];
        rmsd_v = [rmsd_v, sqrt(rmsd_v(1)^2+rmsd_v(2)^2 + rmsd_v(3)^2)];
        rmsd_pos(rowcount, :) = rmsd_p;
        rmsd_vel(rowcount, :) = rmsd_v;
    end
    
    if filterSimulation.INSMEKFDSTVarianceChange
        rmsd_gamma_i_INSMEKFDSTVarianceChange = [RMSD_angle(truthdata.gamma_inertia(1,:), gamma_i_INSMEKFDSTVarianceChange(1,:)), ...
                               RMSD_angle(truthdata.gamma_inertia(2,:), gamma_i_INSMEKFDSTVarianceChange(2,:)), ...
                               RMSD_angle(truthdata.gamma_inertia(3,:), gamma_i_INSMEKFDSTVarianceChange(3,:))];
        rowNames{end+1} = 'INSMEKFDSTVarianceChange';
        rowcount = rowcount +1;
        rmsd(rowcount,:) = rmsd_gamma_i_INSMEKFDSTVarianceChange;
        temp = abs(calculateAngleDifferenceDeg(truthdata.gamma_inertia,gamma_i_INSMEKFDSTVarianceChange));
        temp = max(temp(:,indexMaxError) , [], 2);
        max_error(rowcount, :) = temp;
        rmsd_p = [RMSD(truthdata.p_i(1,:), p_INSMEKFDSTVarianceChange(1,:)), ...
                  RMSD(truthdata.p_i(2,:), p_INSMEKFDSTVarianceChange(2,:)), ...
                  RMSD(truthdata.p_i(3,:), p_INSMEKFDSTVarianceChange(3,:))];
        rmsd_p = [rmsd_p, sqrt(rmsd_p(1)^2+rmsd_p(2)^2 + rmsd_p(3)^2)];
        rmsd_v = [RMSD(truthdata.v_i(1,:), v_INSMEKFDSTVarianceChange(1,:)), ...
                  RMSD(truthdata.v_i(2,:), v_INSMEKFDSTVarianceChange(2,:)), ...
                  RMSD(truthdata.v_i(3,:), v_INSMEKFDSTVarianceChange(3,:))];
        rmsd_v = [rmsd_v, sqrt(rmsd_v(1)^2+rmsd_v(2)^2 + rmsd_v(3)^2)];
        rmsd_pos(rowcount, :) = rmsd_p;
        rmsd_vel(rowcount, :) = rmsd_v;
    end
    
    if filterSimulation.INSMEKFDST_tetherSag
        rmsd_gamma_i_INSMEKFDST_tetherSag = [RMSD_angle(truthdata.gamma_inertia(1,:), gamma_i_INSMEKFDST_tetherSag(1,:)), ...
                                               RMSD_angle(truthdata.gamma_inertia(2,:), gamma_i_INSMEKFDST_tetherSag(2,:)), ...
                                               RMSD_angle(truthdata.gamma_inertia(3,:), gamma_i_INSMEKFDST_tetherSag(3,:))];
        rowNames{end+1} = 'INSMEKFDST_tetherSag';
        rowcount = rowcount +1;
        rmsd(rowcount,:) = rmsd_gamma_i_INSMEKFDST_tetherSag;
        temp = abs(calculateAngleDifferenceDeg(truthdata.gamma_inertia,gamma_i_INSMEKFDST_tetherSag));
        temp = max(temp(:,indexMaxError) , [], 2);
        max_error(rowcount, :) = temp;
        rmsd_p = [RMSD(truthdata.p_i(1,:), p_INSMEKFDST_tetherSag(1,:)), ...
                  RMSD(truthdata.p_i(2,:), p_INSMEKFDST_tetherSag(2,:)), ...
                  RMSD(truthdata.p_i(3,:), p_INSMEKFDST_tetherSag(3,:))];
        rmsd_p = [rmsd_p, sqrt(rmsd_p(1)^2+rmsd_p(2)^2 + rmsd_p(3)^2)];
        rmsd_v = [RMSD(truthdata.v_i(1,:), v_INSMEKFDST_tetherSag(1,:)), ...
                  RMSD(truthdata.v_i(2,:), v_INSMEKFDST_tetherSag(2,:)), ...
                  RMSD(truthdata.v_i(3,:), v_INSMEKFDST_tetherSag(3,:))];
        rmsd_v = [rmsd_v, sqrt(rmsd_v(1)^2+rmsd_v(2)^2 + rmsd_v(3)^2)];
        rmsd_pos(rowcount, :) = rmsd_p;
        rmsd_vel(rowcount, :) = rmsd_v;
    end
    
    clear temp;
      
    if length(rmsd) > 0
         sizeTableError = size(rmsd);
        if sizeTableError(1) ~= length(rowNames)
            error('Number of rows from rmsd does not match number of rowNames');
        end
        sizeTableMaxError = size(max_error);
        if sizeTableError(2)+sizeTableMaxError(2) ~= length(colNames)
            %error('Number of columns from rmsd does not match number of colNames');
        end


        numRows = length(rowNames);
        numCols = length(colNames);
        errorTable = table([rmsd, max_error, rmsd_pos, rmsd_vel], 'RowNames', rowNames);
        figure;
        uitable('Data', errorTable{:,:}, 'ColumnName', colNames,'RowName', rowNames, 'units', 'normalized', 'position',[0, 0, 1, 1]);
        clear numRows numCols colNames rowCount rowNames rowcount sizeTableError sizeTableMaxError;
        
        tableExportPath = '../../Text/MatlabTableExport/';
        if createDifferentSampleTimeTable
            if simulationSampleRateForSensors
                warning('No latex table created, because simulationSampleRateForSensors is 1');
            else
                if disableGPSatToHighAccelerations
                   warning('No latex table created, because disableGPSatToHighAccelerations is 1'); 
                else
                    % create table
                    % angle
                    data = [rmsd];
                    columnNames = {'Filter & Name in \ref{fig:DST}', '\rmsd(\phi )','\rmsd(\theta )','\rmsd(\psi )'};
                    caption = 'Estimated angle error with different sensor availability handling';
                    rowNames = {'\makecell{\ac{INS} \ac{MEKF} \\ without interpolation} & INSMEKF', ...
                                '\makecell{	\ac{INS} \ac{MEKF} \\ no measurement step \\if no measurement available} & INSMEKFDST', ...
                                '\makecell{	\ac{INS} \ac{MEKF} \\ with interpolation} &\makecell{INSMEKF\\Interpolate}'};
                                
                    label = 'sampleRateSimulationAngle';
                    onlyData = 1;
                    latexTable.DSTAngle = createLatexTable(data, rowNames, columnNames, caption, label, onlyData);
                    fid = fopen([tableExportPath,label, datestr(datetime, 'yyyy-mm-dd_HH-MM-SS'), '.LatexTable'],'wt');
                    fprintf(fid, latexTable.DSTAngle);
                    fclose(fid);
                    % pos
                    data = [rmsd_pos];
                    columnNames = {'Filter & Name in \ref{fig:DST}', '\rmsd(x )','\rmsd(y)','\rmsd(z)', '\rmsd(||v||)'};
                    caption = 'Estimated position error with different sensor availability handling';      
                    label = 'sampleRateSimulationPos';
                    onlyData = 1;
                    latexTable.DSTPos = createLatexTable(data, rowNames, columnNames, caption, label, onlyData);
                    fid = fopen([tableExportPath,label, datestr(datetime, 'yyyy-mm-dd_HH-MM-SS'), '.LatexTable'],'wt');
                    fprintf(fid, latexTable.DSTPos);
                    fclose(fid);
                    
                    % vel
                    data = [rmsd_vel];
                    columnNames = {'Filter & Name in \ref{fig:DST}', '\rmsd(x)','\rmsd(y)','\rmsd(z)' , '\rmsd(||v||)'};
                    caption = 'Estimated position error with different sensor availability handling';
                    label = 'sampleRateSimulationVel';
                    onlyData = 1;
                    latexTable.DSTVel = createLatexTable(data, rowNames, columnNames, caption, label, onlyData);
                    fid = fopen([tableExportPath,label, datestr(datetime, 'yyyy-mm-dd_HH-MM-SS'), '.LatexTable'],'wt');
                    fprintf(fid, latexTable.DSTVel);
                    fclose(fid);
                end
            end
        end
        
        if create4gLimitTable
            if disableGPSatToHighAccelerations
                % create table
                % angle
                data = [rmsd];
                columnNames = {'Filter & Name in \ref{fig:DST}', '\rmsd(\phi )','\rmsd(\theta )','\rmsd(\psi )'};
                caption = 'Estimated angle error with \SI{4}{g} limitation of the GPS module';
                rowNames = {'\makecell{\ac{INS} \ac{MEKF} \\ without interpolation} & INSMEKF', ...
                            '\makecell{	\ac{INS} \ac{MEKF} \\ no measurement step \\if no measurement available} & INSMEKFDST', ...
                            '\makecell{	\ac{INS} \ac{MEKF} \\ with interpolation} &\makecell{INSMEKF\\Interpolate}'};

                label = 'angle4gLimit';
                onlyData = 1;
                latexTable.Limit4g.Angle = createLatexTable(data, rowNames, columnNames, caption, label, onlyData);
                fid = fopen([tableExportPath,label, datestr(datetime, 'yyyy-mm-dd_HH-MM-SS'), '.LatexTable'],'wt');
                fprintf(fid, latexTable.Limit4g.Angle);
                fclose(fid);
                
                % pos
                data = [rmsd_pos];
                columnNames = {'Filter & Name in \ref{fig:DST}', '\rmsd(x )','\rmsd(y)','\rmsd(z)', '\rmsd(||v||)'};
                caption = 'Estimated position error with \SI{4}{g} limitation of the GPS module';
                label = 'pos4gLimit';
                onlyData = 1;
                latexTable.Limit4g.Pos = createLatexTable(data, rowNames, columnNames, caption, label, onlyData);
                fid = fopen([tableExportPath,label, datestr(datetime, 'yyyy-mm-dd_HH-MM-SS'), '.LatexTable'],'wt');
                fprintf(fid, latexTable.Limit4g.Pos);
                fclose(fid);
                    
                % vel
                data = [rmsd_vel];
                columnNames = {'Filter & Name in \ref{fig:DST}', '\rmsd(x)','\rmsd(y)','\rmsd(z)' , '\rmsd(||v||)'};
                caption = 'Estimated velocity error with \SI{4}{g} limitation of the GPS module';
                label = 'vel4gLimit';
                onlyData = 1;
                latexTable.Limit4g.Vel = createLatexTable(data, rowNames, columnNames, caption, label, onlyData);
                fid = fopen([tableExportPath,label, datestr(datetime, 'yyyy-mm-dd_HH-MM-SS'), '.LatexTable'],'wt');
                fprintf(fid, latexTable.Limit4g.Vel);
                fclose(fid);
            else
                warning('No table created, because 4g limit is not turned on');
            end
        end
        
        if createTetherSagTable
            if ~tetherSag
                warning('No table created, because tetherSag=0');
            else
               % create table
                % angle
                data = [rmsd];
                columnNames = {'Filter & Name in \ref{fig:DST}', '\rmsd(\phi )','\rmsd(\theta )','\rmsd(\psi )'};
                caption = 'Estimated angle error with \SI{4}{g} limitation of the GPS module';
                rowNames = {'\makecell{\ac{INS} \ac{MEKF} \\ without tether sag compensation}', ...
                            '\makecell{	\ac{INS} \ac{MEKF} \\ with tether sag compensation}'};

                label = 'AngletetherSagComp';
                onlyData = 1;
                latexTable.tetherSag.Angle = createLatexTable(data, rowNames, columnNames, caption, label, onlyData);
                fid = fopen([tableExportPath,label, datestr(datetime, 'yyyy-mm-dd_HH-MM-SS'), '.LatexTable'],'wt');
                fprintf(fid, latexTable.tetherSag.Angle);
                fclose(fid);
                
                % pos
                data = [rmsd_pos];
                columnNames = {'Filter & Name in \ref{fig:DST}', '\rmsd(x )','\rmsd(y)','\rmsd(z)', '\rmsd(||v||)'};
                caption = 'Estimated position error with \SI{4}{g} limitation of the GPS module';
                label = 'PostetherSagComp';
                onlyData = 1;
                latexTable.tetherSag.Pos = createLatexTable(data, rowNames, columnNames, caption, label, onlyData);
                fid = fopen([tableExportPath,label, datestr(datetime, 'yyyy-mm-dd_HH-MM-SS'), '.LatexTable'],'wt');
                fprintf(fid, latexTable.tetherSag.Pos);
                fclose(fid);
                    
                % vel
                data = [rmsd_vel];
                columnNames = {'Filter & Name in \ref{fig:DST}', '\rmsd(x)','\rmsd(y)','\rmsd(z)' , '\rmsd(||v||)'};
                caption = 'Estimated velocity error with \SI{4}{g} limitation of the GPS module';
                label = 'VeltetherSagComp';
                onlyData = 1;
                latexTable.tetherSag.Vel = createLatexTable(data, rowNames, columnNames, caption, label, onlyData);
                fid = fopen([tableExportPath,label, datestr(datetime, 'yyyy-mm-dd_HH-MM-SS'), '.LatexTable'],'wt');
                fprintf(fid, latexTable.tetherSag.Vel);
                fclose(fid); 
            end
        end
        
    else
        warning('No simulations available to evaluate rmsd. No table created');
    end
    
    clear f numRows numCols skip indexMaxErrro lastIndex fid label caption columnNames rowNames onlyData;
end