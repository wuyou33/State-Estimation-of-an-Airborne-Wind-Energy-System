function sigma_2 = estimateVariance(sensorData)
plotting = 0;
acc = sensorData.meas.a_b'*9.81; % [g] --> [m/s^2]
mag = sensorData.meas.m_b'; % [1]
gyr = sensorData.meas.w_b'; % [rad/s]

if plotting
    figure();
    title('PDF of acc_x');
    subplot(3,1,1);
    histogram(acc(:,1),  'Normalization','pdf'); %plot estimated pdf from the generated data
    subplot(3,1,2);
    histogram(acc(:,2),  'Normalization','pdf'); %plot estimated pdf from the generated data
    subplot(3,1,3);
    histogram(acc(:,3),  'Normalization','pdf'); %plot estimated pdf from the generated data
end

mean = sum(acc)/length(acc);
variance = (sqrt(sum((acc-mean).^2))/length(acc))';
sigma_2.a_b = variance.^2;

mean = sum(mag)/length(mag);
variance = (sqrt(sum((mag-mean).^2))/length(mag))';
sigma_2.m_b = variance.^2;

mean = sum(gyr)/length(gyr);
variance = (sqrt(sum((gyr-mean).^2))/length(gyr))';
sigma_2.w_b = variance.^2;

end