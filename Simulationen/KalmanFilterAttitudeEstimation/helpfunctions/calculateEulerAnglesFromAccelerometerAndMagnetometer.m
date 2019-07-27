% output angles are in degree!
% For the derivation see Systemmatrix.wxm
function [phi, theta, psi] = calculateEulerAnglesFromAccelerometerAndMagnetometer(acc, mag, delta_d, delta_i)

    a_n = vecnorm(acc);
    ax = acc(1,:)./a_n;
    ay = acc(2,:)./a_n;
    az = acc(3,:)./a_n;

    phi= atan2(ay,az);
    theta= asin(-ax/1);
    
    % use magnetometer only for yaw angle
    m_n = vecnorm(mag);
    m_x = mag(1,:)./m_n;
    m_y = mag(2,:)./m_n;
    m_z = mag(3,:)./m_n;
    
    a=(m_y.*sin(phi)+m_z.*cos(phi)).*sin(theta)+m_x.*cos(theta);
    b = m_y.*cos(phi)-m_z.*sin(phi);
    y = -((a.*sin((pi*delta_d)/180)+b.*cos((pi*delta_d)/180)));
    x = (a.*cos((pi*delta_d)/180)-b.*sin((pi*delta_d)/180));
    psi = atan2(y, x);

    phi = rad2deg(phi);
    theta = rad2deg(theta);
    psi = rad2deg(psi);
end

% test:
% delta_i = 0;
% delta_d = 0;
% gamma_test = [0;0;0];
% gamma_test = gamma_test*pi/180;
% a_b_meas = Re_to_b_euler(gamma_test)*[0;0;1];
% m_b_meas = Re_to_b_euler(gamma_test)*Rmag_to_e([1;0;0], delta_d, delta_i);
% [phi_log, theta_log, psi_log] = calculateEulerAnglesFromAccelerometerAndMagnetometer(a_b_meas, m_b_meas, delta_d, delta_i)