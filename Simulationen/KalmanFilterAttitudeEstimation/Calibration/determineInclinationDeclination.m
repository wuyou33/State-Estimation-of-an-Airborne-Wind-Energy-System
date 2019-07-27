function [delta_d, delta_i] = DetermineInclinationDeclination(mag, acc)

% acc not used. Can be used to calculate from pitch and roll to delta_i and
% delta_d

    % Derivation see: Systemmatrix.wxm
    % Board must direct to geographic north

    mag_norm = mag./vecnorm(mag);
    mag_norm = mean(mag_norm')';
    %mag_norm = median(mag_norm')';

    mx = mag_norm(1);
    my = mag_norm(2);
    mz = mag_norm(3);

    delta_i = -(180*asin(mz))/pi;
    delta_d = (180*acos(mx/sqrt(1-mz^2)))/pi;

end