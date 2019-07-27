% Calculates the variance for the line Angle estimation
% it's quite difficult to find a algebraic function for the variance, due
% to the nonlinearity and the size of the equations
% for different angles the variance varies to due to the nonlinearity, but
% for this configuration I found the highest variance
function sigma2Pos = tetherEstimatePositionVariance(sigma2_lineAngle, sigma2_tetherLength)

    nbrElements = 1000;
    plotting = 0;
    
    theta_base = 1 * pi/180 * ones(1, nbrElements);
    theta_kite = 70 * pi/180 * ones(1, nbrElements);
    phi_sag = 0 * pi/180 * ones(1, nbrElements);
    l = 150 * ones(1, nbrElements);

    theta_base = theta_base + sqrt(sigma2_lineAngle)*randn(1, nbrElements);
    theta_kite = theta_kite + sqrt(sigma2_lineAngle)*randn(1, nbrElements);
    phi_sag = phi_sag + sqrt(sigma2_lineAngle)*randn(1, nbrElements);
    length_sag = l + sqrt(sigma2_tetherLength) * randn(1, nbrElements);

    pos = zeros(3, nbrElements);
    
    %TODO: vectorize this function!
    for i = 1: nbrElements
        [x, y, z, ~, ~] = calculatePositionFromLineAngle(theta_base(i), theta_kite(i), phi_sag(i), length_sag(i));
        pos(:, i) = [x; y; z];
    end
    
    mean = sum(pos, 2)./length(pos);
    variance = sqrt(sum((pos-mean).^2, 2)/length(pos));
    sigma2Pos = max(variance.^2)*ones(3,1);
    
    
    if plotting
        figure();
        title('PDF of position');
        subplot(3,1,1);
        histogram(pos(1, :),  'Normalization','pdf'); %plot estimated pdf from the generated data
        subplot(3,1,2);
        histogram(pos(2, :),  'Normalization','pdf'); %plot estimated pdf from the generated data
        subplot(3,1,3);
        histogram(pos(3, :),  'Normalization','pdf'); %plot estimated pdf from the generated data
    end
    
end