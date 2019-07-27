% Root Mean Square Deviation
% https://en.wikipedia.org/wiki/Root-mean-square_deviation
% used to estimate the quality of the filter
% y: real value
% y_hat: estimated value
function rmsd = RMSD(y, y_hat)
    if length(y) ~= length(y_hat)
        error('RMSD: y and y_hat must have the same length!')
    end
    N = length(y);
    e = y - y_hat;
    e_square = e.^2;
    rmsd = sqrt(sum(e_square)/N);
end