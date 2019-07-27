% finding a circle in a 3D dataset
% source: [1] https://meshlogic.github.io/posts/jupyter/curve-fitting/fitting-a-circle-to-cluster-of-3d-points/

function [center, radius, normalvector] = findCircleIn3DLeastSquare(data)

center = mean(data')';

% not implemented
radius = 0;
normalvector = zeros(3, 0);

end