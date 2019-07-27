function R = Re_to_b_euler(gamma)
    phi = gamma(1, :);
    theta = gamma(2, :);
    psi = gamma(3, :);
    R = [cos(psi).*cos(theta),	sin(psi).*cos(theta),	-sin(theta);
		sin(phi).*cos(psi).*sin(theta)-cos(phi).*sin(psi),	sin(phi).*sin(psi).*sin(theta)+cos(phi).*cos(psi),	sin(phi).*cos(theta);
		cos(phi).*cos(psi).*sin(theta)+sin(phi).*sin(psi),	cos(phi).*sin(psi).*sin(theta)-sin(phi).*cos(psi),	cos(phi).*cos(theta)];
end