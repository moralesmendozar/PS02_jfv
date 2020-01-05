function y = nonlinear(theta,x,a,b)

% nonlinear system of equations

	y = residual(x(1:length(theta)),theta,a,b);