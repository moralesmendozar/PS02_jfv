function y = residual(x,theta,a,b)

% Residual function

	syms v;
	y = diff(cheby_approx(v,theta,a,b),'v');
	y = subs(y,v,x);
	y = double(y);
	y = y+cheby_approx(x,theta,a,b);
	