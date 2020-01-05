% Test script for pseudospectral methods:
% 4 grade versus 6 expansion
%
% Test case:
%
% d'(x)+d(x)=0
%
% with d(0)=1, x=(0,6)
%
% Jesus Fernandez-Villaverde, August 3, 2006

%----------------------------------------------------------------
% 0. Housekeeping
%----------------------------------------------------------------

clear variables
close all
tic

%----------------------------------------------------------------
% 1. Numerical Parameters
%----------------------------------------------------------------

theta0 = [0.0 0.1 0.05 0.05];
a= 0;
b= 6;

%----------------------------------------------------------------
% 2. Solving 4
%----------------------------------------------------------------

x = col_points(a,b,4);
options=optimset('Display','off');
[theta,~] = fsolve('nonlinear',theta0,options,x,a,b);

%----------------------------------------------------------------
% 3. Results
%----------------------------------------------------------------

theta = real(theta);
y = cheby_approx(0,theta,a,b);
theta = theta./y;
y = cheby_approx(0:0.1:6,theta,a,b);
y1 = 0:0.1:6;
y2 = exp(-y1);
plot(y1,y2,'--',y1,y)
hold on

%----------------------------------------------------------------
% 4. Solving 6
%----------------------------------------------------------------

x = col_points(a,b,6);
theta0 = [theta, 0.0325   -0.0094];
options = optimset('Display','off');
[theta,fval] = fsolve('nonlinear',theta0,options,x,a,b);
fval;

%----------------------------------------------------------------
% 5. Results
%----------------------------------------------------------------

theta = real(theta);
y = cheby_approx(0,theta,a,b);
theta = theta./y;
y = cheby_approx(0:0.1:6,theta,a,b);
plot(0:0.1:6,y,'-g')

toc