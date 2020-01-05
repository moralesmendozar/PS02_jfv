% Test script for pseudospectral methods:
% orthogonal points versus uniform grid
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

theta0 = [0.0 0.1 0.05 0.025];
a= 0;
b= 6;

%----------------------------------------------------------------
% 2. Solving
%----------------------------------------------------------------

x = col_points(a,b,4);
options=optimset('Display','off');
[theta,~] = fsolve('nonlinear',theta0,options,x,a,b);

%----------------------------------------------------------------
% 3. Results
%----------------------------------------------------------------

theta = real(theta);
y = cheby_approx(0,theta,a,b);
theta = theta./y;										% Boundary condition
y = cheby_approx(0:0.1:6,theta,a,b);
y1 = 0:0.1:6;
y2 = exp(-y1);
plot(y1,y2,'--',y1,y)
hold on

x =[0.01 2 4 5.99];
options=optimset('Display','off');
[theta,fval] = fsolve('nonlinear',theta,options,x,a,b);

theta = real(theta);
y = cheby_approx(0,theta,a,b);
theta = theta./y;
y = cheby_approx(0:0.1:6,theta,a,b);
plot(0:0.1:6,y,'-r')

toc