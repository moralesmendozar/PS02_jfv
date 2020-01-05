% -------------------------------------------------------------------------
% HOMEWORK 2 - ECON 714 
% Prof. Jesus Fernandez-Villaverde
% Due date: January 10th, 2020
% 
% Author: Juan Castellanos Silvan 
% -------------------------------------------------------------------------

clc;
close all;


%% 1. Economy description

% -----------------
% 1.1. Calibration
% -----------------
p = struct();

% 1st technology 
p.lambda1 = 1.1;
p.theta1 = 0.5;
p.alpha1 = 0.3;
p.gamma1 = 0.6;

% 2nd technology
p.lambda1 = 0.9;
p.theta1 = 0.4;
p.alpha1 = 0.4;
p.gamma1 = 0.5;

% depreciation rates
p.deltaT = 0.10;
p.deltaI = 0.15;

% prices
p.r_firm = 0.025;
p.r = 0.03;
p.w = 1;

% --------------------------
% 1.2. Stochastic processes
% --------------------------

% 1st technology
z1_grid = [0.95; 1.05];
mTransition1 = [0.95 0.05; ...
                0.05 0.95];

% 2nd technology       
z2_grid = [0.90; 1.10];
mTransition2 = [0.90 0.10; ...
                0.10 0.90];

%--------------------------
% 1.3. Declare size vectors
%--------------------------

% Smolyak's algorithm
n = 3;                             % Number of state variables 
mu = 2;                            % Precision of the approximation
q = n + mu;                        % Indexes the size of the grid
nGrid = 1 + 4*n + 4*(n*(n-1)/2);   % Number of points

% Euler errors
grid_num  = 3000;                  % # of grid points for capital (to compute euler errors)

% Simulation
T         = 10000;                 % Number of periods for the simulation of the economy
dropT     = 1000;                  % Burn-in

%% Steady State (need to compute)


%% Smolyak's algorithm
%{
% Transform the domian of the state varaibles
x_min = min(x_grid);
x_max = max(x_grid);    
x_scaled_down = 2 * (x_grid - x_min) ./ (x_max - x_min) - 1;
%}

% Building a Sparse Grid
a = [-1,-1,-1];
b = [1,1,1];

[x, s] = smolyak_grid(n,mu,a,b);

figure(1)
scatter3(x(:,1),x(:,2),x(:,3), 'filled')


