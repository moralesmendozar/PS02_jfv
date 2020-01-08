% -------------------------------------------------------------------------
% HOMEWORK 2 - ECON 714 
% Prof. Jesus Fernandez-Villaverde
% Due date: January 10th, 2020
% 
% Author: Shasha Wang
% Reference: Rodrigo, Juan, Sara's code
% Discussion with: Rodrigo A. Mendoza, Juan
% -------------------------------------------------------------------------

clear;
% close all;

%% 1. Economy description

% ----------------------------------------------------------
% 1.1. Calibration
% ----------------------------------------------------------
% 1st technology 
lambda1 = 1.1;
theta1 = 0.5;
alpha1 = 0.3;
gamma1 = 0.6;

% 2nd technology
lambda2 = 0.9;
theta2 = 0.4;
alpha2 = 0.4;
gamma2 = 0.5;

% depreciation rates
deltaT = 0.10;
deltaI = 0.15;

% prices
rFirm = 0.025;
rBond = 0.03;
wage = 1;

% save parameters in an array
% can be passed to functions later

params = struct();

% 1st technology 
params.lambda1 = lambda1;
params.theta1 = theta1;
params.alpha1 = alpha1;
params.gamma1 = gamma1;

% 2nd technology
params.lambda2 = lambda2;
params.theta2 = theta2;
params.alpha2 = alpha2;
params.gamma2 = gamma2;

% depreciation rates
params.deltaT = deltaT;
params.deltaI = deltaI;

% prices
params.rFirm = rFirm;
params.rBond = rBond;
params.wage = wage;

% ----------------------------------------------------------
% 1.2. Stochastic processes
% ----------------------------------------------------------

% 1st technology
z1_grid = exp([0.95; 1.05]);
mTransition1 =  [0.95 0.05; ...
                            0.05 0.95];

% 2nd technology       
z2_grid = exp([0.90; 1.10]);
mTransition2 =  [0.90 0.10; ...
                            0.10 0.90];

%% Steady State (need to compute)
  
% input_ss = [kT,kT_1,kI,kI_1,labor_1,labor_2]
input_ss_initial=[0.9,0.2,0.5,0.4,0.5,0.5];
opts1 = optimoptions('fsolve','Tolx',1e-9, 'Display','iter');

display('Start Solving the Function');
SS = fsolve(@(input_ss) steadyStateFunction(input_ss,params), input_ss_initial,opts1);
display('Solution Obtained.');

fprintf('The Candidate Solution Is Found to Be: %2.4f \n', SS);
fprintf('The Function Value At this Candidate Solution Is: %2.6f \n', ...
        steadyStateFunction(SS,params));
    
kT_ss = SS(1);
kT_1_ss = SS(2);
kI_ss = SS(3);
kI_1_ss = SS(4);
labor_1_ss = SS(5);
labor_2_ss = SS(6);
bond_ss = (rFirm - rBond)/0.02 + 0.2;

T = table(kT_ss,kT_1_ss,kI_ss,kI_1_ss,labor_1_ss,labor_2_ss,bond_ss)

%       kT_ss        kT_1_ss       kI_ss       kI_1_ss     labor_1_ss    labor_2_ss    bond_ss
%     __________    _________    _________    _________    __________    __________    _______
% 
%     0.00023636    3.622e-08    0.0032974    5.245e-05    7.7954e-05    0.00064264     -0.05 

% 
% 
% % ----------------------------------------------------------
% % 1.3. Declare size vectors
% % ----------------------------------------------------------
% 
% % Smolyak's algorithm
% n = 3;                             % Number of state variables 
% mu = 2;                            % Precision of the approximation
% q = n + mu;                        % Indexes the size of the grid
% nGrid = 1 + 4*n + 4*(n*(n-1)/2);   % Number of points
% 
% % Euler errors
% grid_num  = 3000;                  % # of grid points for capital (to compute euler errors)
% 
% % Simulation
% T         = 10000;                 % Number of periods for the simulation of the economy
% dropT     = 1000;                  % Burn-in
% 
% 
% 
% %% 3. PROJECTION: The Smolyak's algorithm
% 
% %{
% % Transform the domian of the state varaibles
% x_min = min(x_grid);
% x_max = max(x_grid);    
% x_scaled_down = 2 * (x_grid - x_min) ./ (x_max - x_min) - 1;
% %}
% 
% % -------------------------
% % 3.1. Building a Sparse Grid
% % -------------------------
% a = [-1,-1,-1];
% b = [1,1,1];
% 
% [x, s] = smolyak_grid(n,mu,a,b);
% 
% % Graph
% pGrid = figure(1);
% scatter3(x(:,1),x(:,2),x(:,3), 'filled');
% 
% X = 29.7;                  % A4 paper size
% Y = 21.0;                  % A4 paper size
% 
% xMargin = 2;               % left/right margins from page borders
% yMargin = 2;               % bottom/top margins from page borders
% xSize = X - 2*xMargin;     % figure size on paper (widht & hieght)
% ySize = Y - 2*yMargin;     % figure size on paper (widht & hieght)
% 
% set(gcf, 'Units','centimeters', 'Position',[xSize ySize 0 0]/2)
% set(gcf, 'PaperUnits','centimeters')
% set(gcf, 'PaperSize',[X Y])
% set(gcf, 'PaperPosition',[xMargin yMargin xSize ySize])
% set(gcf, 'PaperOrientation','portrait')
% 
% %saveas(pGrid,'/Users/Castesil/Documents/EUI/Year II - PENN/Computational Economics/Homework 2/LaTeX/pSmolyakGrid', 'pdf')
% 
% % ----------------------------------------
% % 3.2. Finding the approximating function
% % ----------------------------------------
% 
% 
% 
% %% 4. PERTURBATION: a third order approx.
% 
% % ------------------
% % 4.1. Housekeeping
% % ------------------
% % clc;
% % close all;
% %addpath '/Applications/dynare/4.5.7';
% 
% % -----------------
% % 4.2. Call dynare
% % ----------------
% 
% %dynare firm_optm

