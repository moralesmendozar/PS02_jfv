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
p.lambda2 = 0.9;
p.theta2 = 0.4;
p.alpha2 = 0.4;
p.gamma2 = 0.5;

% depreciation rates
p.deltaT = 0.10;
p.deltaI = 0.15;

% prices
p.r_i = 0.025;
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

%% 2. STEADY STATE

% Solving for [k1_T, k1_I, k_T, K_I, l1, l2] using FOC + LoM
x0 =[0.00590718, 0.00191017, 0.00598928,  0.00195713, 0.00393336, 9.32542e-05];
options = optimoptions('fsolve','Display','iter');
[x,fval] = fsolve(@(x)ss_foc(x,p),x0,options);

k1_T_ss = x(1);
k1_I_ss = x(2);
k_T_ss = x(3);
k_I_ss = x(4);
l1_ss = x(5);
l2_ss = x(6);

% Solving for b using FOC w.r.t. b'
b_ss = (p.r_i - p.r)/0.02 + 0.2;

results = table(k1_T_ss,k1_I_ss,k_T_ss,k_I_ss,l1_ss,l2_ss,b_ss);
display(results)


%{
                  FSOLVE W/ DYNARE'S VALUES AS INITIAL GUESS
k1_T_ss      k1_I_ss      k_T_ss      k_I_ss        l1_ss        l2_ss       b_ss 
_________    _________    ________    _________    _________    __________    _____
0.0098306    0.0016105    0.010018    0.0016721    0.0045426    0.00010418    -0.05
                  FSOLVE W/ SHASHA'S VALUES AS INITIAL GUESS
k1_T_ss      k1_I_ss      k_T_ss       k_I_ss        l1_ss        l2_ss       b_ss 
_________    _________    _________    _________    _________    __________    _____
0.0086168    0.0012747    0.0090318    0.0014001    0.0039714    0.00023721    -0.05
            DYNARE
value            		 0.115007
k1_T             		 0.00590718
k1_I             		 0.00191017
l1               		 0.00393336
k2_T             		 8.21054e-05
k2_I             		 4.69562e-05
l2               		 9.32542e-05
b                		 -0.05
y1               		 0.0065556
y2               		 0.000195713
s                		 0.00595667
x                		 0.000598928
d                		 0.00280506
k_T              		 0.00598928
k_I              		 0.00195713
z1               		 0
z2               		 0
%}

%% 3. PROJECTION: The Smolyak's algorithm

% ------------------
% 3.1. Housekeeping
% ------------------
d = 3;                                    % Number of state variables 
mu = 2;                                   % Precision of the approximation
q = d + mu;                               % Indexes the size of the grid
nGridSmolyak = 1 + 4*d + 4*(d*(d-1)/2);   % Number of points for Smolyak grid
nIntpPoints = 50;                         % Number of points for Smolyak interp.

a = 2*[-1,-1,-1];
b = 4*[1,1,1];

% ---------------------------
% 3.2. Building a Sparse Grid
% ----------------------------
[x, s] = smolyak_grid(d,mu,a,b);

% Convert the x points to [-1,1]^d via linear transformation
x_scaled_down = NaN(size(x));
for dind = 1:d
    x_scaled_down(:,dind) = 2.d0 * (x(:,dind)- s.a(dind)) / (s.b(dind)-s.a(dind)) - 1.d0;
end

% Graph
pGrid = figure(1);
scatter3(x_scaled_down(:,1),x_scaled_down(:,2),x_scaled_down(:,3), 'filled');

X = 29.7;                  % A4 paper size
Y = 21.0;                  % A4 paper size

xMargin = 2;               % left/right margins from page borders
yMargin = 2;               % bottom/top margins from page borders
xSize = X - 2*xMargin;     % figure size on paper (widht & hieght)
ySize = Y - 2*yMargin;     % figure size on paper (widht & hieght)

set(gcf, 'Units','centimeters', 'Position',[xSize ySize 0 0]/2)
set(gcf, 'PaperUnits','centimeters')
set(gcf, 'PaperSize',[X Y])
set(gcf, 'PaperPosition',[xMargin yMargin xSize ySize])
set(gcf, 'PaperOrientation','portrait')

saveas(pGrid,'/Users/Castesil/Documents/EUI/Year II - PENN/Computational Economics/Homework 2/LaTeX/pSmolyakGrid', 'pdf')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEED TO FIND k1_T k1_I kT' b l1 l2 evaluated at x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% state variables
exo_states = cartprod(z1_grid, z2_grid);
states = cartprod(exo_states, x);

nGrid = size(states,1);
nExoStates = size(exo_states,1);


% ------------------
% Capital functions 
% ------------------
% Pre-allocation
fx_k1_T = ones(nGrid,1);
fx_k1_I = ones(nGrid,1);
fx_k_T = ones(nGrid,1);


% Smolyak => fhat_k1_T, fhat_k1_I, fhat_k_T evaluated at tomorrow's states


% -----------------------
% optimal bond holdings
% -----------------------
b_vec = states(:,3);
%fx_b = optimal_bond();
%fx_b = reshape(fx_b, nGridSmolyak, nExoStates);


% ----------------
% optimal labor 1
% ----------------
% Find f(x)
z1_vec = states(:,1);
fx_l1 = optimal_labor1(fx_k1_T, fx_k1_I, z1_vec, p);
fx_l1 = reshape(fx_l1, nGridSmolyak, nExoStates);

% Smolyak coefficients 
[s_l1] = smolyak_coefficients(fx_l1,s);

% I should evaluate the function at the values of the states tomorrow...
% Evalutate approximated f at whatever values one wishes
xpts = 2*rand(nIntpPoints,d)-1; % get random draw from [-1,1]^d
for dind = 1:d
    xpts(:,dind) = (xpts(:,dind)+1)*(b(dind)-a(dind))/2 + a(dind); % convert to space defined by a,b
end

fhat_l1 = smolyak_fnc(s_l1,xpts);


% ----------------
% optimal labor 2
% ----------------
% Find f(x)
fx_l2 = optimal_labor2(fx_k1_T, fx_k1_I, fhat_k1_T, fhat_k1_I, fhat_l1, states, z1_grid, p);
fx_l2 = reshape(fx_l2, nGridSmolyak, nExoStates);

% Smolyak coefficients 
[s_l2] = smolyak_coefficients(fx_l2,s);

% I should evaluate the function at the values of the states tomorrow...
% Evalutate approximated f at whatever values one wishes
xpts = 2*rand(nIntpPoints,d)-1; % get random draw from [-1,1]^d
for dind = 1:d
    xpts(:,dind) = (xpts(:,dind)+1)*(b(dind)-a(dind))/2 + a(dind); % convert to space defined by a,b
end

fhat_l2 = smolyak_fnc(s_l2,xpts);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example for any given function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcn = @(x) sum(x.^3,2); 

% ----------------------------------------
% 3.2. Finding the approximating function
% ----------------------------------------
% 3.2.1. Smolyak coefficients 
fx = fcn(x);
[s] = smolyak_coefficients(fx,s);

% Evalutate approximated f at whatever values one wishes
xpts = 2*rand(32,d)-1; % get random draw from [-1,1]^d
 for dind = 1:d
      xpts(:,dind) = (xpts(:,dind)+1)*(b(dind)-a(dind))/2 + a(dind); % convert to space defined by a,b
 end

% 3.2.2. Smolyak polynomial
fhat = smolyak_fnc(s,xpts);

% Graph - compare true value with residual
figure(2)
plot(fhat - fnc(xpts))

figure(3)
plot3(xpts(:,1),xpts(:,2),fhat,'x',xpts(:,1),xpts(:,2),fcn(xpts),'o') 
grid on;

return
%% 4. PERTURBATION: a third order approx.

% ------------------
% 4.1. Housekeeping
% ------------------
clc;
close all;

% -----------------
% 4.2. Call dynare
% ----------------

dynare firm_optm