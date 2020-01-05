% This file solves the stochastic neoclassical growth model using
% spectral methods. In particular I use orthogonal collocation with the
% value function approximated by N Chebychev polynomials, with the highest
% polynomial of order N.

% This file calls the following functions:
% (i) tauchen.m
% (ii)residual_fcn.m: solves for the coefficients associated to the
% Chebychev polynomials
% (iii)eulerr_grid.m: computes Euler errors on a grid for capital and the
% exogenous shock.
% (iv)simulation.m

%----------------------------------------------------------------
% 0 Housekeeping
%----------------------------------------------------------------

clc
clear all
close all
format short

%----------------------------------------------------------------
% 1 Calibration
%----------------------------------------------------------------

% Calibration quarterly frequency
alpha= 0.3;                       % Capital Share
beta= 0.991 ;                     % Discount factor
delta= 0.0196;                    % Depreciation

% Preference parameters
gamma=5;
psi=1.5;
theta=(1-gamma)*psi/(psi-1);

% Technology process
lambda=0.95;      % Persistence parameter of the productivity shock
sigma=0.007;     % S.D. of the productivity shock Z

% Declare size vectors
shock_num=17;     % number of nodes for technology process Z
grid_num=3000;    % # of grid points for  capital (to compute euler errors)
listSize=2;  
nodeSizeList = [3,11]; % number of polynomials
T=10000;          % Number of periods for the simulation of the economy
dropT=1000;

%----------------------------------------------------------------
% 1 Steady State + Tauchen 
%----------------------------------------------------------------

% Compute Steady State values

lss = 1/3;
kss = ((1-beta*(1-delta))/(alpha*beta))^(1/(alpha-1))*lss;
css = kss^alpha*lss^(1-alpha)-delta*kss;
nu = css/((1-alpha)*kss^alpha*lss^(-alpha)*(1-lss)+css);
vss = css^nu*(1-lss)^(1-nu);

[Z,PI]=tauchen(shock_num,0,lambda,sigma,3);

%----------------------------------------------------------------
% 2. Spectral Method using Chebychev Polynomials
%----------------------------------------------------------------

% Define boundaries for capital
if(sigma == 0.007)
    cover_grid = 0.25;
    k_min = kss*(1 - cover_grid);
    k_max = kss*(1 + cover_grid);
    interval = kss*2* cover_grid;
elseif( sigma == 0.035)
    k_min = 3;
    k_max = 26.5;
    interval = k_max - k_min;
else
    disp('sigma should be 0.007 or 0.035')
end

tic
for multinode_step=1:listSize
    
    node_num = nodeSizeList(multinode_step);
    M = node_num*shock_num;
    
    % Find Zeros of the Chebychev Polynomial on order M 
    ZC=-cos((2*(1:node_num)'-1)*pi/(2*node_num));

    % Define Chebychev polynomials
    T_k=ones(node_num,node_num);
    T_k(:,2)=ZC;

    for i1=3:node_num
        T_k(:,i1)=2*ZC.*T_k(:,i1-1)-T_k(:,i1-2);
    end

    % Project collocation points in the K space
    grid_k=((ZC+1)*(k_max-k_min))/2+k_min;

    % Initial Guess for Chebyshev coefficients
    rho_guess = zeros(2*M,1);
    if(multinode_step == 1)
        for z_index = 1:shock_num
            rho_guess((z_index-1)*node_num+1) = vss;
            rho_guess((z_index-1)*node_num+M+1) = lss;
        end
    else
        for z_index = 1:2*shock_num
            rho_guess((z_index-1)*node_num+1: (z_index-1)*node_num+node_num_old)...
            = rho_old((z_index-1)*node_num_old+1: z_index*node_num_old);
        end
    end
    
    % Solve for Chebyshev coefficients
    rho = residual_fcn(alpha,beta,delta,nu,gamma,theta,k_min,k_max,rho_guess,grid_k,T_k,Z,PI,node_num,shock_num,M);
    rho_old = rho;
    node_num_old = node_num;
end
toc

%----------------------------------------
% Compute Euler Errors and Decision rules
%----------------------------------------

grid_k_complete = zeros(grid_num,1);
for i = 1:grid_num
    grid_k_complete(i) = k_min + (i-1)*interval/(grid_num-1);
end

[g_k,g_c,g_l,value_fcn,euler_error,max_error]= ...
                    eulerr_grid(alpha,beta,delta,gamma,theta,nu,rho,Z,PI,...
                    k_min,k_max,grid_k_complete,shock_num,node_num,grid_num,M);


[vSeries,kSeries,cSeries,lSeries,ySeries,...
    eeSeries,rbSeries,rkSeries,rkCondSeries] = ...
    simulation(alpha,beta,delta,gamma,theta,nu,kss,rho,Z,PI,k_min,k_max,node_num,shock_num,M,T,dropT);

mean_error = sum(eeSeries)/(T-dropT);
max_error_sim = max(eeSeries);

disp('Integral of Euler Equation Error:')
disp(mean_error)
disp('Max Euler Equation Error Simulation:')
disp(max_error_sim)

% -------
% Figures
% -------

% Decision Rules
figure(1)
subplot(221)
plot(grid_k_complete,value_fcn)
title('Value Function')
subplot(222)
plot(grid_k_complete,g_c)
title('Consumption Decision Rule')
subplot(223)
plot(grid_k_complete,g_l)
title('Labor Decision Rule')
subplot(224)
plot(grid_k_complete,g_k,grid_k_complete,grid_k_complete)
title('Capital Decision Rule')

% Euler Equation Error on the Grid
figure(2)
plot(grid_k_complete,euler_error)
title('Log10 Euler Error')

% Distribution of simulated variables

[f_c,x_c] = ksdensity(cSeries);
[f_k,x_k] = ksdensity(kSeries);
[f_y,x_y] = ksdensity(ySeries);
[f_l,x_l] = ksdensity(lSeries);
[f_rb,x_rb] = ksdensity(rbSeries);
[f_rk,x_rk] = ksdensity(rkSeries);

figure(3)
subplot(3,2,1);
    plot(x_c,f_c)
    title('Density of Consumption')
subplot(3,2,2);
    plot(x_l,f_l)
    title('Density of Consumption')
subplot(3,2,3);
    plot(x_k,f_k)
    title('Density of Capital')
subplot(3,2,4);
    plot(x_y,f_y)
    title('Density of Output')
subplot(3,2,5);
    plot(x_rb,f_rb)
    title('Density of Return on Risk Free bond')
subplot(3,2,6);
    plot(x_rk,f_rk)
    title('Density of Return on Equity')