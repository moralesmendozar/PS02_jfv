function [p] = getParams()
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
    
    % 1st technology
    params.z1_grid = [0.95; 1.05];
    params.mTransition1 =  [0.95 0.05; ...
                                0.05 0.95];

    % 2nd technology       
    params.z2_grid = [0.90; 1.10];
    params.mTransition2 =  [0.90 0.10; ...
                                0.10 0.90];
    
    p = params;
end