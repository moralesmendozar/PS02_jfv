function [kip] = KIP(k1T,k1I,states,z1,z2,params,dl2)
    lambda1 = params.lambda1;
    theta1 = params.theta1;
    alpha1 = params.alpha1;
    gamma1 = params.gamma1;
    % 2nd technology
    lambda2 = params.lambda2;
    theta2 = params.theta2;
    alpha2 = params.alpha2;
    gamma2 = params.gamma2;

    % depreciation rates
    deltaT = params.deltaT;
    deltaI = params.deltaI;

    prices
    rFirm = params.rFirm;
    rBond = params.rBond;
    wage = params.wage;
    
    %technologies
    % 1st technology
    z1_grid = params.z1_grid;
    mTransition1 = params.mTransition1;
    mT1 = mTransition1;
    % 2nd technology       
    z2_grid = params.z2_grid;
    mTransition2 = params.mTransition2;
    mT2 = mTransition2;
    
    kT = states(1);
    kI = states(2);
    b = states(3);
%     z1 = states(4);
%     z2 = states(5);
    
    kip = (theta1.*(kT-k1T).^((lambda2-1)/lambda2) +...
        (1-theta2).*(kI-k1I).^((lambda2-1)/lambda2) ).^(lambda2*alpha2/(lambda2-1))...
        .*(exp(z2)*dl2(k1T,k1I,kpT,bp,states,params)).^gamma2+(1-deltaI).*kI;
end