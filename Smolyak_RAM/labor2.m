function [lab2] = labor2(k1T,k1I,kTP,bp,dk1T,dk1I,states,params)

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
    z1 = states(4);
    z2 = states(5);
    
    if z1 == z1_grid(1)
        z1i =1;
    else
        z1i = 2;
    end
    if z2 == z2_grid(1)
        z2i =1;
    else
        z2i = 2;
    end
    
    kPI = ; 
    
    p11 = mT1(z1i,1)*mT2(z1i,1);% to z11, 
    k1Tp = dk1T(states);
    k1Ip = dk1I(states);
    E11 = (1-theta1)*alpha1.*k1Ip^(-1/lambda1)...
        .*(theta1*k1Tp^((lambda1-1)/lambda1) +...
        (1-theta1).*k1Ip^((lambda1-1)/lambda1) ).^(lambda1*alpha1/(lambda1-1)-1)...
        .*(e^(z1_grid(1))*labor1())
        
    p12 = mT1(z1i,1)*mT2(z1i,2);% to z12
    
    p21 = mT1(z1i,2)*mT2(z1i,1);% to z21
    p22 = mT1(z1i,2)*mT2(z1i,2);% to z22
    EE = p11.*E11 + p12.*E12 + p21.*E21 + p22.*E22;
    EE = EE./(1+rFirm);
    
    lab2 = ( wage/ ( EE.* (theta1*k1T.^((lambda1-1)/lambda1)+...
        (1-theta1)*k1I.^(lambda1/(lambda1-1)) ).^(lambda1*alpha1/(lambda1-1)) .*...
        gamma1*exp(z1*gamma1) ) ) .^( 1/(gamma1-1) );
    
end