function [lab2] = labor2(k1T,k1I,dKT11,dKT12,dKT21,dKT22,dKI11,dKI12,dKI21,dKI22,states,z1,z2,params)

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
    
    %to state z11
    p11 = mT1(z1i,1)*mT2(z1i,1);% to z11, 
    k1Tp = dKT11(states);
    k1Ip = dKI11(states);
    zprime = z1_grid(1);
    l1p = labor1(k1Tp,k1Ip,zprime,params);
    E11 = (1-theta1)*alpha1.*k1Ip.^(-1/lambda1)...
        .*(theta1*k1Tp^((lambda1-1)/lambda1) +...
        (1-theta1).*k1Ip^((lambda1-1)/lambda1) ).^(lambda1*alpha1/(lambda1-1)-1)...
        .*(e^(zprime)*l1p).^gamma1;
    
    %to state z12
    p12 = mT1(z1i,1)*mT2(z1i,2);% to z12
    k1Tp = dKT12(states);
    k1Ip = dKI12(states);
    zprime = z1_grid(1);
    l1p = labor1(k1Tp,k1Ip,zprime,params);
    E12 = (1-theta1)*alpha1.*k1Ip.^(-1/lambda1)...
        .*(theta1*k1Tp^((lambda1-1)/lambda1) +...
        (1-theta1).*k1Ip^((lambda1-1)/lambda1) ).^(lambda1*alpha1/(lambda1-1)-1)...
        .*(e^(zprime)*l1p).^gamma1;
    
    %to state z21
    p21 = mT1(z1i,2)*mT2(z1i,1);% to z21
    k1Tp = dKT21(states);
    k1Ip = dKI21(states);
    zprime = z1_grid(2);
    l1p = labor1(k1Tp,k1Ip,zprime,params);
    E21 = (1-theta1)*alpha1.*k1Ip.^(-1/lambda1)...
        .*(theta1*k1Tp^((lambda1-1)/lambda1) +...
        (1-theta1).*k1Ip^((lambda1-1)/lambda1) ).^(lambda1*alpha1/(lambda1-1)-1)...
        .*(e^(zprime)*l1p).^gamma1;
    
    %to state z21
    p22 = mT1(z1i,2)*mT2(z1i,2);% to z22
    k1Tp = dKT22(states);
    k1Ip = dKI22(states);
    zprime = z1_grid(2);
    l1p = labor1(k1Tp,k1Ip,zprime,params);
    E22 = (1-theta1)*alpha1.*k1Ip.^(-1/lambda1)...
        .*(theta1*k1Tp^((lambda1-1)/lambda1) +...
        (1-theta1).*k1Ip^((lambda1-1)/lambda1) ).^(lambda1*alpha1/(lambda1-1)-1)...
        .*(e^(zprime)*l1p).^gamma1;
    
    EE = p11.*E11 + p12.*E12 + p21.*E21 + p22.*E22;
    EE = EE./(1+rFirm);
    
    lab2 = ( wage/ ( EE.* (theta1*k1T.^((lambda1-1)/lambda1)+...
        (1-theta1)*k1I.^(lambda1/(lambda1-1)) ).^(lambda1*alpha1/(lambda1-1)) .*...
        gamma1*exp(z2*gamma1) ) ) .^( 1/(gamma1-1) );
    
end