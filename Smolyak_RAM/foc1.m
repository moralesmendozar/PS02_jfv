function [ff1] = foc1(k1T,k1I,kTp,bp,l1,l2,states,z1,z2, params,dKT1,dKI1,dl2)

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

    %prices
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
    
    
    
    
    kip = KIP(k1T,k1I,states,z1,z2,params,dl2);
    
    k1Tps = dKT1([kTp,kip,bp]);
    k1Ips = dKI1([kTp,kip,bp]);
    
    %to state z11
    p11 = mT1(z1i,1)*mT2(z1i,1);% to z11, 
    k1Tp = k1Tps(1);
    k1Ip = k1Ips(1);
    zprime = z1_grid(1);
    l1p = labor1(k1Tp,k1Ip,zprime,params);
    E11 = (1-theta1)*alpha1.*k1Ip.^(-1/lambda1)...
        .*(theta1*k1Tp^((lambda1-1)/lambda1) +...
        (1-theta1).*k1Ip^((lambda1-1)/lambda1) ).^(lambda1*alpha1/(lambda1-1)-1)...
        .*(exp(zprime)*l1p).^gamma1;
    
    %to state z12
    p12 = mT1(z1i,1)*mT2(z1i,2);% to z12
    k1Tp = k1Tps(2);
    k1Ip = k1Ips(2);
    %zprime = z1_grid(1);
    l1p = labor1(k1Tp,k1Ip,zprime,params);
    E12 = (1-theta1)*alpha1.*k1Ip.^(-1/lambda1)...
        .*(theta1*k1Tp^((lambda1-1)/lambda1) +...
        (1-theta1).*k1Ip^((lambda1-1)/lambda1) ).^(lambda1*alpha1/(lambda1-1)-1)...
        .*(exp(zprime)*l1p).^gamma1;
    
    %to state z21
    p21 = mT1(z1i,2)*mT2(z1i,1);% to z21
    k1Tp = k1Tps(3);
    k1Ip = k1Ips(3);
    zprime = z1_grid(2);
    l1p = labor1(k1Tp,k1Ip,zprime,params);
    E21 = (1-theta1)*alpha1.*k1Ip.^(-1/lambda1)...
        .*(theta1*k1Tp^((lambda1-1)/lambda1) +...
        (1-theta1).*k1Ip^((lambda1-1)/lambda1) ).^(lambda1*alpha1/(lambda1-1)-1)...
        .*(exp(zprime)*l1p).^gamma1;
    
    %to state z21
    p22 = mT1(z1i,2)*mT2(z1i,2);% to z22
    k1Tp = k1Tps(4);
    k1Ip = k1Ips(4);
    %zprime = z1_grid(2);
    l1p = labor1(k1Tp,k1Ip,zprime,params);
    E22 = (1-theta1)*alpha1.*k1Ip.^(-1/lambda1)...
        .*(theta1*k1Tp^((lambda1-1)/lambda1) +...
        (1-theta1).*k1Ip^((lambda1-1)/lambda1) ).^(lambda1*alpha1/(lambda1-1)-1)...
        .*(exp(zprime)*l1p).^gamma1;
    
    EE = p11.*E11 + p12.*E12 + p21.*E21 + p22.*E22;
    EE = EE./(1+rFirm);
    
    ff1 = theta1*alpha1*k1T.^(-1/lambda1)*...
        (theta1*k1T.^((lambda1-1)/lambda1)+...
        (1-theta1)*k1I.^((lambda1-1)/(lambda1)) ).^(lambda1*alpha1/(lambda1-1)-1) .*...
        (exp(z1).*l1).^gamma1 - EE*theta2*alpha2*(kT-k1T).^(-1/lambda2).*...
        (theta2*(kT -k1T).^((lambda2-1)/lambda2)+...
        (1-theta2)*(kI-k1I).^((lambda2-1)/(lambda2)) ).^(lambda2*alpha2/(lambda2-1)-1) .*...
        (exp(z2).*l2).^gamma2;
        
    
    %( wage/ ( EE.* (theta2*(kT-k1T).^((lambda2-1)/lambda2)+...
    %    (1-theta2)*(kI-k1I).^(lambda2/(lambda2-1)) ).^(lambda2*alpha2/(lambda2-1)) .*...
    %    gamma2*exp(z2*gamma2) ) ) .^( 1/(gamma2-1) );
    
end