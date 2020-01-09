function [lab2] = labor2(k1T,k1I,kTp,bp,dKT11,dKT12,dKT21,dKT22,dKI11,dKI12,dKI21,dKI22,states,z1,z2,params)

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
    kip = KIP(k1T,k1I,states,z1,z2,params,dl2);
    k1Tp = dKT1(1,[kTp,kip,bp]);
    k1Ip = dKI1(1,[kTp,kip,bp]);
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
    
    
    
    % Expectation
    Evalue = 0;
    for iZ1next = 1:size(mTransition1,2)
         Evalue = Evalue + mTransition1(iZ1, iZ1next) * p.alpha1 * (1-p.theta1) .* (fhat_k1_I).^(-1/p.lambda1)*(p.theta1.*(fhat_k1_T).^((p.lambda1-1)/p.lambda1)* ...
          (1-p.theta1).*(fhat_k1_I).^((p.lambda1-1)/p.lambda1)).^((p.lambda1*p.alpha1/(p.lambda1-1))-1).* ...
            (exp(z1_grid(iZ1next)).*fhat_l1).^(p.gamma1);
    end
    
    % Result
    l2 = (1/(1+p.r_i)* p.gamma1 .* exp(z2).^p.gamma2 * (p.theta2.*(k_T - k1_T).^((p.lambda2-1)/p.lambda2) + ...
         (1-p.theta2).*(k_I - k1_I).^(lambda2/(lambda2-1))).^(lambda2*alpha2/(lambda2-1)) ./ w .* Evalue).^(1/(1-p.gamma2));
    
    
    lab2 = ( wage/ ( EE.* (theta2*(kT-k1T).^((lambda2-1)/lambda2)+...
        (1-theta2)*(kI-k1I).^(lambda2/(lambda2-1)) ).^(lambda2*alpha2/(lambda2-1)) .*...
        gamma2*exp(z2*gamma2) ) ) .^( 1/(gamma2-1) );
    
end