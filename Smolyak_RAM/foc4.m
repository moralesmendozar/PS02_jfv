function [ff4] = foc4(k1T,k1I,kTp,bp,states,z1,z2, params,dl2,dB)

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
    
    bpps = dB([kTp,kip,bp]);
    
    %to state z11
    p11 = mT1(z1i,1)*mT2(z1i,1);% to z11, 
    bpp = bpps(1);
    E11 = - 0.04.* (bpp - bp);
    
    %to state z12
    p12 = mT1(z1i,1)*mT2(z1i,2);% to z12
    bpp = bpps(2);
    E12 = - 0.04.* (bpp - bp);
    
    %to state z21
    p21 = mT1(z1i,2)*mT2(z1i,1);% to z21
    bpp = bpps(3);
    E21 = - 0.04.* (bpp - bp);
    
    %to state z21
    p22 = mT1(z1i,2)*mT2(z1i,2);% to z22
    bpp = bpps(4);
    E22 = - 0.04.* (bpp - bp);
    
    EE = p11.*E11 + p12.*E12 + p21.*E21 + p22.*E22;
    EE = EE + 1+rBond +  0.02*(bp - 0.2) ;
    EE = EE./(1+rFirm);
    
    ff4 = 1 - 0.04*(bp-b) - EE;
    
end