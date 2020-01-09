function [focs l1 l2] = FOCS(d,states,z1,z2, params,dKT11,dKT12,dKT21,dKT22,dKI11,dKI12,dKI21,dKI22)
    
    k1T = d(1);
    k1I = d(2);
    kTp = d(3);
    bp = d(4);
    l1 = labor1(k1T,k1I,z1,params);
    l2 = labor2(k1T,k1I,kTp,bp,dKT11,dKT12,dKT21,dKT22,dKI11,dKI12,dKI21,dKI22,states,z1,z2,params dl2);
    
    f1 = foc1(k1T,k1I,kTp,bp,l1,l2,states,z1,z2, params);
    f2 = foc2(k1T,k1I,kTp,bp,l1,l2,states,z1,z2, params);
    f3 = foc3(k1T,k1I,kTp,bp,l1,l2,states,z1,z2, params);
    f4 = foc4(k1T,k1I,kTp,bp,l1,l2,states,z1,z2, params);
    
    focs = [f1; f2; f3; f4];
    
end