function [focs] = FOCS(d,l1,l2,states, params,dk1T,dk1I)
    
    k1T = d(1);
    k1I = d(2);
    kTp = d(3);
    bp = d(4);
    l1 = labor1(k1T,k1I,kTp,bp,states,params);
    l2 = labor2(k1T,k1I,kTp,bp,dk1T,dk1I,states,params);
    f1 = foc1(k1T,k1I,kTp,bp,l1,l2,states, params);
    f2 = foc2(k1T,k1I,kTp,bp,l1,l2,states, params);
    f3 = foc3(k1T,k1I,kTp,bp,l1,l2,states, params);
    f4 = foc4(k1T,k1I,kTp,bp,l1,l2,states, params);
    
    focs = [f1; f2; f3; f4];
    
end