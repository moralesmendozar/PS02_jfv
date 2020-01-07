function [focs] = FOCS(d,states, params)
    
    k1T = d(1);
    k1I = d(2);
    kpT = d(3);
    bp = d(4);
    
    l1 = ;
    l2 = ;
    
    f1 = foc1(k1T,k1I,kpT,bp,l1,l2,states, params);
    f2 = foc2(k1T,k1I,kpT,bp,l1,l2,states, params);
    f3 = foc3(k1T,k1I,kpT,bp,l1,l2,states, params);
    f4 = foc4(k1T,k1I,kpT,bp,l1,l2,states, params);
    
    focs = [f1; f2; f3; f4];
    
end