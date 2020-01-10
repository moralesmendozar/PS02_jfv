function [F] = ss_foc(x,p) 
    
    % ---------------------------------------------------------------------
    % Note:
    % x(1) = k1_T, x(2) = k1_I, x(3) = kT, x(4)= kI, x(5) = l1, x(6) = l2
    % --------------------------------------------------------------------- 
    r_i = p.rFirm;
    w = p.wage;
    % 1. FOC w.r.t. k1_T
    F(1) =  p.theta1/(1-p.theta1)*(x(1)/x(2))^(-1/p.lambda1) - ...
        (1/(1+r_i))*p.theta2*p.alpha2*(x(3)-x(1))^(-1/p.lambda2)*(p.theta2*(x(3)-x(1))^((p.lambda2 - 1)/p.lambda2) +...
        (1-p.theta2)*(x(4)-x(2))^((p.lambda2 - 1)/p.lambda2))^(p.lambda2*p.alpha2/(p.lambda2-1)-1)* x(6)^p.gamma2;
           
    
    % 2. FOC w.r.t. k1_I
    F(2) = 1 - (1/(1+r_i))*(1-p.theta2)*p.alpha2*(x(4)-...
        x(2))^(-1/p.lambda2)*(p.theta2*(x(3)-x(1))^((p.lambda2 - 1)/p.lambda2) +...
        (1-p.theta2)*(x(4)-x(2))^((p.lambda2 - 1)/p.lambda2))^(p.lambda2*p.alpha2/(p.lambda2-1)-1)* x(6)^p.gamma2;
        
    
    % 3. FOC w.r.t. kT' 
    F(3) = 1 - 1/(1+r_i) * (p.alpha1 * p.theta1 * x(1)^(-1/p.lambda1) * (p.theta1*x(1)^((p.lambda1-1)/p.lambda1) + (1-p.theta1)*x(2)^((p.lambda1-1)/p.lambda1))^(p.lambda1*p.alpha1/(p.lambda1-1)-1)* x(5)^p.gamma1 + (1-p.deltaT));
    
    
    % 4. LoM of intangible capital
    F(4) = p.deltaI*x(4) - (p.theta2*(x(3)-x(1))^((p.lambda2 - 1)/p.lambda2) + ...
        (1-p.theta2)*(x(4)-x(2))^((p.lambda2 - 1)/p.lambda2))^(p.lambda2*p.alpha2/(p.lambda2-1)) * x(6)^p.gamma2;
    
    
    % 5. FOC w.r.t. l1
    F(5) = w - p.gamma1*x(5)^(p.gamma1-1)*(p.theta1*x(1)^((p.lambda1-1)/p.lambda1) + (1-p.theta1)*x(2)^((p.lambda1-1)/p.lambda1))^(p.lambda1*p.alpha1/(p.lambda1-1));
    
    
    % 6. FOC w.r.t. l2
    F(6) = w - (1/(1 + r_i)) * p.gamma2 * x(6)^(p.gamma2-1)*(p.theta2*(x(3)-x(1))^((p.lambda2-1)/p.lambda2) + (1-p.theta2)*(x(4)-x(2))^((p.lambda2-1)/p.lambda2))^(p.lambda2*p.alpha2/(p.lambda2-1))...
        * p.alpha1*(1-p.theta1)*x(2)^(-1/p.lambda1)*(p.theta1*x(1)^((p.lambda1-1)/p.lambda1) + (1-p.theta1)*x(2)^((p.lambda1-1)/p.lambda1))^(p.lambda1*p.alpha1/(p.lambda1-1)-1)* x(5)^p.gamma1;

   
end