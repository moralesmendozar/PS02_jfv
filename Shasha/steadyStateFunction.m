function [y] = steady_state_function(input_ss,params) % input is kT,kT_1,kI,kI_1,labor_1,labor_2
    % 1st technology 
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

    % prices
    rFirm = params.rFirm;
    rBond = params.rBond;
    wage = params.wage;

%  input_ss(1)=kT, (2)kT_1, (3)kI, (4)kI_1, (5)labor_1, (6)labor_2
% k1i
    y(1) = 1 + rFirm - (1/(1+rFirm))...
        *((1-theta1)*alpha1*input_ss(4)^(-1/lambda1)*(theta1*input_ss(2)^((lambda1-1)/lambda1) +...
        (1-theta1)*input_ss(4)^((lambda1-1)/lambda1))^(lambda1*alpha1/(lambda1-1)-1)*(exp(1)*input_ss(5))^gamma1)...
        *((theta2)*alpha2*(input_ss(1)-input_ss(2))^(-1/lambda2)*(theta2*(input_ss(1)-input_ss(2))^((lambda2-1)/lambda2) + (1-theta2)*(input_ss(3)-input_ss(4))^((lambda2-1)/lambda2))^(lambda2*alpha2/(lambda2-1)-1)*(exp(1)*input_ss(6))^gamma2);
%     y(2) = (1-theta1)*alpha1* input_ss(2)^(-1/lambda1) * (exp(1)*input_ss(5))^gamma1 * (theta1*input_ss(2)^((lambda1-1)/lambda1) + (1-theta1)*input_ss(4)^((lambda1-1)/lambda1))^(lambda1*alpha1/(lambda1-1)-1)...
%         -(1/(1+rFirm)) * ((1-theta1)*alpha1*)...
%         * ();
    y(2) = 1 - (1/(1+rFirm)) * ( 1 - deltaI + (exp(1) * input_ss(6))^gamma2 * (1-theta2) * alpha2 * (input_ss(3)-input_ss(4))^(-1/lambda2) * (theta2*(input_ss(1)-input_ss(2))^((lambda2-1)/lambda2) + (1-theta2)*(input_ss(3)-input_ss(4))^((lambda2-1)/lambda2))^(lambda2*alpha2/(lambda2-1)-1));
 
    y(3) = 1 - (1/(1+rFirm)) * (1 - deltaT + (exp(1)*input_ss(5))^gamma1 * theta1 * alpha1 *input_ss(1)^(-1/lambda1)*(theta1*input_ss(2)^((lambda1-1)/lambda1) + (1-theta1)*input_ss(4)^((lambda1-1)/lambda1))^(lambda1*alpha1/(lambda1-1)-1));
    y(4) = wage - gamma1*input_ss(5)^(gamma1-1) * exp(gamma1) * (theta1*input_ss(2)^((lambda1-1)/lambda1) + (1-theta1)*input_ss(4)^((lambda1-1)/lambda1))^(lambda1*alpha1/(lambda1-1));
    y(5) = wage - (1/(1+rFirm)) * gamma2 * input_ss(6)^(gamma2-1)*exp(gamma2)* (theta2*(input_ss(1)-input_ss(2))^((lambda2-1)/lambda2) + (1-theta2)*(input_ss(3)-input_ss(4))^((lambda2-1)/lambda2))^(lambda2*alpha2/(lambda2-1))...
        * ((1-theta1)*alpha1*input_ss(4)^(-1/lambda1)*(theta1*input_ss(2)^((lambda1-1)/lambda1) + (1-theta1)*input_ss(4)^((lambda1-1)/lambda1))^(lambda1*alpha1/(lambda1-1)-1)*(exp(1)*input_ss(5))^gamma1);
    y(6) = deltaI*input_ss(3) - (exp(1)*input_ss(6))^gamma2 * (theta2*(input_ss(1)-input_ss(2))^((lambda2-1)/lambda2) + (1-theta2)*(input_ss(3)-input_ss(4))^((lambda2-1)/lambda2))^(lambda2*alpha2/(lambda2-1));
end