function [residual, g1, g2, g3] = rbc_log_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [nperiods by M_.exo_nbr] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   steady_state  [M_.endo_nbr by 1] double       vector of steady state values
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations.
%                                          Dynare may prepend auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(8, 1);
T30 = params(2)*y(6)^params(3);
T44 = y(1)^params(4);
T48 = (y(6)*exp(y(10)))^(1-params(4));
lhs =1;
rhs =params(1)*y(4)/y(11)*(1+y(12)-params(5));
residual(1)= lhs-rhs;
lhs =y(8);
rhs =params(4)*y(3)/y(1);
residual(2)= lhs-rhs;
lhs =y(4)*T30;
rhs =y(9);
residual(3)= lhs-rhs;
lhs =y(9);
rhs =y(3)*(1-params(4))/y(6);
residual(4)= lhs-rhs;
lhs =y(7);
rhs =y(1)*(1-params(5))+y(5);
residual(5)= lhs-rhs;
lhs =y(3);
rhs =T44*T48;
residual(6)= lhs-rhs;
lhs =y(4)+y(5);
rhs =y(3);
residual(7)= lhs-rhs;
lhs =y(10);
rhs =params(6)*y(2)+params(7)*x(it_, 1);
residual(8)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(8, 13);

  %
  % Jacobian matrix
  %

T76 = params(2)*getPowerDeriv(y(6),params(3),1);
T82 = getPowerDeriv(y(6)*exp(y(10)),1-params(4),1);
T91 = getPowerDeriv(y(1),params(4),1);
  g1(1,4)=(-((1+y(12)-params(5))*params(1)*1/y(11)));
  g1(1,11)=(-((1+y(12)-params(5))*params(1)*(-y(4))/(y(11)*y(11))));
  g1(1,12)=(-(params(1)*y(4)/y(11)));
  g1(2,3)=(-(params(4)/y(1)));
  g1(2,1)=(-((-(params(4)*y(3)))/(y(1)*y(1))));
  g1(2,8)=1;
  g1(3,4)=T30;
  g1(3,6)=y(4)*T76;
  g1(3,9)=(-1);
  g1(4,3)=(-((1-params(4))/y(6)));
  g1(4,6)=(-((-(y(3)*(1-params(4))))/(y(6)*y(6))));
  g1(4,9)=1;
  g1(5,5)=(-1);
  g1(5,1)=(-(1-params(5)));
  g1(5,7)=1;
  g1(6,3)=1;
  g1(6,6)=(-(T44*exp(y(10))*T82));
  g1(6,1)=(-(T48*T91));
  g1(6,10)=(-(T44*y(6)*exp(y(10))*T82));
  g1(7,3)=(-1);
  g1(7,4)=1;
  g1(7,5)=1;
  g1(8,2)=(-params(6));
  g1(8,10)=1;
  g1(8,13)=(-params(7));

if nargout >= 3,
  %
  % Hessian matrix
  %

  v2 = zeros(25,3);
T124 = params(2)*getPowerDeriv(y(6),params(3),2);
T135 = getPowerDeriv(y(6)*exp(y(10)),1-params(4),2);
T136 = exp(y(10))*T135;
T142 = getPowerDeriv(y(1),params(4),2);
  v2(1,1)=1;
  v2(1,2)=134;
  v2(1,3)=(-((1+y(12)-params(5))*params(1)*(-1)/(y(11)*y(11))));
  v2(2,1)=1;
  v2(2,2)=50;
  v2(2,3)=  v2(1,3);
  v2(3,1)=1;
  v2(3,2)=141;
  v2(3,3)=(-((1+y(12)-params(5))*params(1)*(-((-y(4))*(y(11)+y(11))))/(y(11)*y(11)*y(11)*y(11))));
  v2(4,1)=1;
  v2(4,2)=147;
  v2(4,3)=(-(params(1)*1/y(11)));
  v2(5,1)=1;
  v2(5,2)=51;
  v2(5,3)=  v2(4,3);
  v2(6,1)=1;
  v2(6,2)=154;
  v2(6,3)=(-(params(1)*(-y(4))/(y(11)*y(11))));
  v2(7,1)=1;
  v2(7,2)=142;
  v2(7,3)=  v2(6,3);
  v2(8,1)=2;
  v2(8,2)=3;
  v2(8,3)=(-((-params(4))/(y(1)*y(1))));
  v2(9,1)=2;
  v2(9,2)=27;
  v2(9,3)=  v2(8,3);
  v2(10,1)=2;
  v2(10,2)=1;
  v2(10,3)=(-((-((-(params(4)*y(3)))*(y(1)+y(1))))/(y(1)*y(1)*y(1)*y(1))));
  v2(11,1)=3;
  v2(11,2)=69;
  v2(11,3)=T76;
  v2(12,1)=3;
  v2(12,2)=45;
  v2(12,3)=  v2(11,3);
  v2(13,1)=3;
  v2(13,2)=71;
  v2(13,3)=y(4)*T124;
  v2(14,1)=4;
  v2(14,2)=68;
  v2(14,3)=(-((-(1-params(4)))/(y(6)*y(6))));
  v2(15,1)=4;
  v2(15,2)=32;
  v2(15,3)=  v2(14,3);
  v2(16,1)=4;
  v2(16,2)=71;
  v2(16,3)=(-((-((-(y(3)*(1-params(4))))*(y(6)+y(6))))/(y(6)*y(6)*y(6)*y(6))));
  v2(17,1)=6;
  v2(17,2)=71;
  v2(17,3)=(-(T44*exp(y(10))*T136));
  v2(18,1)=6;
  v2(18,2)=6;
  v2(18,3)=(-(exp(y(10))*T82*T91));
  v2(19,1)=6;
  v2(19,2)=66;
  v2(19,3)=  v2(18,3);
  v2(20,1)=6;
  v2(20,2)=1;
  v2(20,3)=(-(T48*T142));
  v2(21,1)=6;
  v2(21,2)=123;
  v2(21,3)=(-(T44*(exp(y(10))*T82+y(6)*exp(y(10))*T136)));
  v2(22,1)=6;
  v2(22,2)=75;
  v2(22,3)=  v2(21,3);
  v2(23,1)=6;
  v2(23,2)=118;
  v2(23,3)=(-(T91*y(6)*exp(y(10))*T82));
  v2(24,1)=6;
  v2(24,2)=10;
  v2(24,3)=  v2(23,3);
  v2(25,1)=6;
  v2(25,2)=127;
  v2(25,3)=(-(T44*(y(6)*exp(y(10))*T82+y(6)*exp(y(10))*y(6)*exp(y(10))*T135)));
  g2 = sparse(v2(:,1),v2(:,2),v2(:,3),8,169);
if nargout >= 4,
  %
  % Third order derivatives
  %

  v3 = zeros(52,3);
T206 = getPowerDeriv(y(6)*exp(y(10)),1-params(4),3);
T207 = exp(y(10))*T206;
T208 = exp(y(10))*T207;
  v3(1,1)=1;
  v3(1,2)=1824;
  v3(1,3)=(-((1+y(12)-params(5))*params(1)*(y(11)+y(11))/(y(11)*y(11)*y(11)*y(11))));
  v3(2,1)=1;
  v3(2,2)=648;
  v3(2,3)=  v3(1,3);
  v3(3,1)=1;
  v3(3,2)=1740;
  v3(3,3)=  v3(1,3);
  v3(4,1)=1;
  v3(4,2)=1831;
  v3(4,3)=(-((1+y(12)-params(5))*params(1)*(y(11)*y(11)*y(11)*y(11)*(-(2*(-y(4))))-(-((-y(4))*(y(11)+y(11))))*(y(11)*y(11)*(y(11)+y(11))+y(11)*y(11)*(y(11)+y(11))))/(y(11)*y(11)*y(11)*y(11)*y(11)*y(11)*y(11)*y(11))));
  v3(5,1)=1;
  v3(5,2)=1993;
  v3(5,3)=(-(params(1)*(-1)/(y(11)*y(11))));
  v3(6,1)=1;
  v3(6,2)=649;
  v3(6,3)=  v3(5,3);
  v3(7,1)=1;
  v3(7,2)=661;
  v3(7,3)=  v3(5,3);
  v3(8,1)=1;
  v3(8,2)=1741;
  v3(8,3)=  v3(5,3);
  v3(9,1)=1;
  v3(9,2)=1837;
  v3(9,3)=  v3(5,3);
  v3(10,1)=1;
  v3(10,2)=1909;
  v3(10,3)=  v3(5,3);
  v3(11,1)=1;
  v3(11,2)=2000;
  v3(11,3)=(-(params(1)*(-((-y(4))*(y(11)+y(11))))/(y(11)*y(11)*y(11)*y(11))));
  v3(12,1)=1;
  v3(12,2)=1832;
  v3(12,3)=  v3(11,3);
  v3(13,1)=1;
  v3(13,2)=1844;
  v3(13,3)=  v3(11,3);
  v3(14,1)=2;
  v3(14,2)=3;
  v3(14,3)=(-((-((-params(4))*(y(1)+y(1))))/(y(1)*y(1)*y(1)*y(1))));
  v3(15,1)=2;
  v3(15,2)=27;
  v3(15,3)=  v3(14,3);
  v3(16,1)=2;
  v3(16,2)=339;
  v3(16,3)=  v3(14,3);
  v3(17,1)=2;
  v3(17,2)=1;
  v3(17,3)=(-((y(1)*y(1)*y(1)*y(1)*(-(2*(-(params(4)*y(3)))))-(-((-(params(4)*y(3)))*(y(1)+y(1))))*(y(1)*y(1)*(y(1)+y(1))+y(1)*y(1)*(y(1)+y(1))))/(y(1)*y(1)*y(1)*y(1)*y(1)*y(1)*y(1)*y(1))));
  v3(18,1)=3;
  v3(18,2)=914;
  v3(18,3)=T124;
  v3(19,1)=3;
  v3(19,2)=578;
  v3(19,3)=  v3(18,3);
  v3(20,1)=3;
  v3(20,2)=890;
  v3(20,3)=  v3(18,3);
  v3(21,1)=3;
  v3(21,2)=916;
  v3(21,3)=y(4)*params(2)*getPowerDeriv(y(6),params(3),3);
  v3(22,1)=4;
  v3(22,2)=913;
  v3(22,3)=(-((-((-(1-params(4)))*(y(6)+y(6))))/(y(6)*y(6)*y(6)*y(6))));
  v3(23,1)=4;
  v3(23,2)=409;
  v3(23,3)=  v3(22,3);
  v3(24,1)=4;
  v3(24,2)=877;
  v3(24,3)=  v3(22,3);
  v3(25,1)=4;
  v3(25,2)=916;
  v3(25,3)=(-((y(6)*y(6)*y(6)*y(6)*(-(2*(-(y(3)*(1-params(4))))))-(-((-(y(3)*(1-params(4))))*(y(6)+y(6))))*(y(6)*y(6)*(y(6)+y(6))+y(6)*y(6)*(y(6)+y(6))))/(y(6)*y(6)*y(6)*y(6)*y(6)*y(6)*y(6)*y(6))));
  v3(26,1)=6;
  v3(26,2)=916;
  v3(26,3)=(-(T44*exp(y(10))*T208));
  v3(27,1)=6;
  v3(27,2)=71;
  v3(27,3)=(-(T91*exp(y(10))*T136));
  v3(28,1)=6;
  v3(28,2)=851;
  v3(28,3)=  v3(27,3);
  v3(29,1)=6;
  v3(29,2)=911;
  v3(29,3)=  v3(27,3);
  v3(30,1)=6;
  v3(30,2)=6;
  v3(30,3)=(-(exp(y(10))*T82*T142));
  v3(31,1)=6;
  v3(31,2)=66;
  v3(31,3)=  v3(30,3);
  v3(32,1)=6;
  v3(32,2)=846;
  v3(32,3)=  v3(30,3);
  v3(33,1)=6;
  v3(33,2)=1;
  v3(33,3)=(-(T48*getPowerDeriv(y(1),params(4),3)));
  v3(34,1)=6;
  v3(34,2)=1592;
  v3(34,3)=(-(T44*(exp(y(10))*T136+exp(y(10))*T136+y(6)*exp(y(10))*T208)));
  v3(35,1)=6;
  v3(35,2)=920;
  v3(35,3)=  v3(34,3);
  v3(36,1)=6;
  v3(36,2)=968;
  v3(36,3)=  v3(34,3);
  v3(37,1)=6;
  v3(37,2)=1527;
  v3(37,3)=(-(T91*(exp(y(10))*T82+y(6)*exp(y(10))*T136)));
  v3(38,1)=6;
  v3(38,2)=75;
  v3(38,3)=  v3(37,3);
  v3(39,1)=6;
  v3(39,2)=123;
  v3(39,3)=  v3(37,3);
  v3(40,1)=6;
  v3(40,2)=855;
  v3(40,3)=  v3(37,3);
  v3(41,1)=6;
  v3(41,2)=963;
  v3(41,3)=  v3(37,3);
  v3(42,1)=6;
  v3(42,2)=1587;
  v3(42,3)=  v3(37,3);
  v3(43,1)=6;
  v3(43,2)=1522;
  v3(43,3)=(-(y(6)*exp(y(10))*T82*T142));
  v3(44,1)=6;
  v3(44,2)=10;
  v3(44,3)=  v3(43,3);
  v3(45,1)=6;
  v3(45,2)=118;
  v3(45,3)=  v3(43,3);
  v3(46,1)=6;
  v3(46,2)=1644;
  v3(46,3)=(-(T44*(exp(y(10))*T82+y(6)*exp(y(10))*T136+exp(y(10))*y(6)*exp(y(10))*T135+y(6)*exp(y(10))*(T136+y(6)*exp(y(10))*T207))));
  v3(47,1)=6;
  v3(47,2)=972;
  v3(47,3)=  v3(46,3);
  v3(48,1)=6;
  v3(48,2)=1596;
  v3(48,3)=  v3(46,3);
  v3(49,1)=6;
  v3(49,2)=1639;
  v3(49,3)=(-(T91*(y(6)*exp(y(10))*T82+y(6)*exp(y(10))*y(6)*exp(y(10))*T135)));
  v3(50,1)=6;
  v3(50,2)=127;
  v3(50,3)=  v3(49,3);
  v3(51,1)=6;
  v3(51,2)=1531;
  v3(51,3)=  v3(49,3);
  v3(52,1)=6;
  v3(52,2)=1648;
  v3(52,3)=(-(T44*(y(6)*exp(y(10))*T82+y(6)*exp(y(10))*y(6)*exp(y(10))*T135+y(6)*exp(y(10))*y(6)*exp(y(10))*T135+y(6)*exp(y(10))*(y(6)*exp(y(10))*T135+y(6)*exp(y(10))*y(6)*exp(y(10))*T206))));
  g3 = sparse(v3(:,1),v3(:,2),v3(:,3),8,2197);
end
end
end
end
