// Basic RBC Model with log preferences and levels
//
// Jesus Fernandez-Villaverde
// Philadelphia, April 9, 2009

//----------------------------------------------------------------
// 0. Housekeeping
//----------------------------------------------------------------

close all

//----------------------------------------------------------------
// 1. Endogenous variables (8=5+2+1)
//----------------------------------------------------------------

var 

// Allocation (5)
y c i l k 

// Prices (2)
r w

// Productivity (1)
z;

//----------------------------------------------------------------
// 2. Exogenous variables (1)
//----------------------------------------------------------------
 
varexo 

ez;

//----------------------------------------------------------------
// 3. Parameters
//----------------------------------------------------------------

parameters 

// Preferences
bbeta ppsi zzeta

// Technology
aalpha ddelta rrho

// S.D.'s stochastic processes
ssigmaz;

//----------------------------------------------------------------
// 4. Calibration
//----------------------------------------------------------------

// Preferences
bbeta      = 0.99;
zzeta      = 0.5;

// Technology
aalpha     = 0.33;
ddelta     = 0.025;
rrho       = 0.95;  

// Stochastic processes
ssigmaz    = 0.007;

l_ss  = 1/3;
r_ss  = 1/bbeta-1+ddelta;
k_ss  = ((r_ss/(aalpha*(l_ss^(1-aalpha))))^(1/(aalpha-1)));
y_ss  = (k_ss^aalpha)*(l_ss^(1-aalpha));
w_ss  = (1-aalpha)*y_ss/l_ss;
i_ss  = ddelta*k_ss;
c_ss  = y_ss-i_ss;
ppsi  = w_ss/(c_ss*(l_ss^zzeta));

//----------------------------------------------------------------
// 5. Model
//----------------------------------------------------------------

model; 
  
  // 1. Euler equation return on private capital
  1 = bbeta*(c/c(+1))*(1+r(+1)-ddelta);
  
  // 2. Rental rate of capital
  r = aalpha*y/k(-1);
  
  // 3. Static condition leisure-consumption
  ppsi*(l^zzeta)*c = w;

  // 4. Wages
  w = (1-aalpha)*y/l;

  // 5. Law of motion for private capital
  k = (1-ddelta)*k(-1)+i;

  // 6. Production function
  y = (k(-1)^aalpha)*((exp(z)*l)^(1-aalpha));

  // 7. Resource constraint of the economy 
  c+i = y;

  // 8. Productivity process 
  z = rrho*z(-1)+ssigmaz*ez;

end;

//----------------------------------------------------------------
// 6. Computation
//----------------------------------------------------------------

initval;
  y = y_ss;
  c = c_ss;
  i = i_ss;
  l = l_ss;
  k = k_ss;
  r = r_ss;
  w = w_ss;
  z = 0;
  ez = 0;
end;

shocks;
  var ez = 1;
end;

steady;

check;

stoch_simul(hp_filter = 1600, irf = 50, order = 3);