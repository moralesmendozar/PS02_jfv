%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

if isoctave || matlab_ver_less_than('8.6')
    clear all
else
    clearvars -global
    clear_persistent_variables(fileparts(which('dynare')), false)
end
tic0 = tic;
% Save empty dates and dseries objects in memory.
dates('initialize');
dseries('initialize');
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'rbc_log';
M_.dynare_version = '4.5.4';
oo_.dynare_version = '4.5.4';
options_.dynare_version = '4.5.4';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('rbc_log.log');
M_.exo_names = 'ez';
M_.exo_names_tex = 'ez';
M_.exo_names_long = 'ez';
M_.endo_names = 'y';
M_.endo_names_tex = 'y';
M_.endo_names_long = 'y';
M_.endo_names = char(M_.endo_names, 'c');
M_.endo_names_tex = char(M_.endo_names_tex, 'c');
M_.endo_names_long = char(M_.endo_names_long, 'c');
M_.endo_names = char(M_.endo_names, 'i');
M_.endo_names_tex = char(M_.endo_names_tex, 'i');
M_.endo_names_long = char(M_.endo_names_long, 'i');
M_.endo_names = char(M_.endo_names, 'l');
M_.endo_names_tex = char(M_.endo_names_tex, 'l');
M_.endo_names_long = char(M_.endo_names_long, 'l');
M_.endo_names = char(M_.endo_names, 'k');
M_.endo_names_tex = char(M_.endo_names_tex, 'k');
M_.endo_names_long = char(M_.endo_names_long, 'k');
M_.endo_names = char(M_.endo_names, 'r');
M_.endo_names_tex = char(M_.endo_names_tex, 'r');
M_.endo_names_long = char(M_.endo_names_long, 'r');
M_.endo_names = char(M_.endo_names, 'w');
M_.endo_names_tex = char(M_.endo_names_tex, 'w');
M_.endo_names_long = char(M_.endo_names_long, 'w');
M_.endo_names = char(M_.endo_names, 'z');
M_.endo_names_tex = char(M_.endo_names_tex, 'z');
M_.endo_names_long = char(M_.endo_names_long, 'z');
M_.endo_partitions = struct();
M_.param_names = 'bbeta';
M_.param_names_tex = 'bbeta';
M_.param_names_long = 'bbeta';
M_.param_names = char(M_.param_names, 'ppsi');
M_.param_names_tex = char(M_.param_names_tex, 'ppsi');
M_.param_names_long = char(M_.param_names_long, 'ppsi');
M_.param_names = char(M_.param_names, 'zzeta');
M_.param_names_tex = char(M_.param_names_tex, 'zzeta');
M_.param_names_long = char(M_.param_names_long, 'zzeta');
M_.param_names = char(M_.param_names, 'aalpha');
M_.param_names_tex = char(M_.param_names_tex, 'aalpha');
M_.param_names_long = char(M_.param_names_long, 'aalpha');
M_.param_names = char(M_.param_names, 'ddelta');
M_.param_names_tex = char(M_.param_names_tex, 'ddelta');
M_.param_names_long = char(M_.param_names_long, 'ddelta');
M_.param_names = char(M_.param_names, 'rrho');
M_.param_names_tex = char(M_.param_names_tex, 'rrho');
M_.param_names_long = char(M_.param_names_long, 'rrho');
M_.param_names = char(M_.param_names, 'ssigmaz');
M_.param_names_tex = char(M_.param_names_tex, 'ssigmaz');
M_.param_names_long = char(M_.param_names_long, 'ssigmaz');
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 8;
M_.param_nbr = 7;
M_.orig_endo_nbr = 8;
M_.aux_vars = [];
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
M_.det_shocks = [];
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
M_.hessian_eq_zero = 0;
erase_compiled_function('rbc_log_static');
erase_compiled_function('rbc_log_dynamic');
M_.orig_eq_nbr = 8;
M_.eq_nbr = 8;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./' M_.fname '_set_auxiliary_variables.m'], 'file') == 2;
M_.lead_lag_incidence = [
 0 3 0;
 0 4 11;
 0 5 0;
 0 6 0;
 1 7 0;
 0 8 12;
 0 9 0;
 2 10 0;]';
M_.nstatic = 4;
M_.nfwrd   = 2;
M_.npred   = 2;
M_.nboth   = 0;
M_.nsfwrd   = 2;
M_.nspred   = 2;
M_.ndynamic   = 4;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(8, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(7, 1);
M_.NNZDerivatives = [25; 25; 52];
close all
M_.params( 1 ) = 0.99;
bbeta = M_.params( 1 );
M_.params( 3 ) = 0.5;
zzeta = M_.params( 3 );
M_.params( 4 ) = 0.33;
aalpha = M_.params( 4 );
M_.params( 5 ) = 0.025;
ddelta = M_.params( 5 );
M_.params( 6 ) = 0.95;
rrho = M_.params( 6 );
M_.params( 7 ) = 0.007;
ssigmaz = M_.params( 7 );
l_ss  = 1/3;
r_ss  = 1/bbeta-1+ddelta;
k_ss  = ((r_ss/(aalpha*(l_ss^(1-aalpha))))^(1/(aalpha-1)));
y_ss  = (k_ss^aalpha)*(l_ss^(1-aalpha));
w_ss  = (1-aalpha)*y_ss/l_ss;
i_ss  = ddelta*k_ss;
c_ss  = y_ss-i_ss;
M_.params( 2 ) = w_ss/(c_ss*l_ss^M_.params(3));
ppsi = M_.params( 2 );
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 1 ) = y_ss;
oo_.steady_state( 2 ) = c_ss;
oo_.steady_state( 3 ) = i_ss;
oo_.steady_state( 4 ) = l_ss;
oo_.steady_state( 5 ) = k_ss;
oo_.steady_state( 6 ) = r_ss;
oo_.steady_state( 7 ) = w_ss;
oo_.steady_state( 8 ) = 0;
oo_.exo_steady_state( 1 ) = 0;
if M_.exo_nbr > 0
	oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
	oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = 1;
steady;
oo_.dr.eigval = check(M_,options_,oo_);
options_.k_order_solver = 1;
options_.hp_filter = 1600;
options_.irf = 50;
options_.order = 3;
var_list_ = char();
info = stoch_simul(var_list_);
save('rbc_log_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('rbc_log_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('rbc_log_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('rbc_log_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('rbc_log_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('rbc_log_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('rbc_log_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
