% RANDOM MODELS FOR THE JOINT TREATMENT OF RISK AND TIME PREFERENCES
% by Jose Apesteguia, Miguel A. Ballester and Angelo Gutierrez
%
% Script to estimate DEU model using data from:
% "Eliciting Risk And Time Preferences"
% by Steffen Andersen, Glenn W. Harrison, Morten I. Lau, And E. Elisabet Rutstrom
% Econometrica, 2008
%
% This script:
% - Estimate risk and time separatelly
%
% November 2019
%
% Tested using Matlab 2019b

clear; close all; clc; rng(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare the data for estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Input/Data_AHLR.mat');

IndividualID = subjectID;
menuList = struct2table(menuList);

% Task Identifier
taskIDX = uniqueLotteryTable.taskID;

taskSetID = nan(100,1);
% Risk tasks
taskSetID( taskIDX<=4 ) = 1;
riskTaskIdx  = obsTable.taskID<=4;

% Delay aversion tasks
taskSetID( taskIDX>=5 ) = 2;
timeTaskIdx  = obsTable.taskID>=5;

nTasks = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NUMERICAL SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Settings for integration

% Options
min_ra  = -2          ; % Lower bound on risk aversion ra
max_ra  =  1          ; % Upper bound on risk aversion ra
min_da  =  0.9        ; % Lower bound on delay aversion da
max_da  =  1          ; % Upper bound on delay aversion da
nPoints =  50000      ; % Number of points used to compute the integrals

%%% Initial values used in the estimation
mu_ra_0     =  75.092     ; % Initial value: mu_ra
sigma_ra_0  =  6.4249     ; % Initial value: sigma_ra
alpha_da_0  =  74.545     ; % Initial value: alpha_da
beta_da_0   =  1.5992     ; % Initial value: beta_da

%%% Compuation of standard errors
std_type = 'robust'; % Cluster standard errors by individual

%%% Settings for optimization routines used

% Estimation
options_fminunc = optimoptions('fminunc',...
    'Display','iter',...
    'OptimalityTolerance',1e-4,...
    'StepTolerance',1e-4,...
    'MaxFunctionEvaluations',5000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GET THINGS READY FOR THE ESTIMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Draws on support of ra and da

% Halton sequences
ld_sequences = haltonset(1,'Skip',1e3,'Leap',1e2);
ld_sequences = scramble(ld_sequences,'RR2');

% Get the points used for evaluation
ld_sequences = net(ld_sequences,nPoints);

% Distribute integration points uniforaly across area
raPoints = ld_sequences*(max_ra-min_ra) + min_ra;
daPoints = ld_sequences*(max_da-min_da) + min_da;

% Integration weights
intWeights_ra =  ( (max_ra-min_ra)/nPoints )  .* ones(nPoints,1);
intWeights_da =  ( (max_da-min_da)/nPoints )  .* ones(nPoints,1);

% Store
algoList.min_ra = min_ra ;
algoList.max_ra = max_ra ;
algoList.raPoints = raPoints;
algoList.intWeights_ra = intWeights_ra;

algoList.min_da = min_da ;
algoList.max_da = max_da ;
algoList.daPoints = daPoints;
algoList.intWeights_da = intWeights_da;

algoList.menuList = menuList;
algoList.taskSetID = taskSetID;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE AUXILIARY OBJECTS USED IN THE ESTIMATION (DEU)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vectors describing each menu
c1_L1 = menuList.c1_L1; c2_L1 = menuList.c2_L1;
c1_L2 = menuList.c1_L2; c2_L2 = menuList.c2_L2;
t_L1  = menuList.t_L1 ; t_L2  = menuList.t_L2 ;
p_L1  = menuList.p_L1 ; p_L2  = menuList.p_L2 ;
nL = length(c1_L1);

% Create appropiate matrices for integration
ra = raPoints   ; % Draws on the support of risk aversion
da = daPoints   ; % Draws on the support of delay aversion
df = da.^(1-ra) ; % Implied discount factor (df)

% Compute allocation chosen for each combination of preferences and lotteries
chosenX = nan(nPoints,nL);
for iMenu = 1:nL
    
    % Utility in each outcome
    U1_L1 = (c1_L1(iMenu).^(1-ra))./(1-ra);
    U2_L1 = (c2_L1(iMenu).^(1-ra))./(1-ra);
    U1_L2 = (c1_L2(iMenu).^(1-ra))./(1-ra);
    U2_L2 = (c2_L2(iMenu).^(1-ra))./(1-ra);
    
    % Expected utility of each lottery
    EU_L1 = p_L1(iMenu).*U1_L1 + (1-p_L1(iMenu)).*U2_L1;
    EU_L2 = p_L2(iMenu).*U1_L2 + (1-p_L2(iMenu)).*U2_L2;
    
    % Choice on this menu given preferences (ra,dr)
    if t_L1(iMenu) == t_L2(iMenu)
        chosenX(:,iMenu) = EU_L1>=EU_L2;
    else
        DEU_L1 = ( df.^(t_L1(iMenu)/30) ).*EU_L1;
        DEU_L2 = ( df.^(t_L2(iMenu)/30) ).*EU_L2;
        chosenX(:,iMenu) = DEU_L1>=DEU_L2 ;
    end
    
end


% Store what we need for estimation
algoList.nL = nL;
algoList.chosenX = chosenX;
algoList.t_L1 = t_L1;
algoList.t_L2 = t_L2;

clearvars ra da df U1_L1 U2_L1 U1_L2 U2_L2 EU_L1 EU_L2 DEU_L1 DEU_L2 chosenX


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ESTIMATION: RISK AVERSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data with risky tasks
Y_ra = Y(riskTaskIdx);
menuID_ra = menuID(riskTaskIdx);
subjectID_ra = subjectID(riskTaskIdx);
nObs_ra = length(Y_ra);

% Define the objective function to minimize
obj_fun = @(x) -loglike_MPL_RISK(x,Y_ra,menuID_ra,algoList);

% Maximize log-likelihood function using a quasi-newton method
X0 = [mu_ra_0; sigma_ra_0];
[theta_hat1,fval,~,~,~,hessian_loglike] = fminunc(obj_fun,X0,options_fminunc);

% Estimated parameters
mu_ra     = theta_hat1(1);
sigma_ra  = theta_hat1(2);

% Moments of the estimated distributions
MOM = par2mom_ra(theta_hat1,algoList);
mean_ra   = MOM(1)  ;
median_ra = MOM(2)  ;
mode_ra   = MOM(3)  ;
sd_ra     = MOM(4)  ;
Z_ra      = MOM(5)  ;

% Compute robust standard errors of estimated parameters
cluster_var   = subjectID_ra;
Cov_theta_hat1 = RobustVarCov_MPL(@loglike_MPL_RISK, ...
    theta_hat1,Y_ra,menuID_ra,-hessian_loglike,cluster_var,algoList,std_type);

% Use delta method to find std. error. of moments of the estimated distributions
gradient_MOM_hat = jacob_prec(@par2mom_ra,theta_hat1,1e-5,algoList);
Cov_MOM_hat1 = gradient_MOM_hat*Cov_theta_hat1*gradient_MOM_hat';

% Std. dev. of estimators
std_theta_hat = real( diag( sqrt(Cov_theta_hat1) ) )  ;
std_MOM_hat   = real( diag( sqrt(Cov_MOM_hat1  ) ) )  ;

% Parameters
std_mu_ra    = std_theta_hat(1);
std_sigma_ra = std_theta_hat(2);

% Moments
std_mean_ra   = std_MOM_hat(1)  ;
std_median_ra = std_MOM_hat(2)  ;
std_mode_ra   = std_MOM_hat(3)  ;
std_sd_ra     = std_MOM_hat(4)  ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ESTIMATION: DELAY AVERSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data with risky tasks
Y_da = Y(timeTaskIdx);
menuID_da = menuID(timeTaskIdx);
subjectID_da = subjectID(timeTaskIdx);
nObs_da = length(Y_da);

% Define the objective function to minimize
obj_fun = @(x) -loglike_MPL_TIME(x,Y_da,menuID_da,algoList);

% Maximize log-likelihood function using a quasi-newton method
X0 = [alpha_da_0; beta_da_0];
[theta_hat2,fval,~,~,~,hessian_loglike] = fminunc(obj_fun,X0,options_fminunc);

% Estimated parameters
alpha_da = theta_hat2(1);
beta_da  = theta_hat2(2);

% Moments of the estimated distributions
MOM = par2mom_da(theta_hat2,algoList);
mean_da   = MOM(1)  ;
median_da = MOM(2)  ;
mode_da   = MOM(3)  ;
sd_da     = MOM(4)  ;
Z_da      = MOM(5)  ;

% Compute robust standard errors of estimated parameters
cluster_var   = subjectID_da;
Cov_theta_hat2 = RobustVarCov_MPL(@loglike_MPL_TIME, ...
    theta_hat2,Y_da,menuID_da,-hessian_loglike,cluster_var,algoList,std_type);

% Use delta method to find std. error. of moments of the estimated distributions
gradient_MOM_hat = jacob_prec(@par2mom_da,theta_hat2,1e-5,algoList);
Cov_MOM_hat2 = gradient_MOM_hat*Cov_theta_hat2*gradient_MOM_hat';

% Std. dev. of estimators
std_theta_hat = real( diag( sqrt(Cov_theta_hat2) ) )  ;
std_MOM_hat   = real( diag( sqrt(Cov_MOM_hat2  ) ) )  ;

% Parameters
std_alpha_da = std_theta_hat(1);
std_beta_da  = std_theta_hat(2);

% Moments
std_mean_da   = std_MOM_hat(1)  ;
std_median_da = std_MOM_hat(2)  ;
std_mode_da   = std_MOM_hat(3)  ;
std_sd_da     = std_MOM_hat(4)  ;

% Point estimates and STD of discount factor
MOM = par2mom_df([median_ra;median_da]);
median_df_DELTA    = MOM(1) ;
median_df_rr_DELTA = MOM(2) ;

% Use delta method to find std. error. of moments of the estimated parameters
Cov_theta_hat = [std_median_ra^2,0;0,std_median_da^2];
gradient_MOM_hat = jacob_prec(@par2mom_df,[median_ra;median_da],1e-5);
Cov_MOM_hat = gradient_MOM_hat*Cov_theta_hat*gradient_MOM_hat';
std_MOM_hat = real( diag( sqrt(Cov_MOM_hat  ) ) )  ;

std_median_df_DELTA    = std_MOM_hat(1) ;
std_median_df_rr_DELTA = std_MOM_hat(2) ;

% Loglikelihood at MLE
log_like_hat_ra = loglike_MPL_RISK([mu_ra;sigma_ra],Y_ra,menuID_ra,algoList);
log_like_hat_da = loglike_MPL_TIME([alpha_da;beta_da],Y_da,menuID_da,algoList);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('Output/results_separate.mat');
