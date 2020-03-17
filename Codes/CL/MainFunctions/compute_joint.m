% RANDOM MODELS FOR THE JOINT TREATMENT OF RISK AND TIME PREFERENCES
% by Jose Apesteguia, Miguel A. Ballester and Angelo Gutierrez
%
% Script to estimate DEU model using data from:
% "At the nexus of risk and time preferences: An experimental investigation" (2010)
% by Keith H. Coble & Jayson L. Lusk
% Journal of Risk and Uncertainty
%
% This script:
% - Estimate risk and time pooling all observations
%
% March 2020
%
% Tested using Matlab 2019b

clear; close all; clc; rng(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare the data for estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Input/Data_CobleLusk.mat');

% Prepare structure with unique menus in the dataset
menuList.c1_L1 = obsTable.hia(1:94)  ;
menuList.c2_L1 = obsTable.loa(1:94)  ;
menuList.c1_L2 = obsTable.hib(1:94)  ;
menuList.c2_L2 = obsTable.lob(1:94)  ;
menuList.p_L1  = obsTable.phia(1:94) ;
menuList.p_L2  = obsTable.phib(1:94) ;

% Time is in weeks so multiply per 7 days to transform in days
menuList.t_L1  = obsTable.t1(1:94)*7  ;
menuList.t_L2  = obsTable.t2(1:94)*7  ;

% Data on menu faced by each individual and their choice
menuID       = obsTable.menuID;
subjectID    = obsTable.id;
Y            = obsTable.cha;

% Identify each task
taskIDX  = obsTable.session(1:94);
taskSetIDX = nan(94,1);
taskSetIDX(taskIDX==1 | taskIDX ==2)  = 1;
taskSetIDX(taskIDX==3 | taskIDX ==4 | taskIDX ==5 | taskIDX ==6 )  = 2;
taskSetIDX(taskIDX==7 | taskIDX ==8 | taskIDX ==9 )  = 3;

taskID = obsTable.session;
taskSetID = nan(94,1);
taskSetID(taskID==1 | taskID ==2)  = 1;
taskSetID(taskID==3 | taskID ==4 | taskID ==5 | taskID ==6 )  = 2;
taskSetID(taskID==7 | taskID ==8 | taskID ==9 )  = 3;

nTasks = 9;

% Risk tasks
riskTaskMenuIdx = taskSetIDX == 2;
riskTaskIdx = taskSetID == 2;

% Delay aversion tasks
timeTaskMenuIdx  = taskSetIDX == 1;
timeTaskIdx = taskSetID == 1;

% Joint tasks
jointTaskMenuIdx = taskSetIDX == 3;
jointTaskIdx = taskSetID == 3;

% Store to use in the log-likelihood functions
algoList.menuList         = menuList         ;
algoList.riskTaskMenuIdx  = riskTaskMenuIdx  ;
algoList.timeTaskMenuIdx  = timeTaskMenuIdx  ;
algoList.jointTaskMenuIdx = jointTaskMenuIdx ;

nObs_all = length(Y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NUMERICAL SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Settings for integration

% Options
min_ra  = -3          ; % Lower bound on risk aversion ra
max_ra  =  1          ; % Upper bound on risk aversion ra
min_da  =  0.5        ; % Lower bound on delay aversion da
max_da  =  1          ; % Upper bound on delay aversion da
nPoints =  50000      ; % Number of points used to compute the integrals

%%% Initial values used in the estimation
mu_ra_0     =  2.2934      ; % Initial value: mu_ra
sigma_ra_0  =  1.3041      ; % Initial value: sigma_ra
alpha_da_0  =  10.851      ; % Initial value: alpha_da
beta_da_0   =  1.3066      ; % Initial value: beta_da

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

% Halton sequences
ld_sequences = haltonset(2,'Skip',1e3,'Leap',1e2);
ld_sequences = scramble(ld_sequences,'RR2');

% Get the points used for evaluation
ld_sequences = net(ld_sequences,nPoints);

% Distribute integration points uniforaly across area
raPoints = ld_sequences(:,1)*(max_ra-min_ra) + min_ra;
daPoints = ld_sequences(:,2)*(max_da-min_da) + min_da;

% Integration weights
intWeights    =  ( (max_ra-min_ra)*(max_da-min_da)/nPoints )  .* ones(nPoints,1);
intWeights_ra =  ( (max_ra-min_ra)/nPoints )  .* ones(nPoints,1);

% Store integration points
algoList.raPoints = raPoints;
algoList.daPoints = daPoints;
algoList.intWeights = intWeights;
algoList.intWeights_ra = intWeights_ra;

% Store bounds on distributions
algoList.min_ra = min_ra ;
algoList.max_ra = max_ra ;
algoList.min_da = min_da ;
algoList.max_da = max_da ;

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

% Define the objective function to minimize
obj_fun = @(x) -loglike_MPL_DEU(x,Y,menuID,algoList);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ESTIMATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Maximize log-likelihood function using a quasi-newton method
X0 = [mu_ra_0; sigma_ra_0; alpha_da_0; beta_da_0];
[theta_hat,fval,~,~,~,hessian_loglike] = fminunc(obj_fun,X0,options_fminunc);

% Estimated parameters
mu_ra     = theta_hat(1);
sigma_ra  = theta_hat(2);
alpha_da  = theta_hat(3);
beta_da   = theta_hat(4);

% Moments of the estimated distributions
MOM = par2mom(theta_hat,algoList);
mean_ra   = MOM(1)  ;
median_ra = MOM(2)  ;
mode_ra   = MOM(3)  ;
sd_ra     = MOM(4)  ;
Z_ra      = MOM(5)  ;
mean_da   = MOM(6)  ;
median_da = MOM(7)  ;
mode_da   = MOM(8)  ;
sd_da     = MOM(9)  ;
Z_da      = MOM(10) ;
median_df_DELTA= MOM(11) ;
median_df_rr_DELTA = MOM(12) ;

% Compute robust standard errors of estimated parameters
cluster_var   = subjectID;
Cov_theta_hat = RobustVarCov_MPL(@loglike_MPL_DEU, ...
    theta_hat,Y,menuID,-hessian_loglike,cluster_var,algoList,std_type);

% Use delta method to find std. error. of moments of the estimated distributions
gradient_MOM_hat = jacob_prec(@par2mom,theta_hat,1e-5,algoList);
Cov_MOM_hat= gradient_MOM_hat*Cov_theta_hat*gradient_MOM_hat';

% Std. dev. of estimators
std_theta_hat = real( diag( sqrt(Cov_theta_hat) ) )  ;
std_MOM_hat   = real( diag( sqrt(Cov_MOM_hat  ) ) )  ;

% Parameters
std_mu_ra    = std_theta_hat(1);
std_sigma_ra = std_theta_hat(2);
std_alpha_da = std_theta_hat(3);
std_beta_da  = std_theta_hat(4);

% Moments
std_mean_ra   = std_MOM_hat(1)  ;
std_median_ra = std_MOM_hat(2)  ;
std_mode_ra   = std_MOM_hat(3)  ;
std_sd_ra     = std_MOM_hat(4)  ;
std_Z_ra      = std_MOM_hat(5)  ;
std_mean_da   = std_MOM_hat(6)  ;
std_median_da = std_MOM_hat(7)  ;
std_mode_da   = std_MOM_hat(8)  ;
std_sd_da     = std_MOM_hat(9)  ;
std_Z_da      = std_MOM_hat(10) ;
std_median_df_DELTA    = std_MOM_hat(11) ;
std_median_df_rr_DELTA = std_MOM_hat(12) ;

% Loglikelihood at MLE
log_like_hat = loglike_MPL_DEU(theta_hat,Y,menuID,algoList);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('Output/results_joint.mat');
