% RANDOM MODELS FOR THE JOINT TREATMENT OF RISK AND TIME PREFERENCES
% by Jose Apesteguia, Miguel A. Ballester and Angelo Gutierrez
%
% Script to estimate DEU-H, PVCE-H and CEPV-H models using data from:
% "Estimating Time Preferences from Convex Budgets"
%  by James Andreoni and Charles Sprenger.
%  American Economic Review, (2012)
%
% This file:
% - Estimate PVCE-H with independent preferences by pooling all observations
%
% November 2019
%
% Tested using Matlab 2019b

clear; close all; clc; rng(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare the data for estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Input/Data_AS.mat');

%%% Choose rounding (multiples of 5 vs multiples of 10)
% Y = Y_5;
% xList = xList_5';

Y = Y_10;
xList = xList_10';

%%% Uncomment following lines to choose a sub-sample
chosenDataset = 2; % (Experiment 1 with no risk, Experiment 2 with risk)
Y = Y(datasetID==chosenDataset);
menuID = obsTable.menuID_session(obsTable.datasetID==chosenDataset);
menuList.a1 = menuList.a1(menuList.datasetID==chosenDataset);
menuList.a2 = menuList.a2(menuList.datasetID==chosenDataset);
menuList.p1 = menuList.p1(menuList.datasetID==chosenDataset);
menuList.p2 = menuList.p2(menuList.datasetID==chosenDataset);
menuList.t1 = menuList.t1(menuList.datasetID==chosenDataset);
menuList.t2 = menuList.t2(menuList.datasetID==chosenDataset);
menuList.menuID = menuID(1:84);
menuList.datasetID = menuList.datasetID(menuList.datasetID==chosenDataset);
subjectID = subjectID(datasetID==chosenDataset);
menuList.xList = xList;

% Identifier for risk conidtions
taskSetID = nan(size(menuID));
taskSetID(menuID<=14) = 1;
taskSetID(menuID>=15 & menuID <=28) = 2;
taskSetID(menuID>=29 & menuID <=42) = 3;
taskSetID(menuID>=43 & menuID <=56) = 4;
taskSetID(menuID>=57 & menuID <=70) = 5;
taskSetID(menuID>=71) = 6;

menuList.taskSetID = nan(84,1);
menuList.taskSetID(1:14)  = 1;
menuList.taskSetID(15:28) = 2;
menuList.taskSetID(29:42) = 3;
menuList.taskSetID(43:56) = 4;
menuList.taskSetID(57:70) = 5;
menuList.taskSetID(71:84) = 6;

nObs_all = length(Y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NUMERICAL SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Settings for integration

% Options
min_ra = -6       ; % Lower limit of the integral on ra
max_ra =  1       ; % Upper limit of the integral on ra
min_da =  0.5     ; % Lower limit of the integral on dr
max_da =  1       ; % Upper limit of the integral on dr
min_is = -1       ; % Lower limit of the integral on is
max_is =  1       ; % Upper limit of the integral on is
nPoints = 200000  ; % Number of evaluation points used to compute the integrals

%%% Initial values used in the estimation
mu_ra_0     =  9.8597    ; % Initial value: mu_ra
sigma_ra_0  =  4.0779    ; % Initial value: sigma_ra
alpha_da_0  =  8.0159    ; % Initial value: mu_dr
beta_da_0   =  0.53372   ; % Initial value: sigma_dr
mu_is_0     =  0.095861  ; % Initial value: mu_is
sigma_is_0  =  0.36075   ; % Initial value: sigma_is

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
ld_sequences = haltonset(3,'Skip',1e3,'Leap',1e2);
ld_sequences = scramble(ld_sequences,'RR2');

% Get the points used for evaluation
ld_sequences = net(ld_sequences,nPoints);

% Distribute integration points uniformly across area
raPoints = ld_sequences(:,1)*(max_ra-min_ra) + min_ra;
daPoints = ld_sequences(:,2)*(max_da-min_da) + min_da;
isPoints = ld_sequences(:,3)*(max_is-min_is) + min_is;

% Integration weights
intWeights =  ( (max_ra-min_ra)*(max_da-min_da)*(max_is-min_is)/nPoints )  .* ones(nPoints,1);

% Store integration points
algoList.raPoints   = raPoints;
algoList.daPoints   = daPoints;
algoList.isPoints   = isPoints;
algoList.intWeights = intWeights;
algoList.max_ra = max_ra ;
algoList.min_ra = min_ra ;
algoList.max_da = max_da ;
algoList.min_da = min_da ;
algoList.max_is = max_is ;
algoList.min_is = min_is ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE AUXILIARY OBJECTS USED IN THE ESTIMATION (PVCE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vectors describing each menu
a1 = menuList.a1 ; a2 = menuList.a2 ; t1 = menuList.t1 ;
p1 = menuList.p1 ; p2 = menuList.p2 ; t2 = menuList.t2 ;
xList = menuList.xList; endowment = menuList.endowment;
nX = length(xList); nL = length(a1);

% Create appropiate matrices for integration
ra = raPoints*ones(1,nX);
da = daPoints*ones(1,nX);
is = isPoints*ones(1,nX);
df = da.^(1-is) ; % Implied discount factor (df)

% Set background consumption
bc = 0;

% Set time (t_da=1: Daily, t=30: Monthly, t=365: Yearly)
t_da = 30;

% Compute allocation chosen for each combination of preferences and lotteries
chosenX = cell(nL,1);
for iMenu = 1:nL
    
    % PVCE-P of the lottery for each combination of (ra,dr,is)
    
    % Payoff of each allocation
    C1 = ones(nPoints,1)*a1(iMenu)*xList;
    C2 = ones(nPoints,1)*(endowment-xList)*a2(iMenu);
    
    % Expected Value for each allocation
    EV1 = p1(iMenu).*(C1+bc).^(1-ra) + (1-p1(iMenu)).*bc.^(1-ra);
    EV2 = p2(iMenu).*(C2+bc).^(1-ra) + (1-p2(iMenu)).*bc.^(1-ra);
    
    % Certainty Equivalent
    CE1 = EV1.^( (1-is)./(1-ra) );
    CE2 = EV2.^( (1-is)./(1-ra) );
    
    % Present value of certainty equivalent
    PVCE1 = ( df.^( t1(iMenu)/ t_da) ).*CE1;
    PVCE2 = ( df.^( t2(iMenu)/ t_da) ).*CE2;
    PVCE  = (PVCE1 + PVCE2).^(1./(1-is));
    
    % Allocation that maximizes value for each parameter configuration (ra,dr)
    [~,idxMax] = max(PVCE,[],2);
    
    % Indicator matrix defining choice in this menu for each parameter configuration (ra,dr)
    chosenX{iMenu} = sparse((1:nPoints)',idxMax,1,nPoints,nX);
    
end

% Store what we need for estimation
algoList.nX = nX;
algoList.nL = nL;
algoList.chosenX = chosenX;

clearvars C1 C2 EV1 EV2 CE1 CE2 PVCE1 PVCE2 PVCE idxMax chosenX
clearvars ra da is df raPoints daPoints isPoints intWeights ld_sequences

% Define the objective function to minimize
obj_fun = @(x) -loglike_CTB_PVCE(x,Y,menuID,algoList);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ESTIMATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Maximize log-likelihood function using a quasi-newton method
X0 = [mu_ra_0; sigma_ra_0; alpha_da_0; beta_da_0; mu_is_0; sigma_is_0];
[theta_hat,fval,exitflag,output,grad,hessian_loglike] = fminunc(obj_fun,X0,options_fminunc);

% Estimated parameters
mu_ra     = theta_hat(1);
sigma_ra  = theta_hat(2);
alpha_da  = theta_hat(3);
beta_da   = theta_hat(4);
mu_is     = theta_hat(5);
sigma_is  = theta_hat(6);

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
mean_is   = MOM(11) ;
median_is = MOM(12) ;
mode_is   = MOM(13) ;
sd_is     = MOM(14) ;
Z_is      = MOM(15) ;
median_df_DELTA= MOM(16) ;
median_df_rr_DELTA = MOM(17) ;

% Compute robust standard errors of estimated parameters
cluster_var   = subjectID;
Cov_theta_hat = RobustVarCov(@loglike_CTB_PVCE, ...
    theta_hat,Y,menuID,-hessian_loglike,cluster_var,algoList,std_type);

% Use delta method to find std. error. of moments of the estimated distributions
gradient_MOM_hat = jacob_prec(@par2mom,theta_hat,1e-5,algoList);
Cov_MOM_hat= gradient_MOM_hat*Cov_theta_hat*gradient_MOM_hat';

% Std. dev. of estimators
std_theta_hat = real( diag( sqrt(Cov_theta_hat) ) )  ;
std_MOM_hat   = real( diag( sqrt(Cov_MOM_hat  ) ) )  ;

% Estimated parameters
std_mu_ra     = std_theta_hat(1);
std_sigma_ra  = std_theta_hat(2);
std_alpha_da  = std_theta_hat(3);
std_beta_da   = std_theta_hat(4);
std_mu_is     = std_theta_hat(5);
std_sigma_is  = std_theta_hat(6);

% Moments of the estimated distributions
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
std_mean_is   = std_MOM_hat(11)  ;
std_median_is = std_MOM_hat(12)  ;
std_mode_is   = std_MOM_hat(13)  ;
std_sd_is     = std_MOM_hat(14)  ;
std_Z_is      = std_MOM_hat(15)  ;
std_median_df_DELTA    = std_MOM_hat(16) ;
std_median_df_rr_DELTA = std_MOM_hat(17) ;


% Loglikelihood at MLE
[log_like_hat,rhoX,rhoCHECK,pdfCHECK] = loglike_CTB_PVCE(theta_hat,Y,menuID,algoList);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('Output/results_baseline_PVCE.mat');
