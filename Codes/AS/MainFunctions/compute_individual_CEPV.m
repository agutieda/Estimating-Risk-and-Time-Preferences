% RANDOM MODELS FOR THE JOINT TREATMENT OF RISK AND TIME PREFERENCES
% by Jose Apesteguia, Miguel A. Ballester and Angelo Gutierrez
%
% Script to estimate DEU-H, PVCE-H and CEPV-H models using data from:
% "Estimating Time Preferences from Convex Budgets"
%  by James Andreoni and Charles Sprenger.
%  American Economic Review, (2012)
%
% This file:
% - Estimate CEPV-H with independent preferences by individual
%
% November 2019
%
% Tested using Matlab 2019b

clear; close all; clc; rng(1);
warning off;

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

% Options
min_ra = -5       ; % Lower limit of the integral on ra
max_ra =  1       ; % Upper limit of the integral on ra
min_da =  0.5     ; % Lower limit of the integral on dr
max_da =  1       ; % Upper limit of the integral on dr
min_is = -3       ; % Lower limit of the integral on is
max_is =  1       ; % Upper limit of the integral on is
nPoints = 200000  ; % Number of evaluation points used to compute the integrals

%%% Initial values used in the estimation
mu_ra_0     =  27.919   ; % Initial value: mu_ra
sigma_ra_0  =  5.5052   ; % Initial value: sigma_ra
alpha_da_0  =  7.2193   ; % Initial value: mu_dr
beta_da_0   =  0.51752  ; % Initial value: sigma_dr
mu_is_0     = -0.13344  ; % Initial value: mu_is
sigma_is_0  =  0.68189  ; % Initial value: sigma_is

%%% Settings for optimization routines used

% Estimation
options_fminunc = optimoptions('fminunc',...
    'Display','off',...
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
%% COMPUTE AUXILIARY OBJECTS USED IN THE ESTIMATION (CEPV)
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
    
    % CEPV-P of the lottery for each combination of (ra,dr,is)
    
    % Payoff of each allocation
    C1 = ones(nPoints,1)*a1(iMenu)*xList;
    C2 = ones(nPoints,1)*(endowment-xList)*a2(iMenu);
    
    % Value of each possible outcome
    V0 = (    bc   ).^(1-is) ; % No lottery pays
    V1 = ( C1 + bc ).^(1-is) ; % Lottery in t1 pays
    V2 = ( C2 + bc ).^(1-is) ; % Lottery in t2 pays
    
    % Present value of each possible outcome
    PV11 = ( df.^(t1(iMenu)/t_da) ).*V1 + ( df.^(t2(iMenu)/t_da) ).*V2 ; % Both lotteries pay
    PV10 = ( df.^(t1(iMenu)/t_da) ).*V1 + ( df.^(t2(iMenu)/t_da) ).*V0 ; % Only lottery in t1 pays
    PV01 = ( df.^(t1(iMenu)/t_da) ).*V0 + ( df.^(t2(iMenu)/t_da) ).*V2 ; % Only lottery in t2 pays
    PV00 = ( df.^(t1(iMenu)/t_da) ).*V0 + ( df.^(t2(iMenu)/t_da) ).*V0 ; % No lottery pays
    
    % Certainty equivalent of the present value
    CEPV11 = (  p1(iMenu)  ).*(  p2(iMenu)  ).*PV11.^( (1-ra)./(1-is) ) ;
    CEPV10 = (  p1(iMenu)  ).*( 1-p2(iMenu) ).*PV10.^( (1-ra)./(1-is) ) ;
    CEPV01 = ( 1-p1(iMenu) ).*(  p2(iMenu)  ).*PV01.^( (1-ra)./(1-is) ) ;
    CEPV00 = ( 1-p1(iMenu) ).*( 1-p2(iMenu) ).*PV00.^( (1-ra)./(1-is) ) ;
    CEPV   = (CEPV11 + CEPV10 + CEPV01 + CEPV00).^(1./(1-ra)) ;
    
    % Allocation that maximizes value for each parameter configuration (ra,dr)
    [~,idxMax] = max(CEPV,[],2);
    
    % Indicator matrix defining choice in this menu for each parameter configuration (ra,dr)
    chosenX{iMenu} = sparse((1:nPoints)',idxMax,1,nPoints,nX);
    
end

% Store what we need for estimation
algoList.nX = nX;
algoList.nL = nL;
algoList.chosenX = chosenX;

clearvars C1 C2 V0 V1 V2 PV11 PV10 PV01 PV00 CEPV11 CEPV10 CEPV01 CEPV00 CEPV idxMax chosenX
clearvars ra da is df raPoints daPoints isPoints intWeights ld_sequences

% Define the objective function to minimize
obj_fun = @(x) -loglike_CTB_CEPV(x,Y,menuID,algoList);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load pooled estimates and use as initial values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Output/results_baseline_CEPV.mat','mu_ra','sigma_ra','alpha_da','beta_da','mu_is','sigma_is');
X0 = [mu_ra; sigma_ra; alpha_da; beta_da; mu_is; sigma_is ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEU: BY INDIVIDUAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% List of unique IDs
subjectList = unique(subjectID);

% Prelocate
nIndividuals = length(subjectList);
med_ra_list  = nan(nIndividuals,1);
sd_ra_list   = nan(nIndividuals,1);
med_da_list  = nan(nIndividuals,1);
sd_da_list   = nan(nIndividuals,1);
med_is_list  = nan(nIndividuals,1);
sd_is_list   = nan(nIndividuals,1);

std_med_ra_list  = nan(nIndividuals,1);
std_sd_ra_list   = nan(nIndividuals,1);
std_med_da_list  = nan(nIndividuals,1);
std_sd_da_list   = nan(nIndividuals,1);
std_med_is_list  = nan(nIndividuals,1);
std_sd_is_list   = nan(nIndividuals,1);

theta_list     = nan(nIndividuals,6);
max_error      = nan(nIndividuals,1);
pdfCHECK       = nan(nIndividuals,1);
exit_flag      = nan(nIndividuals,1);
fval           = nan(nIndividuals,1);
noSwitch_flag  = zeros(nIndividuals,1);

% Set a waitbar
wait_bar_cluster = waitbar(0,'Computing estimates by individual...',...
    'Name',['Individual Estimates: AS-CEPV']);
update_par = 10;
update_j = round(nIndividuals/update_par);

tic; flag_time=0;
for j = 1:nIndividuals
    
    % Create dataset for j-th individual
    menuID_j = menuID(subjectID==subjectList(j));
    Y_j      = Y(subjectID==subjectList(j));
    
    % Check no-switches
    if all(Y_j==length(xList))
        noSwitch_flag(j,1) = 1; % x=0: Allocates all tokens in the future
    end
    
    if all(Y_j==1)
        noSwitch_flag(j,2) = -1; % x=100: Allocates all tokens in the present
    end
    
    % Define the objective function to minimize
    obj_fun = @(x) -loglike_CTB_CEPV(x,Y_j,menuID_j,algoList);
    
    % Maximize log-likelihood function using a quasi-newton method
    [theta_hat_j,~,exitflag_j,~,~,hessian_loglike_j] = fminunc(obj_fun,X0,options_fminunc);
    
    % Estimated parameters
    mu_ra_j     = theta_hat_j(1);
    sigma_ra_j  = theta_hat_j(2);
    alpha_da_j  = theta_hat_j(3);
    beta_da_j   = theta_hat_j(4);
    mu_is_j     = theta_hat_j(5);
    sigma_is_j  = theta_hat_j(6);
    
    
    % Moments of the estimated distributions
    MOM = par2mom(theta_hat_j,algoList);
    mean_ra_j   = MOM(1)  ;
    median_ra_j = MOM(2)  ;
    mode_ra_j   = MOM(3)  ;
    sd_ra_j     = MOM(4)  ;
    Z_ra_j      = MOM(5)  ;
    mean_da_j   = MOM(6)  ;
    median_da_j = MOM(7)  ;
    mode_da_j   = MOM(8)  ;
    sd_da_j     = MOM(9)  ;
    Z_da_j      = MOM(10) ;
    mean_is_j   = MOM(11) ;
    median_is_j = MOM(12) ;
    mode_is_j   = MOM(13) ;
    sd_is_j     = MOM(14) ;
    Z_is_j      = MOM(15) ;
    median_df_DELTA_j = MOM(16) ;
    median_df_rr_DELTA_j = MOM(17) ;
    
    % Compute robust standard errors of estimated parameters
    cluster_var   = subjectID(subjectID==subjectList(j)); % Not used here
    Cov_theta_hat = RobustVarCov(@loglike_CTB_CEPV, ...
        theta_hat_j,Y_j,menuID_j,-hessian_loglike_j,cluster_var,algoList,'simple');
    
    % Use delta method to find std. error. of moments of the estimated distributions
    gradient_MOM_hat = jacob_prec(@par2mom,theta_hat_j,1e-5,algoList);
    Cov_MOM_hat= gradient_MOM_hat*Cov_theta_hat*gradient_MOM_hat';
    
    % Std. dev. of estimators
    std_theta_hat = real( diag( sqrt(Cov_theta_hat) ) )  ;
    std_MOM_hat   = real( diag( sqrt(Cov_MOM_hat  ) ) )  ;
    
    % Parameters
    std_mu_ra_j    = std_theta_hat(1);
    std_sigma_ra_j = std_theta_hat(2);
    std_alpha_da_j = std_theta_hat(3);
    std_beta_da_j  = std_theta_hat(4);
    
    % Moments
    std_mean_ra_j   = std_MOM_hat(1)  ;
    std_median_ra_j = std_MOM_hat(2)  ;
    std_mode_ra_j   = std_MOM_hat(3)  ;
    std_sd_ra_j     = std_MOM_hat(4)  ;
    std_Z_ra_j      = std_MOM_hat(5)  ;
    std_mean_da_j   = std_MOM_hat(6)  ;
    std_median_da_j = std_MOM_hat(7)  ;
    std_mode_da_j   = std_MOM_hat(8)  ;
    std_sd_da_j     = std_MOM_hat(9)  ;
    std_Z_da_j      = std_MOM_hat(10) ;
    std_mean_is_j   = std_MOM_hat(11) ;
    std_median_is_j = std_MOM_hat(12) ;
    std_mode_is_j   = std_MOM_hat(13) ;
    std_sd_is_j     = std_MOM_hat(14) ;
    std_Z_is_j      = std_MOM_hat(15) ;
    std_median_df_DELTA_j    = std_MOM_hat(11) ;
    std_median_df_rr_DELTA_j = std_MOM_hat(12) ;
    
    % Store
    theta_list(j,:)    = theta_hat_j;
    med_ra_list(j)     = median_ra_j     ;
    sd_ra_list(j)      = sd_ra_j         ;
    med_da_list(j)     = median_da_j     ;
    sd_da_list(j)      = sd_da_j         ;
    med_is_list(j)     = median_is_j     ;
    sd_is_list(j)      = sd_is_j         ;
    std_med_ra_list(j) = std_median_ra_j ;
    std_sd_ra_list(j)  = std_sd_ra_j     ;
    std_med_da_list(j) = std_median_da_j ;
    std_sd_da_list(j)  = std_sd_da_j     ;
    std_med_is_list(j) = std_median_is_j ;
    std_sd_is_list(j)  = std_sd_is_j     ;
    
    % Check probabilities add to one
    [fval_j,~,~,pdfCHECK_j] = loglike_CTB_PVCE(theta_hat_j,Y_j,menuID_j,algoList);
    
    % Other
    exit_flag(j,1) = exitflag_j;
    fval(j,1)      = fval_j;
    pdfCHECK(j,1)  = pdfCHECK_j;
    
    % Update waitbar
    if mod(j,update_j)==0
        
        if flag_time==0
            time_by_round = toc;
            flag_time = 1;
        end
        
        update_par = update_par-1;
        timeLeft = round(time_by_round*update_par/60,2);
        
        waitbar(j/nIndividuals,wait_bar_cluster,['Estimated time left: ' num2str(timeLeft),' min']);
        
    end
    
    
end
close(wait_bar_cluster);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('Output/results_individual_CEPV.mat');
