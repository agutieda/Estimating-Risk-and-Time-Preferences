% RANDOM MODELS FOR THE JOINT TREATMENT OF RISK AND TIME PREFERENCES
% by Jose Apesteguia, Miguel A. Ballester and Angelo Gutierrez
%
% Script to estimate DEU-H model using data from:
% "Estimating Time Preferences from Convex Budgets"
%  by James Andreoni and Charles Sprenger.
%  American Economic Review, (2012)
%
% This file:
% - Estimate DEU-H with independent preferences by individual
%
% March 2020
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
min_ra = -2       ; % Lower limit of the integral on ra
max_ra =  1       ; % Upper limit of the integral on ra
min_da =  0       ; % Lower limit of the integral on da
max_da =  1       ; % Upper limit of the integral on da
nPoints = 50000   ; % Number of evaluation points used to compute the integrals

%%% Initial values used in the estimation
mu_ra_0     =  0.044798   ; % Initial value: mu_ra
sigma_ra_0  =  0.4139     ; % Initial value: sigma_ra
alpha_da_0  =  5.0009     ; % Initial value: alpha_da
beta_da_0   =  0.50126    ; % Initial value: beta_da

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
ld_sequences = haltonset(2,'Skip',1e3,'Leap',1e2);
ld_sequences = scramble(ld_sequences,'RR2');

% Get the points used for evaluation
ld_sequences = net(ld_sequences,nPoints);

% Distribute integration points uniformly across area
raPoints = ld_sequences(:,1)*(max_ra-min_ra) + min_ra;
daPoints = ld_sequences(:,2)*(max_da-min_da) + min_da;

% Integration weights
intWeights =  ( (max_ra-min_ra)*(max_da-min_da)/nPoints )  .* ones(nPoints,1);

% Store integration points
algoList.raPoints   = raPoints;
algoList.daPoints   = daPoints;
algoList.intWeights = intWeights;
algoList.max_ra = max_ra ;
algoList.min_ra = min_ra ;
algoList.max_da = max_da ;
algoList.min_da = min_da ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE AUXILIARY OBJECTS USED IN THE ESTIMATION (DEU)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vectors describing each menu
a1 = menuList.a1 ; a2 = menuList.a2 ; t1 = menuList.t1 ;
p1 = menuList.p1 ; p2 = menuList.p2 ; t2 = menuList.t2 ;
xList = menuList.xList; endowment = menuList.endowment;
nX = length(xList); nL = length(a1);

% Create appropiate matrices for integration
ra = raPoints*ones(1,nX);
da = daPoints*ones(1,nX);
df = da.^(1-ra) ; % Implied discount factor (df)

% Set background consumption
bc = 0;

% Set time (t_da=1: Daily, t=30: Monthly, t=365: Yearly)
t_da = 30;

% Compute allocation chosen for each combination of preferences and lotteries
chosenX = cell(nL,1);
for iMenu = 1:nL
    
    % DEU-P of the lottery for each combination of (ra,da)
    
    % Payoff of each allocation
    C1 = ones(nPoints,1)*a1(iMenu)*xList;
    C2 = ones(nPoints,1)*(endowment-xList)*a2(iMenu);
    
    % Expected utility of each allocation
    EU1 = p1(iMenu) .* ( (C1+bc).^(1-ra) )./(1-ra) + (1-p1(iMenu)) .* ( bc.^(1-ra) )./(1-ra);
    EU2 = p2(iMenu) .* ( (C2+bc).^(1-ra) )./(1-ra) + (1-p2(iMenu)) .* ( bc.^(1-ra) )./(1-ra);
    
    % Discounted expected utility of each allocation
    DEU1 = ( df.^( t1(iMenu) / t_da ) ).*EU1;
    DEU2 = ( df.^( t2(iMenu) / t_da ) ).*EU2;
    DEU  = DEU1 + DEU2;
    
    % Allocation that maximizes value for each parameter configuration (ra,da)
    [~,idxMax] = max(DEU,[],2);
    
    % Indicator matrix defining choice in this menu for each parameter configuration (ra,da)
    chosenX{iMenu} = sparse((1:nPoints)',idxMax,1,nPoints,nX);
    
end

% Store what we need for estimation
algoList.nX = nX;
algoList.nL = nL;
algoList.chosenX = chosenX;

clearvars C1 C2 EU1 EU2 DEU1 DEU2 DEU idxMax chosenX
clearvars ra da df raPoints daPoints intWeights ld_sequences

% Define the objective function to minimize
obj_fun = @(x) -loglike_CTB_DEU(x,Y,menuID,algoList);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load pooled estimates and use as initial values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Output/results_baseline_DEU.mat','mu_ra','sigma_ra','alpha_da','beta_da');
X0 = [mu_ra; sigma_ra; alpha_da; beta_da];

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

std_med_ra_list  = nan(nIndividuals,1);
std_sd_ra_list   = nan(nIndividuals,1);
std_med_da_list  = nan(nIndividuals,1);
std_sd_da_list   = nan(nIndividuals,1);

theta_list     = nan(nIndividuals,4);
max_error      = nan(nIndividuals,1);
pdfCHECK       = nan(nIndividuals,1);
exit_flag      = nan(nIndividuals,1);
fval           = nan(nIndividuals,1);
noSwitch_flag  = zeros(nIndividuals,1);

% Set a waitbar
wait_bar_cluster = waitbar(0,'Computing estimates by individual...',...
    'Name',['Individual Estimates: AS-DEU']);
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
    obj_fun = @(x) -loglike_CTB_DEU(x,Y_j,menuID_j,algoList);
    
    % Maximize log-likelihood function using a quasi-newton method
    [theta_hat_j,~,exitflag_j,~,~,hessian_loglike_j] = fminunc(obj_fun,X0,options_fminunc);
    
    % Estimated parameters
    mu_ra_j     = theta_hat_j(1);
    sigma_ra_j  = theta_hat_j(2);
    alpha_da_j  = theta_hat_j(3);
    beta_da_j   = theta_hat_j(4);
    
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
    median_df_DELTA_j = MOM(11) ;
    median_df_rr_DELTA_j = MOM(12) ;
    
    % Compute robust standard errors of estimated parameters
    cluster_var   = subjectID(subjectID==subjectList(j)); % Not used here
    Cov_theta_hat = RobustVarCov(@loglike_CTB_DEU, ...
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
    std_median_df_DELTA_j    = std_MOM_hat(11) ;
    std_median_df_rr_DELTA_j = std_MOM_hat(12) ;
    
    % Store
    theta_list(j,:)    = theta_hat_j;
    med_ra_list(j)     = median_ra_j     ;
    sd_ra_list(j)      = sd_ra_j         ;
    med_da_list(j)     = median_da_j     ;
    sd_da_list(j)      = sd_da_j         ;
    std_med_ra_list(j) = std_median_ra_j ;
    std_sd_ra_list(j)  = std_sd_ra_j     ;
    std_med_da_list(j) = std_median_da_j ;
    std_sd_da_list(j)  = std_sd_da_j     ;
    
    % Check probabilities add to one
    [fval_j,~,~,pdfCHECK_j] = loglike_CTB_DEU(theta_hat_j,Y_j,menuID_j,algoList);
    
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

save('Output/results_individual_DEU.mat');
