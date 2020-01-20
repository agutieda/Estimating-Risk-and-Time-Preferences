% RANDOM MODELS FOR THE JOINT TREATMENT OF RISK AND TIME PREFERENCES
% by Jose Apesteguia, Miguel A. Ballester and Angelo Gutierrez
%
% Script to estimate DEU model using data from:
% "At the nexus of risk and time preferences: An experimental investigation" (2010)
% by Keith H. Coble & Jayson L. Lusk
% Journal of Risk and Uncertainty
%
% This script:
% - Estimate risk and time by individual
%
% November 2019
%
% Tested using Matlab 2019b

clear; close all; clc; rng(1);
warning off;

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

% Time is in weeks so multiply per 7 to transform in days
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

nTasks = 9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NUMERICAL SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Options
min_ra  = -3          ; % Lower bound on risk aversion ra
max_ra  =  1          ; % Upper bound on risk aversion ra
min_da  =  0.5        ; % Lower bound on delay aversion da
max_da  =  1          ; % Upper bound on delay aversion da
nPoints =  50000      ; % Number of points used to compute the integrals
                 
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

% Distribute integration points uniforaly across area
raPoints = ld_sequences(:,1)*(max_ra-min_ra) + min_ra;
daPoints = ld_sequences(:,2)*(max_da-min_da) + min_da;

% Integration weights
intWeights   =  ( (max_ra-min_ra)*(max_da-min_da)/nPoints )  .* ones(nPoints,1);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load pooled estimates and use as initial values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Output/results_joint.mat','mu_ra','sigma_ra','alpha_da','beta_da');

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

max_error      = nan(nIndividuals,1);
pdfCHECK       = nan(nIndividuals,1);
exit_flag      = nan(nIndividuals,1);
fval           = nan(nIndividuals,1);
noSwitch_flag  = zeros(nIndividuals,3);

% Set a waitbar
wait_bar_cluster = waitbar(0,'Computing estimates by individual...',...
    'Name',['Individual Estimates: AHLR']);
update_par = 10;
update_j = round(nIndividuals/update_par);

tic; flag_time=0;
for j = 1:nIndividuals
    
    % Create dataset for j-th individual
    menuID_j  = menuID(subjectID==subjectList(j));
    Y_j       = Y(subjectID==subjectList(j));
    Y_j_risk  = Y_j(riskTaskMenuIdx==1);
    Y_j_time  = Y_j(timeTaskMenuIdx==1);
    Y_j_joint = Y_j(jointTaskMenuIdx==1);
    
    % Check no-switches
    if all(Y_j_risk==1)
        noSwitch_flag(j,1) = 1;
    end
    if all(Y_j_risk==0)
        noSwitch_flag(j,1) = -1;
    end
    
    if all(Y_j_time==1)
        noSwitch_flag(j,2) = 1;
    end
    if all(Y_j_time==0)
        noSwitch_flag(j,2) = -1;
    end
    
    if all(Y_j_joint==1)
        noSwitch_flag(j,3) = 1;
    end
    if all(Y_j_joint==0)
        noSwitch_flag(j,3) = -1;
    end
    
    % Define the objective function to minimize
    obj_fun = @(x) -loglike_MPL_DEU(x,Y_j,menuID_j,algoList);
    
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
    cluster_var   = subjectID; % Not used here
    Cov_theta_hat = RobustVarCov_MPL(@loglike_MPL_DEU, ...
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
    med_ra_list(j)     = median_ra_j     ;
    sd_ra_list(j)      = sd_ra_j        ;
    med_da_list(j) = median_da_j;
    sd_da_list(j)  = sd_da_j;
    std_med_ra_list(j) = std_median_ra_j ;
    std_sd_ra_list(j)  = std_sd_ra_j    ;
    std_med_da_list(j) = std_median_da_j;
    std_sd_da_list(j)  = std_sd_da_j;
    
    % Check probabilities add to one
    [fval_j,~,~,pdfCHECK_j] = loglike_MPL_DEU(theta_hat_j,Y_j,menuID_j,algoList);
    
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

save('Output/results_individual.mat');
