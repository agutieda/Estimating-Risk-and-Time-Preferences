% RANDOM MODELS FOR THE JOINT TREATMENT OF RISK AND TIME PREFERENCES
% by Jose Apesteguia, Miguel A. Ballester and Angelo Gutierrez
%
% Script to estimate DEU model using data from:
% "Eliciting Risk And Time Preferences"
% by Steffen Andersen, Glenn W. Harrison, Morten I. Lau, And E. Elisabet Rutstrom
% Econometrica, 2008
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

load('Input/Data_AHLR.mat');

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

% Data with risky tasks
Y_ra = Y(riskTaskIdx);
menuID_ra = menuID(riskTaskIdx);
subjectID_ra = subjectID(riskTaskIdx);

% Data with discount tasks
Y_da = Y(timeTaskIdx);
menuID_da = menuID(timeTaskIdx);
subjectID_da = subjectID(timeTaskIdx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NUMERICAL SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Options
min_ra  = -2          ; % Lower bound on risk aversion ra
max_ra  =  1          ; % Upper bound on risk aversion ra
min_da  =  0.9        ; % Lower bound on delay aversion da
max_da  =  1          ; % Upper bound on delay aversion da
nPoints =  100000     ; % Number of points used to compute the integrals

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
%% Load pooled estimates and use as initial values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Output/results_separate.mat','mu_ra','sigma_ra','alpha_da','beta_da');
X0_ra = [mu_ra; sigma_ra];
X0_da = [alpha_da; beta_da];

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
pdfCHECK       = nan(nIndividuals,2);
exit_flag      = nan(nIndividuals,2);
fval           = nan(nIndividuals,2);
noSwitch_flag  = zeros(nIndividuals,2);

% Set a waitbar
wait_bar_cluster = waitbar(0,'Computing estimates by individual...',...
    'Name',['Individual Estimates: AHLR']);
update_par = 10;
update_j = round(nIndividuals/update_par);

tic; flag_time=0;
for j = 1:nIndividuals
    
    % Create dataset for j-th individual
    menuID_j = menuID(subjectID==subjectList(j));
    Y_j      = Y(subjectID==subjectList(j));
    riskTaskIdx_j = riskTaskIdx(subjectID==subjectList(j));
    timeTaskIdx_j = timeTaskIdx(subjectID==subjectList(j));
    fval_DCE_j_risk = nan;
    exitflag_j_risk = nan;
    fval_DCE_j_time = nan;
    exitflag_j_time = nan;
    
    %% Estimate risk preferences
    
    % Create dataset with subset of time tasks only
    menuID_j_risk = menuID_j(riskTaskIdx_j);
    Y_j_risk      = Y_j(riskTaskIdx_j);
    
    if ~isempty(menuID_j_risk)
        
        % Check no switches
        
        auxY_j = Y_j_risk(Y_j_risk~=-1);
        
        if all(auxY_j==1)
            noSwitch_flag(j,1) = 1;
        end
        
        if all(auxY_j==0)
            noSwitch_flag(j,1) = -1;
        end
        
        % Define the objective function to minimize
        obj_fun = @(x) -loglike_MPL_RISK(x,Y_j_risk,menuID_j_risk,algoList);
        
        % Estimate
        [theta_hat_j,fval_j_risk,exitflag_j_risk,~,~,hessian_loglike_j] = fminunc(obj_fun,X0_ra,options_fminunc);
        
        % Estimated parameters
        mu_ra_j    = theta_hat_j(1);
        sigma_ra_j = theta_hat_j(2);
        
        % Moments of the estimated distributions
        MOM_j = par2mom_ra(theta_hat_j,algoList);
        mean_ra_j   = MOM_j(1)  ;
        median_ra_j = MOM_j(2)  ;
        mode_ra_j   = MOM_j(3)  ;
        sd_ra_j     = MOM_j(4)  ;
        Z_ra_j      = MOM_j(5)  ;
        
        % Compute standard errors
        cluster_var = []; % Not used here so it doesn't matter
        Cov_theta_hat_j = RobustVarCov_MPL(@loglike_MPL_RISK, ...
            theta_hat_j,Y_j_risk,menuID_j_risk,-hessian_loglike_j,cluster_var,algoList,'simple');
        
        % Get the standard deviations from the Covariance matrix
        std_x_hat_j = diag(sqrt(Cov_theta_hat_j));
        
        % Use delta method to find std. error. of moments of the estimated distributions
        gradient_MOM_hat_j = jacob_prec(@par2mom_ra,theta_hat_j,1e-5,algoList);
        Cov_MOM_hat_j = gradient_MOM_hat_j*Cov_theta_hat_j*gradient_MOM_hat_j';
        
        % Std. dev. of estimators
        std_theta_hat_j = real( diag( sqrt(Cov_theta_hat_j) ) )  ;
        std_MOM_hat_j   = real( diag( sqrt(Cov_MOM_hat_j  ) ) )  ;
        
        % Parameters
        std_mu_ra_j    = std_theta_hat_j(1);
        std_sigma_ra_j = std_theta_hat_j(2);
        
        % Moments
        std_mean_ra_j   = std_MOM_hat_j(1)  ;
        std_median_ra_j = std_MOM_hat_j(2)  ;
        std_mode_ra_j   = std_MOM_hat_j(3)  ;
        std_sd_ra_j     = std_MOM_hat_j(4)  ;
        
        % Store estimated parameters
        med_ra_list(j)     = median_ra_j     ;
        sd_ra_list(j)      = sd_ra_j        ;
        std_med_ra_list(j) = std_median_ra_j ;
        std_sd_ra_list(j)  = std_sd_ra_j    ;
        
        % Check probabilities add to one
        [fval_j_risk,~,~,pdfCHECK_j] = loglike_MPL_RISK(theta_hat_j,Y_j_risk,menuID_j_risk,algoList);
        
        % Other
        exit_flag(j,1) = exitflag_j_risk;
        fval(j,1)      = fval_j_risk;
        pdfCHECK(j,1)  = pdfCHECK_j;
        
    end
    
    
    
    %% Estimate time preferences
    
    % Create dataset with subset of time tasks only
    % Create dataset with subset of time tasks only
    menuID_j_time = menuID_j(timeTaskIdx_j);
    Y_j_time      = Y_j(timeTaskIdx_j);
    
    if ~isempty(menuID_j_time)
        
        auxY_j = Y_j_time(Y_j_time~=-1);
        
        % Check no switches
        if all(auxY_j==1)
            noSwitch_flag(j,2) = 1 ;
        end
        
        if all(auxY_j==0)
            noSwitch_flag(j,2) = -1 ;
        end
        
        % Define the objective function to minimize
        obj_fun = @(x) -loglike_MPL_TIME(x,Y_j_time,menuID_j_time,algoList);
        
        % Estimate        
        [theta_hat_j,fval_j_time,exitflag_j_time,~,~,hessian_loglike_j] = fminunc(obj_fun,X0_da,options_fminunc);
        
        % Estimated parameters
        alpha_da_j  = theta_hat_j(1);
        beta_da_j   = theta_hat_j(2);
        
        % Moments of the estimated distributions
        MOM_j = par2mom_da(theta_hat_j,algoList);
        mean_da_j   = MOM_j(1)  ;
        median_da_j = MOM_j(2)  ;
        mode_da_j   = MOM_j(3)  ;
        sd_da_j     = MOM_j(4)  ;
        Z_da_j      = MOM_j(5)  ;
        
        % Compute robust standard errors of estimated parameters
        cluster_var   = []; % Not used here so it doesn't matter
        Cov_theta_hat_j = RobustVarCov_MPL(@loglike_MPL_TIME, ...
            theta_hat_j,Y_da,menuID_da,-hessian_loglike_j,cluster_var,algoList,'simple');
        
        % Use delta method to find std. error. of moments of the estimated distributions
        gradient_MOM_hat_j = jacob_prec(@par2mom_da,theta_hat_j,1e-5,algoList);
        Cov_MOM_hat_j= gradient_MOM_hat_j*Cov_theta_hat_j*gradient_MOM_hat_j';
        
        % Std. dev. of estimators
        std_theta_hat_j = real( diag( sqrt(Cov_theta_hat_j) ) )  ;
        std_MOM_hat_j   = real( diag( sqrt(Cov_MOM_hat_j  ) ) )  ;
        
        % Parameters
        std_alpha_da_j = std_theta_hat_j(1);
        std_beta_da_j  = std_theta_hat_j(2);
        
        % Moments
        std_mean_da_j   = std_MOM_hat_j(1) ;
        std_median_da_j = std_MOM_hat_j(2) ;
        std_mode_da_j   = std_MOM_hat_j(3) ;
        std_sd_da_j     = std_MOM_hat_j(4) ;
        
        % Estimated parameters
        med_da_list(j) = median_da_j;
        sd_da_list(j)  = sd_da_j;
        
        % Standard deviations
        std_med_da_list(j) = std_median_da_j;
        std_sd_da_list(j)  = std_sd_da_j;
        
        % Check probabilities add to one
        [fval_j_time,~,~,pdfCHECK_j] = loglike_MPL_TIME(theta_hat_j,Y_j_time,menuID_j_time,algoList);
        
        % Other
        exit_flag(j,2) = exitflag_j_time;
        fval(j,2)      = fval_j_time;
        pdfCHECK(j,2)  = pdfCHECK_j;
        
    end
    
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
