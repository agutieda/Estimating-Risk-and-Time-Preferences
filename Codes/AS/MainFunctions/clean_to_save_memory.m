% RANDOM MODELS FOR THE JOINT TREATMENT OF RISK AND TIME PREFERENCES
% by Jose Apesteguia, Miguel A. Ballester and Angelo Gutierrez
%
% Script to estimate DEU-H, PVCE-H and CEPV-H models using data from:
% "Estimating Time Preferences from Convex Budgets"
%  by James Andreoni and Charles Sprenger.
%  American Economic Review, (2012)
%
% This file:
% - Save memory by keeping only the variables used in
%   table_results and plot_results
%
% November 2019
%
% Tested using Matlab 2019b


%%% DEU %%%

clear;
load('Output/results_baseline_DEU.mat');
clearvars -except mu_ra sigma_ra alpha_da beta_da ...
    median_ra sd_ra median_da sd_da ...
    std_median_ra std_sd_ra std_median_da std_sd_da ...
    taskSetID Y menuID subjectID rhoX nX xList log_like_hat;
save('Output/results_baseline_DEU.mat');

clear;
load('Output/results_correlated_DEU.mat');
clearvars -except mu_ra sigma_ra alpha_da beta_da ...
    rho_rada std_rho_rada ...
    median_ra sd_ra median_da sd_da ...
    std_median_ra std_sd_ra std_median_da std_sd_da ...
    taskSetID Y menuID subjectID rhoX nX xList log_like_hat;
save('Output/results_correlated_DEU.mat');

clear;
load('Output/results_individual_DEU.mat');
clearvars -except med_ra_list sd_ra_list med_da_list sd_da_list fval;
save('Output/results_individual_DEU.mat');


%%% PVCE %%%

clear;
load('Output/results_baseline_PVCE.mat');
clearvars -except mu_ra sigma_ra alpha_da beta_da mu_is sigma_is ...
    median_ra sd_ra median_da sd_da median_is sd_is ...
    std_median_ra std_sd_ra std_median_da std_sd_da std_median_is std_sd_is ...
    taskSetID Y menuID subjectID rhoX nX xList log_like_hat;
save('Output/results_baseline_PVCE.mat');

clear;
load('Output/results_correlated_PVCE.mat');
clearvars -except mu_ra sigma_ra alpha_da beta_da mu_is sigma_is ...
    median_ra sd_ra median_da sd_da median_is sd_is ...
    rho_rada rho_rais rho_dais std_rho_rada std_rho_rais std_rho_dais ...
    std_median_ra std_sd_ra std_median_da std_sd_da std_median_is std_sd_is ...
    taskSetID Y menuID subjectID rhoX nX xList log_like_hat;
save('Output/results_correlated_PVCE.mat');

clear;
load('Output/results_individual_PVCE.mat');
clearvars -except med_ra_list sd_ra_list med_da_list sd_da_list ...
    med_is_list sd_is_list fval;
save('Output/results_individual_PVCE.mat');


%%% CEPV %%%

load('Output/results_baseline_CEPV.mat');
clearvars -except mu_ra sigma_ra alpha_da beta_da mu_is sigma_is ...
    median_ra sd_ra median_da sd_da median_is sd_is ...
    std_median_ra std_sd_ra std_median_da std_sd_da std_median_is std_sd_is ...
    taskSetID Y menuID subjectID rhoX nX xList log_like_hat;
save('Output/results_baseline_CEPV.mat');
clear;

clear;
load('Output/results_correlated_CEPV.mat');
clearvars -except mu_ra sigma_ra alpha_da beta_da mu_is sigma_is ...
    median_ra sd_ra median_da sd_da median_is sd_is ...
    rho_rada rho_rais rho_dais std_rho_rada std_rho_rais std_rho_dais ...
    std_median_ra std_sd_ra std_median_da std_sd_da std_median_is std_sd_is ...
    taskSetID Y menuID subjectID rhoX nX xList log_like_hat;
save('Output/results_correlated_CEPV.mat');

clear;
load('Output/results_individual_CEPV.mat');
clearvars -except med_ra_list sd_ra_list med_da_list sd_da_list ...
    med_is_list sd_is_list fval;
save('Output/results_individual_CEPV.mat');

