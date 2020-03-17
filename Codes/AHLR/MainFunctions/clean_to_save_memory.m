% RANDOM MODELS FOR THE JOINT TREATMENT OF RISK AND TIME PREFERENCES
% by Jose Apesteguia, Miguel A. Ballester and Angelo Gutierrez
%
% Script to estimate DEU model using data from:
% "Eliciting Risk And Time Preferences"
% by Steffen Andersen, Glenn W. Harrison, Morten I. Lau, And E. Elisabet Rutstrom
% Econometrica, 2008
%
% This file:
% - Save memory by keeping only the variables used in
%   table_results and plot_results
%
% March 2020
%
% Tested using Matlab 2019b

clear;
load('Output/results_joint.mat');
clearvars -except mu_ra sigma_ra alpha_da beta_da ...
    median_ra sd_ra median_da sd_da ...
    std_median_ra std_sd_ra std_median_da std_sd_da ...
    min_ra max_ra min_da max_da ...
    log_like_hat nObs_all;
save('Output/results_joint.mat');

clear;
load('Output/results_separate.mat');
clearvars -except mu_ra sigma_ra alpha_da beta_da ...
    median_ra sd_ra median_da sd_da ...
    std_median_ra std_sd_ra std_median_da std_sd_da ...
    min_ra max_ra min_da max_da ...
    log_like_hat_ra log_like_hat_da nObs_ra nObs_da;
save('Output/results_separate.mat');

clear;
load('Output/results_individual.mat');
clearvars -except med_ra_list sd_ra_list med_da_list sd_da_list ...
    min_ra max_ra min_da max_da ...
    nIndividuals fval noSwitch_flag;
save('Output/results_individual.mat');
