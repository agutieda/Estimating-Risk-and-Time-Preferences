% RANDOM MODELS FOR THE JOINT TREATMENT OF RISK AND TIME PREFERENCES
% by Jose Apesteguia, Miguel A. Ballester and Angelo Gutierrez
%
% Script to estimate DEU-H model using data from:
% "Estimating Time Preferences from Convex Budgets"
%  by James Andreoni and Charles Sprenger.
%  American Economic Review, (2012)
%
% This file:
% - Summarize estimates in a table and export to latex and csv
%
% March 2020
%
% Tested using Matlab 2019b

clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimates from compute_baseline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% DEU %%%%

% Load results
load('Output/results_baseline_DEU.mat','median_ra','sd_ra','median_da','sd_da');
load('Output/results_baseline_DEU.mat','std_median_ra','std_sd_ra','std_median_da','std_sd_da');
load('Output/results_baseline_DEU.mat','log_like_hat');

% Store in column
results_baseline_DEU = [
    median_ra;std_median_ra;
    sd_ra;std_sd_ra;
    median_da;std_median_da;
    sd_da;std_sd_da;
    nan(2,1);
    log_like_hat;
    ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimates from compute_correlated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% DEU %%%%

% Load results
load('Output/results_correlated_DEU.mat','median_ra','sd_ra','median_da','sd_da');
load('Output/results_correlated_DEU.mat','std_median_ra','std_sd_ra','std_median_da','std_sd_da');
load('Output/results_correlated_DEU.mat','rho_rada','std_rho_rada');
load('Output/results_correlated_DEU.mat','log_like_hat');

% Store in column
results_correlated_DEU = [
    median_ra;std_median_ra;
    sd_ra;std_sd_ra;
    median_da;std_median_da;
    sd_da;std_sd_da;
    rho_rada;std_rho_rada;    
    log_like_hat;
    ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimates from compute_individual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% DEU %%%%

%%%% Load results from compute_individual
load('Output/results_individual_DEU.mat','med_ra_list','sd_ra_list','med_da_list','sd_da_list');
load('Output/results_individual_DEU.mat','fval');

% Compute summary statistics of the individual estimates
med_md_ra_bySubject = median(med_ra_list);
sd_md_ra_bySubject  = std(med_ra_list);
med_md_da_bySubject = median(med_da_list);
sd_md_da_bySubject  = std(med_da_list);

med_xd_ra_bySubject = median(sd_ra_list);
sd_xd_ra_bySubject  = std(sd_ra_list);
med_xd_da_bySubject = median(sd_da_list);
sd_xd_da_bySubject  = std(sd_da_list);

rho_rada_bySubject  = corr(med_ra_list,med_da_list);
loglike_bySubject   = mean(fval(:,1));

% Store in column
results_individual_DEU = [
    med_md_ra_bySubject;sd_md_ra_bySubject;
    med_xd_ra_bySubject;sd_xd_ra_bySubject;
    med_md_da_bySubject;sd_md_da_bySubject;
    med_xd_da_bySubject;sd_xd_da_bySubject;    
    rho_rada_bySubject;nan;    
    loglike_bySubject;
    ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DISPLAY RESULTS IN A TABLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Round to three decimal places
results_baseline_DEU  = round(results_baseline_DEU ,3);
results_correlated_DEU  = round(results_correlated_DEU ,3);
results_individual_DEU  = round(results_individual_DEU ,3);

disp(' ')
disp(' ')
disp('********************************************************** Estimated Risk and Time Preferences: Andreoni and Sprenger (2012) ***********************************************************')
Col_Names = {'All-Tasks', 'All-Tasks-Correlated','by-Individual'};
Row_Names = {...
    'med_ra','se_med_ra','sd_ra','se_sd_ra',...
    'med_da','se_med_da','sd_da','se_sd_da',...
    'rho_ra_da','se_rho_ra_da',...
    'loglike'
    };

table(results_baseline_DEU,results_correlated_DEU,results_individual_DEU,...
    'RowNames',Row_Names,'VariableNames',Col_Names)

disp(' ')
disp('***************************************************************************************************************************************************************************************')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXPORT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

table_AS = table(results_baseline_DEU,results_correlated_DEU,results_individual_DEU,...
    'RowNames',Row_Names,'VariableNames',Col_Names);
table2latex(table_AS,'Output/table_AS.tex');

writetable(table_AS,'Output/table_AS.csv','WriteRowNames',1,'WriteVariableNames',1);
