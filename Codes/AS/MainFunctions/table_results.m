% RANDOM MODELS FOR THE JOINT TREATMENT OF RISK AND TIME PREFERENCES
% by Jose Apesteguia, Miguel A. Ballester and Angelo Gutierrez
%
% Script to estimate DEU-H, PVCE-H and CEPV-H models using data from:
% "Estimating Time Preferences from Convex Budgets"
%  by James Andreoni and Charles Sprenger.
%  American Economic Review, (2012)
%
% This file:
% - Summarize estimates in a table and export to latex and csv
%
% November 2019
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
    nan(5*2,1);
    log_like_hat;
    ];

%%%% PVCE %%%%

% Load results
load('Output/results_baseline_PVCE.mat','median_ra','sd_ra','median_da','sd_da','median_is','sd_is');
load('Output/results_baseline_PVCE.mat','std_median_ra','std_sd_ra',...
    'std_median_da','std_sd_da','std_median_is','std_sd_is');
load('Output/results_baseline_PVCE.mat','log_like_hat');

% Store in column
results_baseline_PVCE = [
    median_ra;std_median_ra;
    sd_ra;std_sd_ra;
    median_da;std_median_da;
    sd_da;std_sd_da;
    median_is;std_median_is;
    sd_is;std_sd_is;
    nan(3*2,1);
    log_like_hat;
    ];

%%%% CEPV %%%%

% Load results
load('Output/results_baseline_CEPV.mat','median_ra','sd_ra','median_da','sd_da','median_is','sd_is');
load('Output/results_baseline_CEPV.mat','std_median_ra','std_sd_ra',...
    'std_median_da','std_sd_da','std_median_is','std_sd_is');
load('Output/results_baseline_CEPV.mat','log_like_hat');

% Store in column
results_baseline_CEPV = [
    median_ra;std_median_ra;
    sd_ra;std_sd_ra;
    median_da;std_median_da;
    sd_da;std_sd_da;
    median_is;std_median_is;
    sd_is;std_sd_is;
    nan(3*2,1);
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
    nan(2*2,1);
    rho_rada;std_rho_rada;
    nan(2*2,1);
    log_like_hat;
    ];

%%%% PVCE %%%%

% Load results
load('Output/results_correlated_PVCE.mat','median_ra','sd_ra','median_da','sd_da','median_is','sd_is');
load('Output/results_correlated_PVCE.mat','std_median_ra','std_sd_ra',...
    'std_median_da','std_sd_da','std_median_is','std_sd_is');
load('Output/results_correlated_PVCE.mat','rho_rada','rho_rais','rho_dais');
load('Output/results_correlated_PVCE.mat','std_rho_rada','std_rho_rais','std_rho_dais');
load('Output/results_correlated_PVCE.mat','log_like_hat');

% Store in column
results_correlated_PVCE = [
    median_ra;std_median_ra;
    sd_ra;std_sd_ra;
    median_da;std_median_da;
    sd_da;std_sd_da;
    median_is;std_median_is;
    sd_is;std_sd_is;
    rho_rada;std_rho_rada;
    rho_rais;std_rho_rais;
    rho_dais;std_rho_dais;
    log_like_hat;
    ];

%%%% CEPV %%%%

% Load results
load('Output/results_correlated_CEPV.mat','median_ra','sd_ra','median_da','sd_da','median_is','sd_is');
load('Output/results_correlated_CEPV.mat','std_median_ra','std_sd_ra',...
    'std_median_da','std_sd_da','std_median_is','std_sd_is');
load('Output/results_correlated_CEPV.mat','rho_rada','rho_rais','rho_dais');
load('Output/results_correlated_CEPV.mat','std_rho_rada','std_rho_rais','std_rho_dais');
load('Output/results_correlated_CEPV.mat','log_like_hat');

% Store in column
results_correlated_CEPV = [
    median_ra;std_median_ra;
    sd_ra;std_sd_ra;
    median_da;std_median_da;
    sd_da;std_sd_da;
    median_is;std_median_is;
    sd_is;std_sd_is;
    rho_rada;std_rho_rada;
    rho_rais;std_rho_rais;
    rho_dais;std_rho_dais;
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
    nan(2*2,1);
    rho_rada_bySubject;nan;
    nan(2*2,1);
    loglike_bySubject;
    ];

%%%% PVCE %%%%

%%%% Load results from compute_individual
load('Output/results_individual_PVCE.mat','med_ra_list','sd_ra_list',...
    'med_da_list','sd_da_list','med_is_list','sd_is_list');
load('Output/results_individual_PVCE.mat','fval');

% Compute summary statistics of the individual estimates
med_md_ra_bySubject = median(med_ra_list);
sd_md_ra_bySubject  = std(med_ra_list);
med_md_da_bySubject = median(med_da_list);
sd_md_da_bySubject  = std(med_da_list);

med_xd_ra_bySubject = median(sd_ra_list);
sd_xd_ra_bySubject  = std(sd_ra_list);
med_xd_da_bySubject = median(sd_da_list);
sd_xd_da_bySubject  = std(sd_da_list);

med_md_is_bySubject = median(med_is_list);
sd_md_is_bySubject  = std(med_is_list);
med_xd_is_bySubject = median(sd_is_list);
sd_xd_is_bySubject  = std(sd_is_list);

corrMat = corr([med_ra_list,med_da_list,med_is_list]);
rho_rada_bySubject  = corrMat(1,2);
rho_rais_bySubject  = corrMat(1,3);
rho_dais_bySubject  = corrMat(2,3);

loglike_bySubject   = mean(fval);

% Store in column
results_individual_PVCE = [
    med_md_ra_bySubject;sd_md_ra_bySubject;
    med_xd_ra_bySubject;sd_xd_ra_bySubject;
    med_md_da_bySubject;sd_md_da_bySubject;
    med_xd_da_bySubject;sd_xd_da_bySubject;
    med_md_is_bySubject;sd_md_is_bySubject;
    med_xd_is_bySubject;sd_xd_is_bySubject;
    rho_rada_bySubject;nan;
    rho_rais_bySubject;nan;
    rho_dais_bySubject;nan;
    loglike_bySubject;
    ];

%%%% CEPV %%%%

%%%% Load results from compute_individual
load('Output/results_individual_CEPV.mat','med_ra_list','sd_ra_list',...
    'med_da_list','sd_da_list','med_is_list','sd_is_list');
load('Output/results_individual_CEPV.mat','fval');

% Compute summary statistics of the individual estimates
med_md_ra_bySubject = median(med_ra_list);
sd_md_ra_bySubject  = std(med_ra_list);
med_md_da_bySubject = median(med_da_list);
sd_md_da_bySubject  = std(med_da_list);

med_xd_ra_bySubject = median(sd_ra_list);
sd_xd_ra_bySubject  = std(sd_ra_list);
med_xd_da_bySubject = median(sd_da_list);
sd_xd_da_bySubject  = std(sd_da_list);

med_md_is_bySubject = median(med_is_list);
sd_md_is_bySubject  = std(med_is_list);
med_xd_is_bySubject = median(sd_is_list);
sd_xd_is_bySubject  = std(sd_is_list);

corrMat = corr([med_ra_list,med_da_list,med_is_list]);
rho_rada_bySubject  = corrMat(1,2);
rho_rais_bySubject  = corrMat(1,3);
rho_dais_bySubject  = corrMat(2,3);

loglike_bySubject   = mean(fval);

% Store in column
results_individual_CEPV = [
    med_md_ra_bySubject;sd_md_ra_bySubject;
    med_xd_ra_bySubject;sd_xd_ra_bySubject;
    med_md_da_bySubject;sd_md_da_bySubject;
    med_xd_da_bySubject;sd_xd_da_bySubject;
    med_md_is_bySubject;sd_md_is_bySubject;
    med_xd_is_bySubject;sd_xd_is_bySubject;
    rho_rada_bySubject;nan;
    rho_rais_bySubject;nan;
    rho_dais_bySubject;nan;
    loglike_bySubject;
    ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DISPLAY RESULTS IN A TABLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Round to three decimal places
results_baseline_DEU  = round(results_baseline_DEU ,3);
results_baseline_PVCE = round(results_baseline_PVCE,3);
results_baseline_CEPV = round(results_baseline_CEPV,3);

results_correlated_DEU  = round(results_correlated_DEU ,3);
results_correlated_PVCE = round(results_correlated_PVCE,3);
results_correlated_CEPV = round(results_correlated_CEPV,3);

results_individual_DEU  = round(results_individual_DEU ,3);
results_individual_PVCE = round(results_individual_PVCE,3);
results_individual_CEPV = round(results_individual_CEPV,3);

disp(' ')
disp(' ')
disp('********************************************************** Estimated Risk and Time Preferences: Andreoni and Sprenger (2012) ***********************************************************')
Col_Names = {...
    'Baseline-DEU', 'Correlated-DEU','byIndividual-DEU',...
    'Baseline-PVCE', 'Correlated-PVCE','byIndividual-PVCE',...
    'Baseline-CEPV', 'Correlated-CEPV','byIndividual-CEPV'...
    };

Row_Names = {...
    'med_ra','se_med_ra','sd_ra','se_sd_ra',...
    'med_da','se_med_da','sd_da','se_sd_da',...
    'med_is','se_med_is','sd_is','se_sd_is',...
    'rho_ra_da','se_rho_ra_da',...
    'rho_ra_is','se_rho_ra_is',...
    'rho_da_is','se_rho_da_is',...
    'loglike'
    };

table(results_baseline_DEU,results_correlated_DEU,results_individual_DEU,...
    results_baseline_PVCE,results_correlated_PVCE,results_individual_PVCE,...
    results_baseline_CEPV,results_correlated_CEPV,results_individual_CEPV,...
    'RowNames',Row_Names,'VariableNames',Col_Names)

disp(' ')
disp('***************************************************************************************************************************************************************************************')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXPORT TO LATEX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

table_AS = table(results_baseline_DEU,results_correlated_DEU,results_individual_DEU,...
    results_baseline_PVCE,results_correlated_PVCE,results_individual_PVCE,...
    results_baseline_CEPV,results_correlated_CEPV,results_individual_CEPV,...
    'RowNames',Row_Names,'VariableNames',Col_Names);

table2latex(table_AS,'Output/table_AS.tex');

writetable(table_AS,'Output/table_AS.csv','WriteRowNames',1,'WriteVariableNames',1);


