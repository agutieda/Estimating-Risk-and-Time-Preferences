% RANDOM MODELS FOR THE JOINT TREATMENT OF RISK AND TIME PREFERENCES
% by Jose Apesteguia, Miguel A. Ballester and Angelo Gutierrez
%
% Script to estimate DEU model using data from:
% "At the nexus of risk and time preferences: An experimental investigation" (2010)
% by Keith H. Coble & Jayson L. Lusk
% Journal of Risk and Uncertainty
%
% This script:
% - Summarize estimates in a table and export to latex
%
% March 2020
%
% Tested using Matlab 2019b

clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimates from compute_separate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load results
load('Output/results_separate.mat','median_ra','sd_ra','median_da','sd_da');
load('Output/results_separate.mat','std_median_ra','std_sd_ra','std_median_da','std_sd_da');
load('Output/results_separate.mat','log_like_hat_ra','log_like_hat_da','nObs_ra','nObs_da');
load('Output/results_separate.mat','median_ra_rada','sd_ra_rada','median_da_rada','sd_da_rada');
load('Output/results_separate.mat','std_median_ra_rada','std_sd_ra_rada','std_median_da_rada','std_sd_da_rada');
load('Output/results_separate.mat','log_like_hat_rada','nObs_rada');

% Store in two columns
results_sep_ra = [
    median_ra;std_median_ra;
    sd_ra;std_sd_ra;
    nan;nan;
    nan;nan;
    nan;nan;
    nObs_ra;log_like_hat_ra;
    ];

results_sep_da = [
    nan;nan;
    nan;nan;
    median_da;std_median_da;
    sd_da;std_sd_da;
    nan;nan;
    nObs_da;log_like_hat_da;
    ];

results_sep_rada = [
    median_ra_rada;std_median_ra_rada;
    sd_ra_rada;std_sd_ra_rada;
    median_da_rada;std_median_da_rada;
    sd_da_rada;std_sd_da_rada;
    nan;nan;
    nObs_rada;log_like_hat_rada;
    ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimates from compute_joint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load results
load('Output/results_joint.mat','median_ra','sd_ra','median_da','sd_da');
load('Output/results_joint.mat','std_median_ra','std_sd_ra','std_median_da','std_sd_da');
load('Output/results_joint.mat','log_like_hat','nObs_all');

% Store in a column
results_joint = [
    median_ra;std_median_ra;
    sd_ra;std_sd_ra;
    median_da;std_median_da;
    sd_da;std_sd_da;
    nan;nan;
    nObs_all;log_like_hat;
    ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimates from compute_correlated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load results
load('Output/results_correlated.mat','median_ra','sd_ra','median_da','sd_da');
load('Output/results_correlated.mat','std_median_ra','std_sd_ra','std_median_da','std_sd_da');
load('Output/results_correlated.mat','rho_rada','std_rho_rada');
load('Output/results_correlated.mat','log_like_hat','nObs_all');

% Store in a column
results_correlated = [
    median_ra;std_median_ra;
    sd_ra;std_sd_ra;
    median_da;std_median_da;
    sd_da;std_sd_da;
    rho_rada;std_rho_rada;
    nObs_all;log_like_hat;
    ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimates from compute_individual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Load results from compute_individual
load('Output/results_individual.mat','med_ra_list','sd_ra_list','med_da_list','sd_da_list');
load('Output/results_individual.mat','fval','noSwitch_flag','nIndividuals');

med_md_ra_bySubject = median(med_ra_list);
sd_md_ra_bySubject  = std(med_ra_list);
med_md_da_bySubject = median(med_da_list);
sd_md_da_bySubject  = std(med_da_list);

med_xd_ra_bySubject = median(sd_ra_list);
sd_xd_ra_bySubject  = std(sd_ra_list);
med_xd_da_bySubject = median(sd_da_list);
sd_xd_da_bySubject  = std(sd_da_list);

rho_rada_bySubject   = corr(med_ra_list,med_da_list);
log_like_bySubject   = mean(fval);

% Store in a column
results_individual = [
    med_md_ra_bySubject;sd_md_ra_bySubject;
    med_xd_ra_bySubject;sd_xd_ra_bySubject;
    med_md_da_bySubject;sd_md_da_bySubject;
    med_xd_da_bySubject;sd_xd_da_bySubject;
    rho_rada_bySubject;nan;
    nIndividuals;log_like_bySubject(1);
    ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DISPLAY RESULTS IN A TABLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Round to three decimal places
results_sep_ra     = round(results_sep_ra,3);
results_sep_da     = round(results_sep_da,3);
results_sep_rada   = round(results_sep_rada,3);
results_joint      = round(results_joint,3);
results_correlated = round(results_correlated,3);
results_individual = round(results_individual,3);

disp(' ')
disp(' ')
disp('************* Estimated Risk and Time Preferences: Coble and Lusk (2010) ******************')
Col_Names = {'Risk Only', 'Time Only','Risk+Time Only','Joint','Correlated','by Individual'};
Row_Names = {'med_ra','se_med_ra','sd_ra','se_sd_ra','med_da','se_med_da',...
    'sd_da','se_sd_da','rho_rada','se_rho_rada','nObs','loglike'};
table(results_sep_ra,results_sep_da,results_sep_rada,results_joint,...
    results_correlated,results_individual,...
    'RowNames',Row_Names,'VariableNames',Col_Names)
disp(' ')
disp('********************************************************************************************')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXPORT TO LATEX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

table_CL = table(results_sep_ra,results_sep_da,results_sep_rada,results_joint,...
    results_correlated,results_individual,...
    'RowNames',Row_Names,'VariableNames',Col_Names);

table2latex(table_CL,'Output/table_CL.tex');

writetable(table_CL,'Output/table_CL.csv','WriteRowNames',1,'WriteVariableNames',1);
