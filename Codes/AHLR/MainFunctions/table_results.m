% RANDOM MODELS FOR THE JOINT TREATMENT OF RISK AND TIME PREFERENCES
% by Jose Apesteguia, Miguel A. Ballester and Angelo Gutierrez
%
% Script to estimate DEU model using data from:
% "Eliciting Risk And Time Preferences"
% by Steffen Andersen, Glenn W. Harrison, Morten I. Lau, And E. Elisabet Rutstrom
% Econometrica, 2008
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

% Store in two columns
results_sep_ra = [
    median_ra;std_median_ra;
    sd_ra;std_sd_ra;
    nan;nan;
    nan;nan;
    nObs_ra;log_like_hat_ra;
    ];

results_sep_da = [
    nan;nan;
    nan;nan;
    median_da;std_median_da;
    sd_da;std_sd_da;
    nObs_da;log_like_hat_da;
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
    nObs_all;log_like_hat;
    ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimates from compute_individual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Load results from compute_individual
load('Output/results_individual.mat','med_ra_list','sd_ra_list','med_da_list','sd_da_list');
load('Output/results_individual.mat','fval','noSwitch_flag','nIndividuals');

% Create lists used for further analysis
md_ra_plot = med_ra_list;
xd_ra_plot = sd_ra_list;
md_da_plot = med_da_list;
xd_da_plot = sd_da_list;

% % Identify individuals that never switch in the tasks
idx_hi_ra = noSwitch_flag(:,1) ==  1; % Individuals always choose alternative 1 (safe lottery) in all risky tasks -> Risk averse (High ra)
idx_lo_ra = noSwitch_flag(:,1) == -1; % Individuals always choose alternative 2 (risky lottery)  in all risk tasks -> Risk loving (Low ra)
idx_hi_da = noSwitch_flag(:,2) ==  1; % Individuals always choose alternative 1 (sooner lottery) in all time tasks -> High delay aversion (Low delta_hat)
idx_lo_da = noSwitch_flag(:,2) == -1; % Individuals always choose alternative 2 (later lottery)  in all time tasks -> Low delay aversion (High delta_hata)

% % Imput a level of risk aversion and delay aversion to these individuals
hi_ra =  1            ; % median ra level assigned to subjects that always choose safest lottery
lo_ra = -2            ; % median ra level assigned to subjects that always choose risky lottery
hi_da = 0.60^(1/12)   ; % median da level assigned to subjects that always choose sooner lottery
lo_da = 0.96^(1/12)   ; % median da level assigned to subjects that always choose later lottery

md_ra_plot(idx_hi_ra) = hi_ra;
xd_ra_plot(idx_hi_ra) = 0.01;
md_ra_plot(idx_lo_ra) = lo_ra;
xd_ra_plot(idx_lo_ra) = 0.01;
md_da_plot(idx_hi_da) = hi_da;
xd_da_plot(idx_hi_da) = 0.01;
md_da_plot(idx_lo_da) = lo_da;
xd_da_plot(idx_lo_da) = 0.01;

% Compute summary statistics of the individual estimates
med_md_ra_bySubject = median(md_ra_plot);
std_md_ra_bySubject = std(md_ra_plot);
med_md_da_bySubject = median(md_da_plot);
std_md_da_bySubject = std(md_da_plot);

med_xd_ra_bySubject = median(xd_ra_plot);
std_xd_ra_bySubject = std(xd_ra_plot);
med_xd_da_bySubject = median(xd_da_plot);
std_xd_da_bySubject = std(xd_da_plot);

rho_rada_bySubject   = corr(md_ra_plot,md_da_plot);
rhoXD_rada_bySubject = corr(xd_ra_plot,xd_da_plot);
log_like_bySubject   = mean(fval(:,1))*0.4+mean(fval(:,2))*0.6;

% Store in a column
results_individual = [
    med_md_ra_bySubject;std_md_ra_bySubject;
    med_xd_ra_bySubject;std_xd_ra_bySubject;
    med_md_da_bySubject;std_md_da_bySubject;
    med_xd_da_bySubject;std_xd_da_bySubject;
    nIndividuals;log_like_bySubject;
    ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DISPLAY RESULTS IN A TABLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Round to three decimal places
results_sep_ra = round(results_sep_ra,3);
results_sep_da = round(results_sep_da,3);
results_joint  = round(results_joint,3);
results_individual = round(results_individual,3);

disp(' ')
disp(' ')
disp('************* Estimated Risk and Time Preferences: Andersen et al. (2008) ******************')
Col_Names = {'Risk Only', 'Time Only','Joint','by Individual'};
Row_Names = {'med_ra','se_med_ra','sd_ra','se_sd_ra','med_da','se_med_da','sd_da','se_sd_da','nObs','loglike'};
table(results_sep_ra,results_sep_da,results_joint,results_individual,...
    'RowNames',Row_Names,'VariableNames',Col_Names)
disp(' ')
disp('********************************************************************************************')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXPORT TO LATEX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

table_AHLR = table(results_sep_ra,results_sep_da,results_joint,results_individual,...
    'RowNames',Row_Names,'VariableNames',Col_Names);

table2latex(table_AHLR,'Output/table_AHLR.tex');

writetable(table_AHLR,'Output/table_AHLR.csv','WriteRowNames',1,'WriteVariableNames',1);
