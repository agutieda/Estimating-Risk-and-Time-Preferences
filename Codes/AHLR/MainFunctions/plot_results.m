% RANDOM MODELS FOR THE JOINT TREATMENT OF RISK AND TIME PREFERENCES
% by Jose Apesteguia, Miguel A. Ballester and Angelo Gutierrez
%
% Script to estimate DEU model using data from:
% "Eliciting Risk And Time Preferences"
% by Steffen Andersen, Glenn W. Harrison, Morten I. Lau, And E. Elisabet Rutstrom
% Econometrica, 2008
%
% This script:
% - Plot PDF and CDF of estimated distributions pooling all data
% - Plot histogram of CDF of estimated preferences by individual
% - Plot scatterplot of individual estimates
%
% November 2019
%
% Tested using Matlab 2019b

clear; close all; clc; rng(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PDF of estimated risk and time preferences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load results from pooled estimation
load('Output/results_separate.mat','mu_ra','sigma_ra','alpha_da','beta_da');
load('Output/results_separate.mat','min_ra','max_ra','min_da','max_da');

% Number of draws used to simulate this distribution
nDraws = 1000000;

% Generate draws from distribution of discount related to delay aversion
daDraws = betarnd(alpha_da,beta_da,nDraws,1);

% Generate draws from the distribution of risk aversion
p  = rand(nDraws,1);
Ga = normcdf(min_ra,mu_ra,sigma_ra);
Gb = normcdf(max_ra,mu_ra,sigma_ra);
u  = Ga + p.*(Gb-Ga);
raDraws = mu_ra-sigma_ra.*sqrt(2)*erfcinv(2*u);

% Generate draws of the discount factor
dfDraws = daDraws.^(1-raDraws);

% Trim extreme values to get a better estimate of the PDF and CDF around dr=1
% dfDraws_uptail = prctile(dfDraws,99);
dfDraws_uptail = 3;
dfDraws = dfDraws(dfDraws<dfDraws_uptail);

% Estimate pdf and cdf of df=da^(1-h)
if max_ra<1
    [CDF_df,df_axis] = ksdensity(dfDraws,'Support',[0,1],'Function','cdf');
    PDF_df = ksdensity(dfDraws,df_axis,'Support',[0,1],'Function','pdf');
else
    [CDF_df,df_axis] = ksdensity(dfDraws,'Support','positive','Function','cdf');
    PDF_df = ksdensity(dfDraws,df_axis,'Support','positive','Function','pdf');
end

% Support used to compute the plots
ra_axis = linspace(min_ra,max_ra,10000);
da_axis = linspace(min_da,max_da,10000);

% Axis Limits
min_ra_plot = min_ra;
max_ra_plot = max_ra;
min_da_plot = min_da;
max_da_plot = max_da;
min_df_plot = min_da;
max_df_plot = max_da;

% Plot options
face_alpha_plot = 0.4;
edge_alpha_plot = 0.8;

% PDF and CDF
PDF_ra = pdf_ra(ra_axis, mu_ra   , sigma_ra, min_ra, max_ra);
PDF_da = pdf_da(da_axis, alpha_da, beta_da , min_da, max_da);

fig = figure('Name','Estimated Distributions');
area(ra_axis,PDF_ra,'FaceAlpha',face_alpha_plot,'EdgeAlpha',edge_alpha_plot);
xlim([min_ra_plot,max_ra_plot]);
grid('on');
print(fig,'Output/EstimatedDistributions_AHLR_pdf_ra','-depsc2');

fig = figure('Name','Estimated Distributions');
area(da_axis,PDF_da,'FaceAlpha',face_alpha_plot,'EdgeAlpha',edge_alpha_plot);
xlim([min_da_plot,max_da_plot]);
grid('on');
print(fig,'Output/EstimatedDistributions_AHLR_pdf_da','-depsc2');

fig = figure('Name','Estimated Distributions');
area(df_axis,PDF_df,'FaceAlpha',face_alpha_plot,'EdgeAlpha',edge_alpha_plot);
xlim([min_df_plot,max_df_plot]);
grid('on');
print(fig,'Output/EstimatedDistributions_AHLR_pdf_df','-depsc2');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CDF of estimated risk and time preferences from pooled and individual estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load individual estimates
load('Output/results_individual.mat','med_ra_list','med_da_list','noSwitch_flag');

% Create lists used for further analysis
md_ra_plot = med_ra_list;
md_da_plot = med_da_list;

% df by subject
md_df_plot = md_da_plot.^(1-md_ra_plot);

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
md_ra_plot(idx_lo_ra) = lo_ra;
md_da_plot(idx_hi_da) = hi_da;
md_da_plot(idx_lo_da) = lo_da;

% Support used to compute the plots
ra_axis = linspace(min_ra,max_ra,10000);
da_axis = linspace(min_da,max_da,10000);

% Axis Limits
min_ra_plot = min_ra;
max_ra_plot = max_ra;
min_da_plot = min_da;
max_da_plot = max_da;
min_df_plot = min_da;
max_df_plot = max_da;

% CDF from pooled estimates
CDF_ra = cdf_ra(ra_axis,mu_ra,sigma_ra,min_ra,max_ra);
CDF_da = cdf_da(da_axis,alpha_da,beta_da,min_da,max_da);

% Plot options
face_alpha_plot_1 = 0.4;
edge_alpha_plot_1 = 0.8;
face_alpha_plot_2 = 0.5;
edge_alpha_plot_2 = 0.8;

% Number of bins used to compute histograms
nBinsHist = 30;

% Plot
fig = figure('Name','Estimated Distributions');
area(ra_axis,CDF_ra,'FaceAlpha',face_alpha_plot_1,'EdgeAlpha',edge_alpha_plot_1);
hold on;
% Adding "max_ra" in the following line to get the function plot the cdf
% plot the histogram of the cdf until the axis limit of the plot
histogram([md_ra_plot],'FaceAlpha',face_alpha_plot_2,'EdgeAlpha',edge_alpha_plot_2, ...
    'Normalization','cdf','NumBins',nBinsHist,'BinLimits',[min_ra_plot,max_ra_plot]);
hold off;
ylim([0 1]); xlim([min_ra_plot,max_ra_plot]);
grid('on');
print(fig,'Output/EstimatedDistributions_AHLR_cdf_ra','-depsc2');

fig = figure('Name','Estimated Distributions');
area(da_axis,CDF_da,'FaceAlpha',face_alpha_plot_1,'EdgeAlpha',edge_alpha_plot_1);
hold on;
histogram([md_da_plot],'FaceAlpha',face_alpha_plot_2,'EdgeAlpha',edge_alpha_plot_2,...
    'Normalization','cdf','NumBins',nBinsHist,'BinLimits',[min_da_plot,max_da_plot]);
hold off;
ylim([0 1]); xlim([min_da_plot,max_da_plot]);
grid('on');
print(fig,'Output/EstimatedDistributions_AHLR_cdf_da','-depsc2');

fig = figure('Name','Estimated Distributions');
area(df_axis,CDF_df,'FaceAlpha',face_alpha_plot_1,'EdgeAlpha',edge_alpha_plot_1);
hold on;
histogram([md_df_plot],'FaceAlpha',face_alpha_plot_2,'EdgeAlpha',edge_alpha_plot_2,...
    'Normalization','cdf','NumBins',nBinsHist,'BinLimits',[min_df_plot,max_df_plot]);
hold off;
ylim([0 1]); xlim([min_df_plot,max_df_plot]);
grid('on');
print(fig,'Output/EstimatedDistributions_AHLR_cdf_df','-depsc2');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Scatterplot of estimates by individual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure1 = figure('Name','Risk Aversion vs Delay Aversion by Subject');
set(figure1,'Units','Normalized','Position', [0.1 0.1 0.8 0.8]);
axes1 = axes('Parent',figure1,'XGrid','on','YGrid','on'); hold(axes1,'on');
ylabel('$h$','FontSize',16,'Interpreter','latex');
xlabel('$\hat{\delta}$','FontSize',16,'Interpreter','latex');
scatter(md_da_plot,md_ra_plot,...
    'MarkerFaceColor',[0.301960784313725 0.749019607843137 0.929411764705882],...
    'MarkerEdgeColor',[0 0.450980392156863 0.741176470588235],...
    'LineWidth',1.5);
regLine = lsline;
set(regLine,'Parent',axes1,'LineWidth',1.2,'LineStyle','--',...
    'Color',[0.87058824300766 0.490196079015732 0]);
xlim([0.9,1]); ylim([-1,1]);
print(figure1,'Output/IndividualScatter_AHLR','-depsc2');


