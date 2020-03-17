% RANDOM MODELS FOR THE JOINT TREATMENT OF RISK AND TIME PREFERENCES
% by Jose Apesteguia, Miguel A. Ballester and Angelo Gutierrez
%
% Script to estimate DEU-H model using data from:
% "Estimating Time Preferences from Convex Budgets"
%  by James Andreoni and Charles Sprenger.
%  American Economic Review, (2012)
%
% This file:
% - Plot PDF of estimated preferences under DEU-H by pooling all data
% - Plot observed and predicted distributions of choices
% - Plot scatterplot of individual estimates
%
% March 2020
%
% Tested using Matlab 2019b

clear; close all; clc; rng(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PDF of estimated preferences under each representation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% LOAD DATA

% DEU
load('Output/results_baseline_DEU.mat','mu_ra','sigma_ra','alpha_da','beta_da');
mu_ra_DEU    =  mu_ra    ;
sigma_ra_DEU =  sigma_ra ;
alpha_da_DEU =  alpha_da ;
beta_da_DEU  =  beta_da  ;


% Specify some bounds for the plots
min_ra = -10      ; % Lower limit of the integral on ra
max_ra =  1       ; % Upper limit of the integral on ra
min_da =  0.5     ; % Lower limit of the integral on da
max_da =  1       ; % Upper limit of the integral on da


%%%%% ESTIMATE DISTRIBUTION OF DISCOUNT FACTOR USING ESTIMATED DISTRIBUTION OF DELAY AVERSION

% Number of draws used to simulate this distribution
nDraws = 50000000;

%%% DEU

% Generate draws from distribution of discount related to delay aversion
daDraws = betarnd(alpha_da_DEU,beta_da_DEU,nDraws,1);

% Generate draws from the distribution of risk aversion
p  = rand(nDraws,1);
Ga = normcdf(min_ra,mu_ra_DEU,sigma_ra_DEU);
Gb = normcdf(max_ra,mu_ra_DEU,sigma_ra_DEU);
u  = Ga + p.*(Gb-Ga);
raDraws = mu_ra_DEU-sigma_ra_DEU.*sqrt(2)*erfcinv(2*u);

% Generate draws of the discount factor
dfDraws = daDraws.^(1-raDraws);

% Trim extreme values to get a better estimate of the PDF and CDF around dr=1
% dfDraws_uptail = prctile(dfDraws,99);
dfDraws_uptail = 3;
dfDraws = dfDraws(dfDraws<dfDraws_uptail);

% Estimate pdf and cdf
[CDF_df_DEU,df_axis_DEU] = ksdensity(dfDraws,'Support','positive','Function','cdf');
PDF_df_DEU = ksdensity(dfDraws,df_axis_DEU,'Support','positive','Function','pdf');

clearvars daDraws p Ga Gb u raDraws dfDraws;

%%% PLOT PDFs for DEU

% Axis Limits
min_ra_plot = -2;
max_ra_plot = 1;
min_da_plot = 0.70;
max_da_plot = 1;
min_df_plot = 0.70;
max_df_plot = 1;

% Support used to compute the plots
ra_axis = linspace(min_ra,max_ra,10000);
da_axis = linspace(min_da,max_da,10000);

% Plot options
face_alpha_plot = 0.4;
edge_alpha_plot = 0.8;

% PDF and CDF
PDF_ra = pdf_ra(ra_axis, mu_ra_DEU   , sigma_ra_DEU, min_ra, max_ra);
PDF_da = pdf_da(da_axis, alpha_da_DEU, beta_da_DEU , min_da, max_da);

fig = figure('Name','Estimated Distributions');
area(ra_axis,PDF_ra,'FaceAlpha',face_alpha_plot,'EdgeAlpha',edge_alpha_plot);
xlim([min_ra_plot,max_ra_plot]);
grid('on');
print(fig,'Output/EstimatedDistributions_AS_pdf_ra','-depsc2');

fig = figure('Name','Estimated Distributions');
area(da_axis,PDF_da,'FaceAlpha',face_alpha_plot,'EdgeAlpha',edge_alpha_plot);
xlim([min_da_plot,max_da_plot]);
grid('on');
print(fig,'Output/EstimatedDistributions_AS_pdf_da','-depsc2');

fig = figure('Name','Estimated Distributions');
area(df_axis_DEU,PDF_df_DEU,'FaceAlpha',face_alpha_plot,'EdgeAlpha',edge_alpha_plot);
xlim([min_df_plot,max_df_plot]);
grid('on');
print(fig,'Output/EstimatedDistributions_AS_pdf_df','-depsc2');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CDF of estimated risk and time preferences from pooled and individual estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load individual estimates
load('Output/results_individual_DEU.mat','med_ra_list','med_da_list');

% Axis Limits
min_ra_plot = -2  ;
max_ra_plot =  1  ;
min_da_plot = 0.7 ;
max_da_plot =  1  ;
min_df_plot = min_da;
max_df_plot = max_da;

% Create lists used for further analysis
md_ra_plot = med_ra_list;
md_da_plot = med_da_list;

% df by subject
md_df_plot = md_da_plot.^(1-md_ra_plot);

% Support used to compute the plots
ra_axis = linspace(min_ra,max_ra,10000);
da_axis = linspace(min_da,max_da,10000);

% CDF from pooled estimates
CDF_ra = cdf_ra(ra_axis,mu_ra_DEU,sigma_ra_DEU,min_ra,max_ra);
CDF_da = cdf_da(da_axis,alpha_da_DEU,beta_da_DEU,min_da,max_da);

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
print(fig,'Output/EstimatedDistributions_AS_cdf_ra','-depsc2');

fig = figure('Name','Estimated Distributions');
area(da_axis,CDF_da,'FaceAlpha',face_alpha_plot_1,'EdgeAlpha',edge_alpha_plot_1);
hold on;
histogram([md_da_plot],'FaceAlpha',face_alpha_plot_2,'EdgeAlpha',edge_alpha_plot_2,...
    'Normalization','cdf','NumBins',nBinsHist,'BinLimits',[min_da_plot,max_da_plot]);
hold off;
ylim([0 1]); xlim([min_da_plot,max_da_plot]);
grid('on');
print(fig,'Output/EstimatedDistributions_AS_cdf_da','-depsc2');

fig = figure('Name','Estimated Distributions');
area(df_axis_DEU,CDF_df_DEU,'FaceAlpha',face_alpha_plot_1,'EdgeAlpha',edge_alpha_plot_1);
hold on;
histogram([md_df_plot],'FaceAlpha',face_alpha_plot_2,'EdgeAlpha',edge_alpha_plot_2,...
    'Normalization','cdf','NumBins',nBinsHist,'BinLimits',[min_df_plot,max_df_plot]);
hold off;
ylim([0 1]); xlim([min_df_plot,max_df_plot]);
grid('on');
print(fig,'Output/EstimatedDistributions_AS_cdf_df','-depsc2');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Show predicted distribution of choice by risk condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% DEU %%%%
clear;
load('Output/results_baseline_DEU.mat','taskSetID','Y','menuID','subjectID','rhoX','nX','xList');

rho_Observed  = cell(6,1);
rho_Predicted = cell(6,1);
rho_All_Observed  = cell(6,1);
rho_All_Predicted = cell(6,1);

%%% Compute predicted probabilities of choosing each allocation in each lottery
for jTask = 1:6
    
    chosenData  = taskSetID == jTask ;
    Y_j         = Y(chosenData);
    menuID_j    = menuID(chosenData);
    subjectID_j = subjectID(chosenData);
    menuInTask  = unique(menuID_j)';
    
    %%% Compute predicted probabilities of choosing each allocation in each lottery
    rhoX_j = rhoX(menuInTask,:);
    
    % Predicted probability of choosing each X unconditionally
    predictedDist_X_j = sum( rhoX_j./length(menuInTask) ,1 );
    
    % Observed probability of choosing each X unconditionally
    observedDist_X_j = nan(1,nX);
    for iX = 1:nX
        observedDist_X_j(iX) = sum(Y_j==iX)./length(Y_j);
    end
    
    % Observed probability of choosing each X for each menu
    observedDist_X_Menu_j = nan(length(menuInTask),nX);
    for iMenu = menuInTask
        Y_menu = Y_j(menuID_j==iMenu);
        for iX = 1:nX
            observedDist_X_Menu_j(iMenu,iX) = sum(Y_menu==iX)./length(Y_menu);
        end
    end
    
    rho_Observed{jTask} = observedDist_X_Menu_j(menuInTask,:);
    rho_Predicted{jTask} = rhoX_j;
    rho_All_Observed{jTask} = observedDist_X_j;
    rho_All_Predicted{jTask} = predictedDist_X_j;
    
    % Print individual
    face_alpha_plot = 0.4;
    edge_alpha_plot = 0.8;
    fig = figure('Name',['PredictionByTask_DEU_',num2str(jTask)]);
    bar(xList,rho_All_Observed{jTask},'FaceAlpha',face_alpha_plot,'EdgeAlpha',edge_alpha_plot); hold on;
    bar(xList,rho_All_Predicted{jTask},'FaceAlpha',face_alpha_plot,'EdgeAlpha',edge_alpha_plot); hold off;
    if jTask == 1
        legend({'Observed', 'Predicted'},'Location','best','FontSize',12);
    end
    ylim([0,0.8]);
    grid('on');
    %     xlabel('Allocation \alpha_\tau','FontSize',12);
    ylabel('Probability','FontSize',12);
    print(fig,['Output/DEU_AS_predictedDist_',num2str(jTask)],'-depsc2');
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Scatterplot of estimates by individual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% DEU %%%%
clear;
load('Output/results_individual_DEU.mat','med_ra_list','med_da_list');

% Axis Limits
min_ra_plot = -2  ;
max_ra_plot =  1  ;
min_da_plot = 0.7 ;
max_da_plot =  1  ;

% Scatterplot ra vs da: DEU
figure1 = figure('Name','Estimates by individual: DEU');
set(figure1,'Units','Normalized');
axes1 = axes('Parent',figure1,'XGrid','on','YGrid','on'); hold(axes1,'on');
ylabel('$h$','FontSize',12,'Interpreter','latex');
xlabel('$\hat{\delta}$','FontSize',12,'Interpreter','latex');
scatter(med_da_list,med_ra_list,...
    'MarkerFaceColor',[0.301960784313725 0.749019607843137 0.929411764705882],...
    'MarkerEdgeColor',[0 0.450980392156863 0.741176470588235],...
    'LineWidth',1.5);
regLine = lsline;
set(regLine,'Parent',axes1,'LineWidth',1.2,'LineStyle','--',...
    'Color',[0.87058824300766 0.490196079015732 0]);
xlim([min_da_plot,max_da_plot]);
ylim([min_ra_plot,max_ra_plot]);
print(figure1,'Output/IndividualScatter_AS_DEU','-depsc2');
