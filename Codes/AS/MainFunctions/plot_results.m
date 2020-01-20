% RANDOM MODELS FOR THE JOINT TREATMENT OF RISK AND TIME PREFERENCES
% by Jose Apesteguia, Miguel A. Ballester and Angelo Gutierrez
%
% Script to estimate DEU-H, PVCE-H and CEPV-H models using data from:
% "Estimating Time Preferences from Convex Budgets"
%  by James Andreoni and Charles Sprenger.
%  American Economic Review, (2012)
%
% This file:
% - Plot PDF of estimated preferences under DEU-H, PVCE-H and CEPV-H by pooling all data
% - Plot observed and predicted distributions of choices under each representation
% - Plot scatterplot of individual estimates under each representation
%
% November 2019
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

% PVCE
load('Output/results_baseline_PVCE.mat','mu_ra','sigma_ra','alpha_da','beta_da','mu_is','sigma_is');
mu_ra_PVCE    = mu_ra    ;
sigma_ra_PVCE = sigma_ra ;
alpha_da_PVCE = alpha_da ;
beta_da_PVCE  = beta_da  ;
mu_is_PVCE    = mu_is    ;
sigma_is_PVCE = sigma_is ;

% CEPV
load('Output/results_baseline_CEPV.mat','mu_ra','sigma_ra','alpha_da','beta_da','mu_is','sigma_is');
mu_ra_CEPV    = mu_ra    ;
sigma_ra_CEPV = sigma_ra ;
alpha_da_CEPV = alpha_da ;
beta_da_CEPV  = beta_da  ;
mu_is_CEPV    = mu_is    ;
sigma_is_CEPV = sigma_is ;

% Specify some bounds for the plots
min_ra = -10      ; % Lower limit of the integral on ra
max_ra =  1       ; % Upper limit of the integral on ra
min_da =  0.5     ; % Lower limit of the integral on da
max_da =  1       ; % Upper limit of the integral on da
min_is = -10      ; % Lower limit of the integral on is
max_is =  1       ; % Upper limit of the integral on is


%%%%% ESTIMATE DISTRIBUTION OF DISCOUNT FACTOR USING ESTIMATED DISTRIBUTION OF DELAY AVERSION

% Number of draws used to simulate this distribution
nDraws = 1000000;

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

%%% PVCE

% Generate draws from distribution of discount related to delay aversion
daDraws = betarnd(alpha_da_PVCE,beta_da_PVCE,nDraws,1);

% Generate draws from the distribution of intertemporal substitution
p  = rand(nDraws,1);
Ga = normcdf(min_is,mu_is_PVCE,sigma_is_PVCE);
Gb = normcdf(max_is,mu_is_PVCE,sigma_is_PVCE);
u  = Ga + p.*(Gb-Ga);
isDraws = mu_is_PVCE-sigma_is_PVCE.*sqrt(2)*erfcinv(2*u);

% Generate draws of the discount factor
dfDraws = daDraws.^(1-isDraws);

% Trim extreme values to get a better estimate of the PDF and CDF around dr=1
% dfDraws_uptail = prctile(dfDraws,99);
dfDraws_uptail = 3;
dfDraws = dfDraws(dfDraws<dfDraws_uptail);

% Estimate pdf and cdf
[CDF_df_PVCE,df_axis_PVCE] = ksdensity(dfDraws,'Support','positive','Function','cdf');
PDF_df_PVCE = ksdensity(dfDraws,df_axis_PVCE,'Support','positive','Function','pdf');

clearvars daDraws p Ga Gb u isDraws dfDraws;

%%% CEPV

% Generate draws from distribution of discount related to delay aversion
daDraws = betarnd(alpha_da_CEPV,beta_da_CEPV,nDraws,1);

% Generate draws from the distribution of intertemporal substitution
p  = rand(nDraws,1);
Ga = normcdf(min_is,mu_is_CEPV,sigma_is_CEPV);
Gb = normcdf(max_is,mu_is_CEPV,sigma_is_CEPV);
u  = Ga + p.*(Gb-Ga);
isDraws = mu_is_CEPV-sigma_is_CEPV.*sqrt(2)*erfcinv(2*u);

% Generate draws of the discount factor
dfDraws = daDraws.^(1-isDraws);

% Trim extreme values to get a better estimate of the PDF and CDF around dr=1
% dfDraws_uptail = prctile(dfDraws,99);
dfDraws_uptail = 3;
dfDraws = dfDraws(dfDraws<dfDraws_uptail);

% Estimate pdf and cdf
[CDF_df_CEPV,df_axis_CEPV] = ksdensity(dfDraws,'Support','positive','Function','cdf');
PDF_df_CEPV = ksdensity(dfDraws,df_axis_CEPV,'Support','positive','Function','pdf');

clearvars daDraws p Ga Gb u isDraws dfDraws;

%%%%% PLOT ESTIMATED DISTRIBUTIONS

% Support used to compute the plots
ra_axis = linspace(min_ra,max_ra,10000);
da_axis = linspace(min_da,max_da,10000);
is_axis = linspace(min_is,max_is,10000);

% PDF: DEU
PDF_ra_DEU = pdf_ra( ra_axis , mu_ra_DEU    , sigma_ra_DEU , min_ra , max_ra ) ;
PDF_da_DEU = pdf_da( da_axis , alpha_da_DEU , beta_da_DEU  , min_da , max_da ) ;
PDF_is_DEU = PDF_ra_DEU;

% PDF: PVCE
PDF_ra_PVCE = pdf_ra( ra_axis , mu_ra_PVCE    , sigma_ra_PVCE , min_ra , max_ra ) ;
PDF_da_PVCE = pdf_da( da_axis , alpha_da_PVCE , beta_da_PVCE  , min_da , max_da ) ;
PDF_is_PVCE = pdf_is( is_axis , mu_is_PVCE    , sigma_is_PVCE , min_is , max_is ) ;

% PDF: CEPV
PDF_ra_CEPV = pdf_ra( ra_axis , mu_ra_CEPV    , sigma_ra_CEPV , min_ra , max_ra ) ;
PDF_da_CEPV = pdf_da( da_axis , alpha_da_CEPV , beta_da_CEPV  , min_da , max_da ) ;
PDF_is_CEPV = pdf_is( is_axis , mu_is_CEPV    , sigma_is_CEPV , min_is , max_is ) ;

% Axis Limits
min_ra_plot = -3    ;
max_ra_plot =  1    ;
min_da_plot =  0.70 ;
max_da_plot =  1    ;
min_is_plot = -3    ;
max_is_plot =  1    ;
min_df_plot =  0.70 ;
max_df_plot =  1    ;

% Plot options
face_alpha_plot = 0.3;
edge_alpha_plot = 0.8;

% PDF

fig = figure('Name','Estimated Distributions');
area(ra_axis,PDF_ra_DEU,'FaceAlpha',face_alpha_plot,'EdgeAlpha',edge_alpha_plot);
hold on;
area(ra_axis,PDF_ra_PVCE,'FaceAlpha',face_alpha_plot,'EdgeAlpha',edge_alpha_plot);
hold on;
area(ra_axis,PDF_ra_CEPV,'FaceAlpha',face_alpha_plot,'EdgeAlpha',edge_alpha_plot);
hold off;
xlim([min_ra_plot,max_ra_plot]);
grid('on');
print(fig,'Output/EstimatedDistributions_AS_PDF_ra_BASELINE','-depsc2');

fig = figure('Name','Estimated Distributions');
area(da_axis,PDF_da_DEU,'FaceAlpha',face_alpha_plot,'EdgeAlpha',edge_alpha_plot);
hold on;
area(da_axis,PDF_da_PVCE,'FaceAlpha',face_alpha_plot,'EdgeAlpha',edge_alpha_plot);
hold on;
area(da_axis,PDF_da_CEPV,'FaceAlpha',face_alpha_plot,'EdgeAlpha',edge_alpha_plot);
hold off;
xlim([min_da_plot,max_da_plot]);
grid('on');
legend({'DEU','PVCE','CEPV'},'Location','Best');
print(fig,'Output/EstimatedDistributions_AS_PDF_da_BASELINE','-depsc2');

fig = figure('Name','Estimated Distributions');
area(is_axis,PDF_is_DEU,'FaceAlpha',face_alpha_plot,'EdgeAlpha',edge_alpha_plot);
hold on;
area(is_axis,PDF_is_PVCE,'FaceAlpha',face_alpha_plot,'EdgeAlpha',edge_alpha_plot);
hold on;
area(is_axis,PDF_is_CEPV,'FaceAlpha',face_alpha_plot,'EdgeAlpha',edge_alpha_plot);
hold off;
xlim([min_is_plot,max_is_plot]);
grid('on');
print(fig,'Output/EstimatedDistributions_AS_PDF_is_BASELINE','-depsc2');

fig = figure('Name','Estimated Distributions');
area(df_axis_DEU,PDF_df_DEU,'FaceAlpha',face_alpha_plot,'EdgeAlpha',edge_alpha_plot);
hold on;
area(df_axis_PVCE,PDF_df_PVCE,'FaceAlpha',face_alpha_plot,'EdgeAlpha',edge_alpha_plot);
hold on;
area(df_axis_CEPV,PDF_df_CEPV,'FaceAlpha',face_alpha_plot,'EdgeAlpha',edge_alpha_plot);
hold off;
xlim([min_df_plot,max_df_plot]);
grid('on');
print(fig,'Output/EstimatedDistributions_AS_PDF_df_BASELINE','-depsc2');

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


%%%% PVCE %%%%
clear;
load('Output/results_baseline_PVCE.mat','taskSetID','Y','menuID','subjectID','rhoX','nX','xList');

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
    fig = figure('Name','PredictionByTask');
    bar(xList,rho_All_Observed{jTask},'FaceAlpha',face_alpha_plot,'EdgeAlpha',edge_alpha_plot); hold on;
    bar(xList,rho_All_Predicted{jTask},'FaceAlpha',face_alpha_plot,'EdgeAlpha',edge_alpha_plot); hold off;
    if jTask == 1
        legend({'Observed','Predicted'},'Location','best','FontSize',12);
    end
    ylim([0,0.8]);
    grid('on');
    %     xlabel('Allocation \alpha_\tau','FontSize',12);
    ylabel('Probability','FontSize',12);
    print(fig,['Output/PVCE_AS_predictedDist_',num2str(jTask)],'-depsc2');
    
end



%%%% CEPV %%%%
clear;
load('Output/results_baseline_CEPV.mat','taskSetID','Y','menuID','subjectID','rhoX','nX','xList');

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
    fig = figure('Name','PredictionByTask');
    bar(xList,rho_All_Observed{jTask},'FaceAlpha',face_alpha_plot,'EdgeAlpha',edge_alpha_plot); hold on;
    bar(xList,rho_All_Predicted{jTask},'FaceAlpha',face_alpha_plot,'EdgeAlpha',edge_alpha_plot); hold off;
    if jTask == 1
        legend({'Observed','Predicted'},'Location','best','FontSize',12);
    end
    ylim([0,0.8]);
    grid('on');
    %     xlabel('Allocation \alpha_\tau','FontSize',12);
    ylabel('Probability','FontSize',12);
    print(fig,['Output/CEPV_AS_predictedDist_',num2str(jTask)],'-depsc2');
    
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
set(figure1,'Units','Normalized','Position', [0.3 0.3 0.25 0.5]);
axes1 = axes('Parent',figure1,'XGrid','on','YGrid','on'); hold(axes1,'on');
ylabel('$h$','FontSize',16,'Interpreter','latex');
xlabel('$\hat{\delta}$','FontSize',16,'Interpreter','latex');
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


%%%% PVCE %%%%%
clear;
load('Output/results_individual_PVCE.mat','med_ra_list','med_da_list','med_is_list');

% Axis Limits
min_ra_plot = -3  ;
max_ra_plot =  1  ;
min_da_plot = 0.5 ;
max_da_plot =  1  ;
min_is_plot = -1  ;
max_is_plot =  1  ;

% Plot a scatterplot of the median estimates
figure1 = figure('Name','Estimates by individual: PVCE');
set(figure1,'Units','Normalized','Position', [0.1 0.35 0.8 0.45]);

axes1 = subplot(1,3,1);
scatter(axes1,med_da_list,med_ra_list,...
    'MarkerFaceColor',[0.301960784313725 0.749019607843137 0.929411764705882],...
    'MarkerEdgeColor',[0 0.450980392156863 0.741176470588235],...
    'LineWidth',1.5);
regLine = lsline;
set(regLine,'Parent',axes1,'LineWidth',1.2,'LineStyle','--',...
    'Color',[0.87058824300766 0.490196079015732 0]);
ylabel(axes1,'$r$','FontSize',16,'Interpreter','latex');
ylim([min_ra_plot,max_ra_plot]);
xlabel(axes1,'$\hat{\delta}$','FontSize',16,'Interpreter','latex');
xlim([min_da_plot,max_da_plot]);
grid(axes1,'on')

axes2 = subplot(1,3,2);
scatter(axes2,med_is_list,med_ra_list,...
    'MarkerFaceColor',[0.301960784313725 0.749019607843137 0.929411764705882],...
    'MarkerEdgeColor',[0 0.450980392156863 0.741176470588235],...
    'LineWidth',1.5);
regLine = lsline;
set(regLine,'Parent',axes2,'LineWidth',1.2,'LineStyle','--',...
    'Color',[0.87058824300766 0.490196079015732 0]);
grid(axes2,'on')
ylabel(axes2,'$r$','FontSize',16,'Interpreter','latex');
ylim([min_ra_plot,max_ra_plot]);
xlabel(axes2,'$\eta$','FontSize',16,'Interpreter','latex');
xlim([min_is_plot,max_is_plot]);

axes3 = subplot(1,3,3);
scatter(axes3,med_da_list,med_is_list,...
    'MarkerFaceColor',[0.301960784313725 0.749019607843137 0.929411764705882],...
    'MarkerEdgeColor',[0 0.450980392156863 0.741176470588235],...
    'LineWidth',1.5);
regLine = lsline;
set(regLine,'Parent',axes3,'LineWidth',1.2,'LineStyle','--',...
    'Color',[0.87058824300766 0.490196079015732 0]);
grid(axes3,'on')
ylabel(axes3,'$\eta$','FontSize',16,'Interpreter','latex');
ylim([min_is_plot,max_is_plot]);
xlabel(axes3,'$\hat{\delta}$','FontSize',16,'Interpreter','latex');
xlim([min_da_plot,max_da_plot]);

print(figure1,'Output/IndividualScatter_AS_PVCE','-depsc2');


%%%% CEPV %%%%%
clear;
load('Output/results_individual_CEPV.mat','med_ra_list','med_da_list','med_is_list');

% Axis Limits
min_ra_plot = -3  ;
max_ra_plot =  1  ;
min_da_plot = 0.5 ;
max_da_plot =  1  ;
min_is_plot = -1  ;
max_is_plot =  1  ;

% Plot a scatterplot of the median estimates
figure1 = figure('Name','Estimates by individual: CEPV');
set(figure1,'Units','Normalized','Position', [0.1 0.35 0.8 0.45]);

axes1 = subplot(1,3,1);
scatter(axes1,med_da_list,med_ra_list,...
    'MarkerFaceColor',[0.301960784313725 0.749019607843137 0.929411764705882],...
    'MarkerEdgeColor',[0 0.450980392156863 0.741176470588235],...
    'LineWidth',1.5);
regLine = lsline;
set(regLine,'Parent',axes1,'LineWidth',1.2,'LineStyle','--',...
    'Color',[0.87058824300766 0.490196079015732 0]);
ylabel(axes1,'$r$','FontSize',16,'Interpreter','latex');
ylim([min_ra_plot,max_ra_plot]);
xlabel(axes1,'$\hat{\delta}$','FontSize',16,'Interpreter','latex');
xlim([min_da_plot,max_da_plot]);
grid(axes1,'on')

axes2 = subplot(1,3,2);
scatter(axes2,med_is_list,med_ra_list,...
    'MarkerFaceColor',[0.301960784313725 0.749019607843137 0.929411764705882],...
    'MarkerEdgeColor',[0 0.450980392156863 0.741176470588235],...
    'LineWidth',1.5);
regLine = lsline;
set(regLine,'Parent',axes2,'LineWidth',1.2,'LineStyle','--',...
    'Color',[0.87058824300766 0.490196079015732 0]);
grid(axes2,'on')
ylabel(axes2,'$r$','FontSize',16,'Interpreter','latex');
ylim([min_ra_plot,max_ra_plot]);
xlabel(axes2,'$\eta$','FontSize',16,'Interpreter','latex');
xlim([min_is_plot,max_is_plot]);

axes3 = subplot(1,3,3);
scatter(axes3,med_da_list,med_is_list,...
    'MarkerFaceColor',[0.301960784313725 0.749019607843137 0.929411764705882],...
    'MarkerEdgeColor',[0 0.450980392156863 0.741176470588235],...
    'LineWidth',1.5);
regLine = lsline;
set(regLine,'Parent',axes3,'LineWidth',1.2,'LineStyle','--',...
    'Color',[0.87058824300766 0.490196079015732 0]);
grid(axes3,'on')
ylabel(axes3,'$\eta$','FontSize',16,'Interpreter','latex');
ylim([min_is_plot,max_is_plot]);
xlabel(axes3,'$\hat{\delta}$','FontSize',16,'Interpreter','latex');
xlim([min_da_plot,max_da_plot]);

print(figure1,'Output/IndividualScatter_AS_CEPV','-depsc2');
