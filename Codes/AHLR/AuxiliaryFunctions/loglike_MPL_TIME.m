function [log_like,rho_L1,rho_CHECK,pdfCHECK] = loglike_MPL_TIME(parameters,Y,menuID,algoList)
% Log-likelihood of DEU model with independent preferences and MPL data
% Delay aversion only

% Load vectors describing each menu
c1_L1 = algoList.menuList.c1_L1; c1_L2 = algoList.menuList.c1_L2;
t_L1  = algoList.menuList.t_L1 ; t_L2  = algoList.menuList.t_L2 ;

% Objects used in estimation
nL = algoList.nL;
taskSetID = algoList.taskSetID;
max_da = algoList.max_da;
min_da = algoList.min_da;

% Recover parameters
alpha_da = parameters(1); % Parameter "alpha" of the distribution of da
beta_da  = parameters(2); % Parameter "beta" of the distribution of da

% Set tremble probability
nu = eps;

% Restrictions on parameters
badParameter = alpha_da<=0 || beta_da<=0 || (alpha_da<1 && beta_da<1); % Last condition imposes uni-modality

if badParameter == 1 % If the restrictions on the parameters are not satisfied, return a nan
    
    log_like = nan;
    
else
    
    % Compute the probability of choosing option 1 for each lottery
    rho_L1 = nan(nL,1);
    for iMenu = 1:nL
        if taskSetID(iMenu) == 2 % Check time tasks only
            
            % Compute integral using a Quasi-Montecarlo method
            if t_L1(iMenu)>t_L2(iMenu)
                da_x = (c1_L2(iMenu)/c1_L1(iMenu)).^(1/(t_L1(iMenu)/30-t_L2(iMenu)/30));
                rho_L1(iMenu) = 1-cdf_da(da_x,alpha_da,beta_da,min_da,max_da);
            elseif t_L1(iMenu)<t_L2(iMenu)
                da_x = (c1_L1(iMenu)/c1_L2(iMenu)).^(1/(t_L2(iMenu)/30-t_L1(iMenu)/30));
                rho_L1(iMenu) = cdf_da(da_x,alpha_da,beta_da,min_da,max_da);
            else
                rho_L1(iMenu) = 0;
            end
            
        end
    end
    
    % Make sure each lottery is chosen with positive probability
    rho_CHECK = rho_L1;
    pdfCHECK  = cdf_da(max_da,alpha_da,beta_da,min_da,max_da);
    rho_L1    = rho_L1*(1-nu) + nu*0.5;
    
    % Compute the log-likelihood
    PC = rho_L1(menuID);
    
    % Contribution to the log-likelihood of those who choose...
    loglike1 = log(PC).*(Y==1)                          ; % Lottery 1
    loglike2 = log(1-PC).*(Y==0)                        ; % Lottery 2
    loglikeI = ( 0.5*log(PC) + 0.5*log(1-PC) ).*(Y==-1) ; % Indifferent
    
    % Log-likelihood
    log_like = mean( loglike1+loglike2+loglikeI );
    
    
    
    
end