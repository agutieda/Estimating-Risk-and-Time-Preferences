function [log_like,rho_L1,rho_CHECK,pdfCHECK] = loglike_MPL_RISK(parameters,Y,menuID,algoList)
% Log-likelihood of DEU model with independent preferences and MPL data
% Risk aversion only

% Objects used in estimation
nL = algoList.nL;
chosenX = algoList.chosenX;
taskSetID = algoList.taskSetID;

% Recover parameters
mu_ra    = parameters(1); % Parameter "mu" of the distribution of ra
sigma_ra = parameters(2); % Parameter "sigma" of distribution of ra

% Set tremble probability
nu = eps;

% Restrictions on parameters
badParameter = sigma_ra<=0;

if badParameter == 1 % If the restrictions on the parameters are not satisfied, return a nan
    
    log_like = nan;
    
else
    
    % Variables used for numerical integration
    raPoints = algoList.raPoints ; % Evaluation points: rm
    intWeights_ra = algoList.intWeights_ra ; % Integration weights
    max_ra = algoList.max_ra;
    min_ra = algoList.min_ra;
    
    % Compute value of pdf of joint distribution at given points
    PDF_ra = pdf_ra( raPoints , mu_ra , sigma_ra , min_ra , max_ra) ;
    
    % Make sure integration is working
    pdfCHECK = intWeights_ra'*PDF_ra;
    
    if abs(pdfCHECK-1)>0.1
        
        log_like = nan;
        warning('Integration is not working! Increase number of points!');
        
    else
        
        % Compute the probability of choosing option 1 for each lottery
        rho_L1 = nan(nL,1);
        for iMenu = 1:nL
            if taskSetID(iMenu) == 1 % Check risk tasks only
                integrand_value = chosenX(:,iMenu).*PDF_ra;
                rho_L1(iMenu) = (intWeights_ra')*integrand_value;
            end
        end
        
        % Add tremble
        rho_CHECK = rho_L1;        
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
    
    
end