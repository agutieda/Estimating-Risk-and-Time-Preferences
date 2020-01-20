function [log_like,rho_L1,rho_CHECK,pdfCHECK] = loglike_MPL_DEU_CORR(parameters,Y,menuID,algoList)
% Log-likelihood of DEU model with correlated preferences and MPL data

% Objects used in estimation
nL = algoList.nL;
t_L1 = algoList.t_L1 ;
t_L2 = algoList.t_L2 ;
chosenX = algoList.chosenX;

% Recover parameters
mu_ra    = parameters(1); % Parameter "mu" of the distribution of ra
sigma_ra = parameters(2); % Parameter "sigma" of distribution of ra
alpha_da = parameters(3); % Parameter "alpha" of the distribution of da
beta_da  = parameters(4); % Parameter "beta" of the distribution of da
rho_rada = parameters(5); % Parameter "rho" characterizing correlation between ra and da

% Set tremble probability
nu = eps;

% Restrictions on parameters
badParameter = sigma_ra<=0 || alpha_da<=0 || beta_da<=0 || abs(rho_rada)>=1 || (alpha_da<1 && beta_da<1) ; % Last condition imposes uni-modality;

if badParameter == 1 % If the restrictions on the parameters are not satisfied, return a nan
    
    log_like = nan;
    rho_L1 = nan;
    rho_CHECK = nan;
    pdfCHECK = nan;
    
    
else
    
    % Variables used for numerical integration
    max_ra = algoList.max_ra;
    min_ra = algoList.min_ra;
    max_da = algoList.max_da;
    min_da = algoList.min_da;
    raPoints = algoList.raPoints;
    daPoints = algoList.daPoints;
    intWeights = algoList.intWeights;
    intWeights_ra = algoList.intWeights_ra;
    
    % Compute value of pdf of joint distribution at given points
    PDF_ra = pdf_ra( raPoints , mu_ra    , sigma_ra , min_ra , max_ra) ;
    PDF_da = pdf_da( daPoints , alpha_da , beta_da  , min_da , max_da) ;
    CDF_ra = cdf_ra( raPoints , mu_ra    , sigma_ra , min_ra , max_ra) ;
    CDF_da = cdf_da( daPoints , alpha_da , beta_da  , min_da , max_da) ;
    copulaC = copulapdf('Gaussian',[CDF_ra,CDF_da],rho_rada);
    PDF = PDF_ra.*PDF_da.*copulaC;
    
    % Make sure integration is working
    pdfCHECK = intWeights'*PDF;
    
    if abs(pdfCHECK-1)>0.1
        
        log_like = nan;
        rho_L1 = nan;
        rho_CHECK = nan;
        warning('Integration is not working! Increase number of points!');
        
    else
        
        % Compute the probability of choosing option 1 for each lottery
        rho_L1 = nan(nL,1);
        for iMenu = 1:nL
            
            if t_L1(iMenu) == t_L2(iMenu)
                integrand_value = chosenX(:,iMenu).*PDF_ra;
                rho_L1(iMenu) = (intWeights_ra')*integrand_value;
            else
                integrand_value = chosenX(:,iMenu).*PDF;
                rho_L1(iMenu) = (intWeights')*integrand_value;
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