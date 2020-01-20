function [log_like,rhoX,rhoCHECK,pdfCHECK] = loglike_CTB_DEU(parameters,Y,menuID,algoList)
% Log-likelihood of DEU model with independent preferences and CTB data

% Load
nX = algoList.nX;
nL = algoList.nL;
chosenX = algoList.chosenX;

% Recover parameters
mu_ra    = parameters(1); % Parameter "mu" of the distribution of ra
sigma_ra = parameters(2); % Parameter "sigma" of distribution of ra
alpha_da = parameters(3); % Parameter "alpha" of the distribution of da
beta_da  = parameters(4); % Parameter "beta" of the distribution of da

% Set tremble probability
nu = eps;

% Restrictions on parameters
badParameter = sigma_ra<=0 || alpha_da<=0 || beta_da<=0 || (alpha_da<1 && beta_da<1); % Last condition imposes uni-modality;

if badParameter == 1 % If the restrictions on the parameters are not satisfied, return a nan
    
    log_like = nan;
    rhoX     = nan;
    rhoCHECK = nan;
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
    
    % Compute value of pdf of joint distribution at given points
    PDF_ra = pdf_ra( raPoints , mu_ra    , sigma_ra , min_ra , max_ra) ;
    PDF_da = pdf_da( daPoints , alpha_da , beta_da  , min_da , max_da) ;
    PDF = PDF_ra.*PDF_da;
    
    % Make sure integration is working
    pdfCHECK = intWeights'*PDF;
    
    if abs(pdfCHECK-1)>0.1
        
        log_like = nan;
        rhoX     = nan;
        rhoCHECK = nan;
        pdfCHECK = nan;
        warning('Integration is not working! Increase number of points!');
        
    else
        
        % Create appropiate matrices for integration
        PDF_mat = PDF*ones(1,nX);
        
        % Compute the probability of allocations 1 to nX for each lottery
        rhoX = nan(nL,nX);
        for iMenu = 1:nL
            integrand_value = chosenX{iMenu}.*PDF_mat;
            rhoX(iMenu,:) =  (intWeights')*integrand_value;
        end
        
        % Add tremble
        rhoCHECK = rhoX;
        rhoX     = rhoX*(1-nu) + nu*(1/nX);
        
        % Compute the log-likelihood
        idxLogLike  = sub2ind([nL,nX],menuID,Y);
        rhoX_chosen = rhoX(idxLogLike);
        log_like    = mean( log(rhoX_chosen) );
        
    end
    
    
end
