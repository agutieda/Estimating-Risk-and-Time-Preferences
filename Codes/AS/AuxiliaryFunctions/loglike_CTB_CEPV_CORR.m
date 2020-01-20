function [log_like,rhoX,rhoCHECK,pdfCHECK] = loglike_CTB_CEPV_CORR(parameters,Y,menuID,algoList)
% Log-likelihood of CEPV model with correlated preferences and CTB data

% Load
nX = algoList.nX;
nL = algoList.nL;
chosenX = algoList.chosenX;

% Recover parameters
mu_ra    = parameters(1); % Parameter "mu" of the log-normal distribution of ra
sigma_ra = parameters(2); % Parameter "sigma" of the log-normal distribution of ra
alpha_da = parameters(3); % Parameter "alpha" of the distribution of da
beta_da  = parameters(4); % Parameter "beta" of the distribution of da
mu_is    = parameters(5); % Parameter "mu" of the log-normal distribution of is
sigma_is = parameters(6); % Parameter "sigma" of the log-normal distribution of is
rho_rada = parameters(7); % Correlation coefficient between ra and da
rho_rais = parameters(8); % Correlation coefficient between ra and is
rho_dais = parameters(9); % Correlation coefficient between da and is

% Set tremble probability
nu = eps;

% Correlation matrix
P = [
    1      ,  rho_rada ,  rho_rais ;
    rho_rada  ,      1    ,  rho_dais ;
    rho_rais  ,  rho_dais ,     1
    ];

% Check that P is positive semidefinite
eigenvalP = eig(P);
PSDcheck = all(eigenvalP>=0);

% Restrictions on parameters
badParameter = sigma_ra<=0 || sigma_is<=0 || PSDcheck==0 || ...
    abs(rho_rada)>=1 || abs(rho_rais)>=1 || abs(rho_dais)>=1 || ...
    alpha_da<=0 || beta_da<=0 || (alpha_da<1 && beta_da<1); % Last condition imposes uni-modality;

if badParameter == 1 % If the restrictions on the parameters are not satisfied, return a nan
    
    log_like = nan;
    rhoX = nan;
    rhoCHECK = nan;
    pdfCHECK = nan;
    
    
else
    
    % Variables used for numerical integration
    raPoints   = algoList.raPoints   ; % Evaluation points: ra
    daPoints   = algoList.daPoints   ; % Evaluation points: da
    isPoints   = algoList.isPoints   ; % Evaluation points: is
    intWeights = algoList.intWeights ; % Integration weights
    
    % Integration bounds
    max_ra = algoList.max_ra;
    min_ra = algoList.min_ra;
    max_da = algoList.max_da;
    min_da = algoList.min_da;
    max_is = algoList.max_is;
    min_is = algoList.min_is;
    
    % Compute value of pdf of joint distribution at given points
    PDF_ra = pdf_ra( raPoints , mu_ra    , sigma_ra , min_ra , max_ra) ;
    PDF_da = pdf_da( daPoints , alpha_da , beta_da  , min_da , max_da) ;
    PDF_is = pdf_is( isPoints , mu_is    , sigma_is , min_is , max_is) ;
    CDF_ra = cdf_ra( raPoints , mu_ra    , sigma_ra , min_ra , max_ra) ;
    CDF_da = cdf_da( daPoints , alpha_da , beta_da  , min_da , max_da) ;
    CDF_is = cdf_is( isPoints , mu_is    , sigma_is , min_is , max_is) ;
    copulaC = copulapdf('Gaussian',[CDF_ra,CDF_da,CDF_is],P);
    PDF = PDF_ra.*PDF_da.*PDF_is.*copulaC;
    
    % Make sure integration is working
    pdfCHECK = intWeights'*PDF;
    
    if abs(pdfCHECK-1)> 0.1
        
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
