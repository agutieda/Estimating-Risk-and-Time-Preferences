function [mean_x,median_x,mode_x,std_x,Z_x] = mom_da(alpha_da,beta_da,a,b)
% Function to compute some of the moments of dr
%
% This file: Truncated Beta Distribution
% The distribution is characterized by shape parameters "alpha_da" and "beta_da",
% lower bound "a" and upper bound "b"

% Make sure parameters are positive
if alpha_da<=0 
    alpha_da = 1e-8;
end

if beta_da<=0
    beta_da = 1e-8;
end

Z_x = betacdf(b,alpha_da,beta_da)-betacdf(a,alpha_da,beta_da);

% Mean and Std. Deviation
[mean_x,var_x] = betastat(alpha_da,beta_da);
std_x = sqrt(var_x);

% Median
if alpha_da>1 && beta_da>1
    median_x = ( alpha_da-(1/3) )/( alpha_da+beta_da-(2/3) ) ;
elseif alpha_da==beta_da
    median_x = 0.5;
elseif alpha_da<1 && beta_da<1
    median_x = nan;
elseif alpha_da<0.01 || beta_da > 1e4
    median_x = 0;
elseif  beta_da<0.01 || alpha_da > 1e4
    median_x = 1;
else    
    % Brute force
    auxFun = @(x) betainc(x,alpha_da,beta_da)-0.5;
    options = optimset('Display','off'); 
    median_x = fzero(auxFun,[0,1],options);   
end

% Mode
if alpha_da>1 && beta_da>1
    mode_x = (alpha_da-1) / (alpha_da+beta_da-2);
elseif alpha_da==1 && beta_da==1
    mode_x = 0.5; % Any value in (0,1)
elseif alpha_da<=1 && beta_da>1
    mode_x = 0;
elseif alpha_da>1 && beta_da<=1
    mode_x = 1;
else
    mode_x = nan;
end

end
