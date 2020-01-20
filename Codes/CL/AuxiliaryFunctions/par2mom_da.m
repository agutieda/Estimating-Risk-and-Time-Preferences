function Y = par2mom_da(X,algoList)
% Function that computes the mean, median, mode and standard deviation of a
% truncated normal distribution, given values of scale and location parameters

max_da = algoList.max_da;
min_da = algoList.min_da;

alpha_da     = X(1);
beta_da  = X(2);

% Moments of distribution of discount rate
[mean_da,median_da,mode_da,sd_da,Z_da] = mom_da(alpha_da,beta_da,min_da,max_da);


% Output
Y = [mean_da;median_da;mode_da;sd_da;Z_da];

end
