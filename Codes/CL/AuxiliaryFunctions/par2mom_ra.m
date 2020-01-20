function Y = par2mom_ra(X,algoList)
% Function that computes the mean, median, mode and standard deviation of
% the estimated distirbutions, given values of scale and location parameters

max_ra = algoList.max_ra;
min_ra = algoList.min_ra;

mu_ra     = X(1);
sigma_ra  = X(2);

% Moments of distribution of risk aversion
[mean_ra,median_ra,mode_ra,sd_ra,Z_ra] = mom_ra(mu_ra,sigma_ra,min_ra,max_ra);

% Output
Y = [mean_ra;median_ra;mode_ra;sd_ra;Z_ra];

end
