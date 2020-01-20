function Y = par2mom(X,algoList)
% Function that computes the mean, median, mode and standard deviation of
% the estimated distirbutions, given values of scale and location parameters

max_ra = algoList.max_ra;
min_ra = algoList.min_ra;
max_da = algoList.max_da;
min_da = algoList.min_da;

mu_ra     = X(1);
sigma_ra  = X(2);
alpha_da  = X(3);
beta_da   = X(4);

% Moments of distribution of risk aversion
[mean_ra,median_ra,mode_ra,sd_ra,Z_ra] = mom_ra(mu_ra,sigma_ra,min_ra,max_ra);

% Moments of distribution of discount rate
[mean_da,median_da,mode_da,sd_da,Z_da] = mom_da(alpha_da,beta_da,min_da,max_da);

% Estimated discount factor and it's implied return rate
median_df = median_da^(1-median_ra);
median_df_rr = 100*( (1/median_df)-1 );

% Output
Y = [
    mean_ra;median_ra;mode_ra;sd_ra;Z_ra;
    mean_da;median_da;mode_da;sd_da;Z_da;
    median_df;median_df_rr;
    ];

end
