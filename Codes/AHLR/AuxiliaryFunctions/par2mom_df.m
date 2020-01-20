function Y = par2mom_df(X)
% Function that computes a point estimate of the discount factor and it's
% implied rate of return in annual terms

median_ra = X(1);
median_da = X(2);

% Estimated discount factor and it's implied return rate
median_df = median_da^(1-median_ra);
median_df_rr = 100*( (1/median_df)-1 );

% Output
Y = [median_df;median_df_rr];


end
