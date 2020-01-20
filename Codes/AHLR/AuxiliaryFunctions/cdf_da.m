function cdf_x = cdf_da(x,alpha,beta,a,b)
% Function to evaluate the cdf of dr
%
% This file: Truncated Beta Distribution
% The distribution is characterized by shape parameters "alpha" and "beta",
% lower bound "a" and upper bound "b"

Z = betacdf(b,alpha,beta)-betacdf(a,alpha,beta);
cdf_x = (betacdf(x,alpha,beta)-betacdf(a,alpha,beta))./Z;
cdf_x(x<a)=0;
cdf_x(x>b)=1;

end