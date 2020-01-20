function pdf_x = pdf_da(x,alpha,beta,a,b)
% Function to evaluate the pdf of dr
%
% This file: Truncated Beta Distribution
% The distribution is characterized by shape parameters "alpha" and "beta",
% lower bound "a" and upper bound "b"

Z = betacdf(b,alpha,beta)-betacdf(a,alpha,beta);
pdf_x = betapdf(x,alpha,beta) ./ Z;
pdf_x(x<a | x>b)=0;


end