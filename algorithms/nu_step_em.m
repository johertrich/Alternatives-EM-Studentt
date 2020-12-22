% This code belongs to the paper
%
% M. Hasannasab, J. Hertrich, F. Laus, and G. Steidl. 
% Alternatives to the EM algorithm for ML-estimation of location, scatter
% matrix and degree of freedom of the student-t distribution.
% Numerical Algorithms, 2020.
% DOI: https://doi.org/10.1007/s11075-020-00959-w
%
% If you use this code, please cite the paper.
%
% This function implements the nu-step of the EM algorithm and the
% accelerated EM algorithm
% INPUTS:
%       nu_r    - current nu
%       delta_r - 1 x n array, where delta_r(i)=(x_i-mu_r)'*sigma_r^(-1)*(x_i-mu_r)
%       d       - dimension of the data
%       n       - number of samples
% OUTPUTS:
%       nu_r_plus_one - next nu
%       evals         . number of evaluations of the subfunctions A and B
function [nu_r_plus_one,evals]=nu_step_em(nu_r,delta_r,d,n)
% define functions
A=@(x)psi(x/2)-log(x/2);
B=@(x)sum((x+d)./(x+delta_r) - log((x+d)./(x+delta_r)) -1)/n;
der_A=@(x).5*psi(1,x/2)-1./x;
b_nu_r=B(nu_r);
f=@(nu)A(nu)-A(nu_r+d)+b_nu_r;

[nu_r_plus_one,evals]=newton(nu_r,f,der_A);
end
