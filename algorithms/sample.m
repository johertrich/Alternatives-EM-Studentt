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
% This function generates samples of a multivariate Student-t distribution.
function X=sample(n,mu,nu,sigma)
[a,b]=eig(sigma);
sqrt_sigma=a'*(b.^.5)*a;
d=size(mu,1);
X=repmat(mu,1,n)+(sqrt_sigma*normrnd(0,1,[d,n]))./repmat((gamrnd(nu/2,2/nu,[1,n])).^.5,d,1);
end
