% This code belongs to the paper
%
% M. Hasannasab, J. Hertrich, F. Laus, and G. Steidl. 
% Alternatives to the EM algorithm for ML-estimation of location, scatter
% matrix and degree of freedom of the student-t distribution.
% ArXiv preprint arXiv:1910.06623, 2019.
%
% If you use this code, please cite the paper.
%
% This function implements the nu-step of the GMMF
% INPUTS:
%       nu_r    - current nu
%       delta_r - 1 x n array, where delta_r(i)=(x_i-mu_r)'*sigma_r^(-1)*(x_i-mu_r)
%       d       - dimension of the data
%       n       - number of samples
% OUTPUTS:
%       nu_r_plus_one - next nu
%       evals         . number of evaluations of the subfunctions A and B

function [nu_r_plus_one,evals]=nu_step_gmmf(nu_r,delta_r,d,n)
% define functions
A=@(x)psi(x/2)-log(x/2);
B=@(x)sum((x+d)./(x+delta_r) - log((x+d)./(x+delta_r)) -1)/n;
der_A=@(x).5*psi(1,x/2)-1./x;

% find zero
evals=0;
max_steps=100;steps=0;tol=1e-5;nu_r_plus_one=nu_r;
while steps<max_steps
    nu_k=nu_r_plus_one;
    b_nu_k=B(nu_k);
    f=@(nu)A(nu)-A(nu+d)+b_nu_k;
    if abs(f(nu_k))<tol
        break
    end
    der_f=@(nu)der_A(nu)-der_A(nu+d);
    [nu_r_plus_one,ev]=newton(nu_r_plus_one,f,der_f);
    evals=evals+ev;
    steps=steps+1;
    if nu_r_plus_one>=1e5
        nu_r_plus_one=1e5;
        break;
    end
end
if isnan(nu_r_plus_one)
    disp('...')
end
%disp(steps)
% if eps>=tol
%     disp(['Algorithm did not converge in ' num2str(max_steps) ' steps. Eps: ' num2str(eps)])
% end
end