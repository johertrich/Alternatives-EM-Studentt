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
% This function implements newtons method for determining the zeros in the nu-steps.
%
% INPUTS:
%   start       - initial value
%   f           - function handle of the function
%   der_f       - function handle of the derivative
%   max_steps   - maximum number of steps
%   tol         - stopping criteria. Stop if abs(f(val))<tol.
%
% OUTPUTS:
%   zero        - root of f
%   evals       - number of evaluations of f and its derivative.
%
function [zero,evals]=newton(start,f,der_f,max_steps,tol)
if nargin < 5
    tol=1e-5;
end
if nargin < 4
    max_steps=1000;
end
rep=0;
zero=start;steps=0; 
f_val=f(zero);
eps=abs(f_val);
while steps<max_steps && eps>=tol
    newzero=zero-f_val/der_f(zero);
    zero=newzero;
    if zero<0
        zero=10^(-2);
        rep=rep+1;
        if rep>10
            disp('One nu was reprojected on the interval [1e-2,inf[');
            break;
        end
    end
    f_val=f(zero);
    eps=abs(f_val);
    steps=steps+1;
    if mod(steps,100)==0
        disp(['Reached step ' num2str(steps) ' in Newton. Eps: ' num2str(eps) ' with argument ' num2str(zero) ' and derivative ' num2str(der_f(zero))])
    end
end
if eps>=tol
    disp(['Newton did not converge in ' num2str(max_steps) ' steps. Eps: ' num2str(eps) ' with argument ' num2str(zero) ' and derivative ' num2str(der_f(zero))])
end
evals=2*steps+1;
end
