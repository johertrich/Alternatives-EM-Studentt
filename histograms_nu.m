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
% This script generates the histograms of the nu estimates from Figure 3 in
% the paper

clear all
addpath('algorithms')
niter=1000;
d=2;
nus_gmmf=zeros(1,niter);

n=1000;
delta=@(X,mu,sigma)sum(((sigma^(-1))*(X-repmat(mu,1,n))).*(X-repmat(mu,1,n)),1);
L=@(nu,mu,sigma,w,delta)-2*log(gamma((d+nu)/2))+2*log(gamma(nu/2))-nu.*log(nu)+(d+nu).*sum(w.*log(nu+delta))+log(det(sigma));
mu=zeros(d,1);
sigma=[1,0;0,1];

w=ones(1,n)/n;

anz_steps=300;

for nu=[1,2,5,10,100,200] 
    for iter=1:niter
        
        X=sample(n,mu,nu,sigma);
        
        % Estimating using the GMMF
        [~,nu_gmmf,~]=iterate_studentT(X,w,'GMMF',anz_steps,1);
        
        nus_gmmf(iter)=nu_gmmf;
        if mod(iter,10)==0
            disp(['step ' num2str(iter)])
        end
    end
    
    switch nu
        case 1
            histogram(nus_gmmf,'BinWidth',0.025,'BinLimits',[0.8625,1.1875])
        case 2
            histogram(nus_gmmf,'BinWidth',0.05,'BinLimits',[1.575,2.425]);
        case 5
            histogram(nus_gmmf,'BinWidth',0.2,'BinLimits',[3.7,8.9]);
        case 10
            histogram(nus_gmmf,'BinWidth',1,'BinLimits',[5.5,45.5]);
        case 100
            histogram(nus_gmmf,'BinWidth',5,'BinLimits',[12.5,217.5])
        case 200
            histogram(nus_gmmf,'BinWidth',10,'BinLimits',[15,255]);
    end
    saveas(gcf,['hist_nu' num2str(nu) '.png'])
end
