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
% This script generates the plots of the negative log-likelihood functions
% as in Figure 2 in the paper

close all;
addpath('algorithms')
n=1000;
d=2;
delta=@(X,mu,sigma)sum(((sigma^(-1))*(X-repmat(mu,1,n))).*(X-repmat(mu,1,n)),1);
L=@(nu,mu,sigma,w,delta)-2*log(gamma((d+nu)/2))+2*log(gamma(nu/2))-nu.*log(nu)+(d+nu).*sum(w.*log(nu+delta))+log(det(sigma));
mu=zeros(d,1);
sigma=[1,0;0,1];
w=ones(1,n)/n;

for nu=[1,2,5,10,100,200]
    

    X=sample(n,mu,nu,sigma);
    
    anz_steps=100;


    % Estimating using the GMMF
    [~,~,~,~,~,obj_gmmf]=iterate_studentT(X,w,'GMMF',anz_steps,0,0,0,1);

    % Estimating using the MMF
    [~,~,~,~,~,obj_mmf]=iterate_studentT(X,w,'MMF',anz_steps,0,0,0,1);

    % Estimating using the EM-Algorithm
    [~,~,~,~,~,obj_em]=iterate_studentT(X,w,'EM',anz_steps,0,0,0,1);

    % Estimating using the accelerated EM-Algorithm
    [~,~,~,~,~,obj_aem]=iterate_studentT(X,w,'aEM',anz_steps,0,0,0,1);


    semilogx(obj_em,'b--','LineWidth',1);
    hold on;
    semilogx(obj_aem,'k:','LineWidth',1);
    semilogx(obj_mmf,'-.','LineWidth',1,'Color','#77AC30');
    semilogx(obj_gmmf,'r-','LineWidth',1);
    hold off;
    
    
    legend('EM','aEM','MMF','GMMF')
    
    saveas(gcf,['Like_nu' num2str(nu) '.png'])
end
