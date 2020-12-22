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
% This script performs a simulation study as in Table 1 and 2 in the paper and
% saves the result in a text document called 'results'

clear all
addpath('algorithms')
delete results;
d=2;
n=1000;
delta=@(X,mu,sigma)sum(((sigma^(-1))*(X-repmat(mu,1,n))).*(X-repmat(mu,1,n)),1);
L=@(nu,mu,sigma,w,delta)-2*log(gamma((d+nu)/2))+2*log(gamma(nu/2))-nu.*log(nu)+(d+nu).*sum(w.*log(nu+delta))+log(det(sigma));
mu=zeros(d,1);
w=ones(1,n)/n;
anz_steps=300;
niter=10000;
nus=[1,2,5,10,100];
for sigma_choices=1:4
    switch sigma_choices
        case 1
            sigma=.1*eye(d);
        case 2
            sigma=eye(d);
        case 3
            sigma=10*eye(d);
        case 4
            sigma=[2,-1;-1,2];
    end
    for nu=[1,2,5,10,100]
        disp(nu)
        disp(sigma)
        times_gmmf=zeros(1,niter);
        steps_gmmf=zeros(1,niter);
        
        times_mmf=zeros(1,niter);
        steps_mmf=zeros(1,niter);
        
        times_em=zeros(1,niter);
        steps_em=zeros(1,niter);
        
        times_aem=zeros(1,niter);
        steps_aem=zeros(1,niter);
        
        for iter=1:niter
            
            X=sample(n,mu,nu,sigma);
            
            % Estimating using the GMMF
            [~,~,~,num_steps_gmmf,time_gmmf,~]=iterate_studentT(X,w,'GMMF',anz_steps,1);
            
            % Estimating using the MMF
            [~,~,~,num_steps_mmf,time_mmf,~]=iterate_studentT(X,w,'MMF',anz_steps,1);
            
            % Estimating using the EM-Algorithm
            [~,~,~,num_steps_em,time_em,~]=iterate_studentT(X,w,'EM',anz_steps,1);
            
            % Estimating using the accelerated EM-Algorithm
            [~,~,~,num_steps_aem,time_aem,~]=iterate_studentT(X,w,'aEM',anz_steps,1);

            times_gmmf(iter)=time_gmmf;
            steps_gmmf(iter)=num_steps_gmmf;
            times_mmf(iter)=time_mmf;
            steps_mmf(iter)=num_steps_mmf;
            times_em(iter)=time_em;
            steps_em(iter)=num_steps_em;
            times_aem(iter)=time_aem;
            steps_aem(iter)=num_steps_aem;
            if mod(iter,10)==0
                disp(['step ' num2str(iter)])
            end
        end
        
        
        av_time_gmmf=mean(times_gmmf);
        std_time_gmmf=(var(times_gmmf))^.5;
        av_steps_gmmf=mean(steps_gmmf);
        std_steps_gmmf=(var(steps_gmmf))^.5;
        
        av_time_mmf=mean(times_mmf);
        std_time_mmf=(var(times_mmf))^.5;
        av_steps_mmf=mean(steps_mmf);
        std_steps_mmf=(var(steps_mmf))^.5;
        
        av_time_em=mean(times_em);
        std_time_em=(var(times_em))^.5;
        av_steps_em=mean(steps_em);
        std_steps_em=(var(steps_em))^.5;
        
        av_time_aem=mean(times_aem);
        std_time_aem=(var(times_aem))^.5;
        av_steps_aem=mean(steps_aem);
        std_steps_aem=(var(steps_aem))^.5;
        
        fileID=fopen('results','at');
        fprintf(fileID,'nu=%d sigma=[%d,%d;%d,%d]\n',nu,sigma(1,1),sigma(1,2),sigma(2,1),sigma(2,2));
        fprintf(fileID,'GMMF: Time: %d+-%d Steps: %d+-%d\n', av_time_gmmf,std_time_gmmf,av_steps_gmmf,std_steps_gmmf);
        fprintf(fileID,' MMF: Time: %d+-%d Steps: %d+-%d\n', av_time_mmf,std_time_mmf,av_steps_mmf,std_steps_mmf);
        fprintf(fileID,'  EM: Time: %d+-%d Steps: %d+-%d\n', av_time_em,std_time_em,av_steps_em,std_steps_em);
        fprintf(fileID,' aEM: Time: %d+-%d Steps: %d+-%d\n\n', av_time_aem,std_time_aem,av_steps_aem,std_steps_aem);
        fprintf(fileID,'Latex table steps: EM, aEM, MMF, GMMF:\n');
        fprintf(fileID,'$%d$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$\n',nu,av_steps_em,std_steps_em,av_steps_aem,std_steps_aem,av_steps_mmf,std_steps_mmf,av_steps_gmmf,std_steps_gmmf);
        fprintf(fileID,'Latex table time: EM, aEM, MMF, GMMF:\n');
        fprintf(fileID,'$%d$&$%1.6f\\pm%1.5f$&$%1.6f\\pm%1.5f$&$%1.6f\\pm%1.5f$&$%1.6f\\pm%1.5f$\n\n',nu,av_time_em,std_time_em,av_time_aem,std_time_aem,av_time_mmf,std_time_mmf,av_time_gmmf,std_time_gmmf);
        fclose(fileID);
    end
    
end
