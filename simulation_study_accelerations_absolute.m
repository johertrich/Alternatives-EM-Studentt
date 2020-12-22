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
% This script performs a simulation study as in Table 5 and 6  in the paper and
% saves the result in a text document called 'results_with_accelerations_absolute'

clear all
addpath('algorithms')
delete results_with_accelerations_absolute;
d=2;
n=1000;
delta=@(X,mu,sigma)sum(((sigma^(-1))*(X-repmat(mu,1,n))).*(X-repmat(mu,1,n)),1);
L=@(nu,mu,sigma,w,delta)-2*gammaln((d+nu)/2)+2*gammaln(nu/2)-nu.*log(nu)+(d+nu).*sum(w.*log(nu+delta))+log(det(sigma));
mu=zeros(d,1);
w=ones(1,n)/n;
anz_steps=300;
niter=1000;
sigma=.1*eye(d);
nus=[1,2,5,10,100];
avg_times=zeros(15,size(nus,2));
std_times=zeros(15,size(nus,2));
avg_steps=zeros(15,size(nus,2));
std_steps=zeros(15,size(nus,2));


for i=1:size(nus,2)
    nu=nus(i);
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
    
    times_ecme=zeros(1,niter);
    steps_ecme=zeros(1,niter);
    
    times_daarem_em=zeros(1,niter);
    steps_daarem_em=zeros(1,niter);
    
    times_daarem_ecme=zeros(1,niter);
    steps_daarem_ecme=zeros(1,niter);
    
    times_daarem_aem=zeros(1,niter);
    steps_daarem_aem=zeros(1,niter);
    
    times_daarem_mmf=zeros(1,niter);
    steps_daarem_mmf=zeros(1,niter);
    
    times_daarem_gmmf=zeros(1,niter);
    steps_daarem_gmmf=zeros(1,niter);
    
    times_squarem_gmmf=zeros(1,niter);
    steps_squarem_gmmf=zeros(1,niter);
    
    times_squarem_mmf=zeros(1,niter);
    steps_squarem_mmf=zeros(1,niter);
    
    times_squarem_em=zeros(1,niter);
    steps_squarem_em=zeros(1,niter);
    
    times_squarem_aem=zeros(1,niter);
    steps_squarem_aem=zeros(1,niter);
    
    times_squarem_ecme=zeros(1,niter);
    steps_squarem_ecme=zeros(1,niter);
    
    for iter=1:niter
        X=sample(n,mu,nu,sigma);
        Ld=@(nu,mu,sigma)L(nu,mu,sigma,w,delta(X,mu,sigma));
        
        % Estimating using the GMMF
        [~,~,~,num_steps_gmmf,time_gmmf,~]=iterate_studentT(X,w,'GMMF',anz_steps,1,1);
        
        % Estimating using the MMF
        [~,~,~,num_steps_mmf,time_mmf,~]=iterate_studentT(X,w,'MMF',anz_steps,1,1);
        
        % Estimating using the EM-Algorithm
        [~,~,~,num_steps_em,time_em,~]=iterate_studentT(X,w,'EM',anz_steps,1,1);
        
        % Estimating using the accelerated EM-Algorithm
        [~,~,~,num_steps_aem,time_aem,~]=iterate_studentT(X,w,'aEM',anz_steps,1,1);
        
        % Estimating using the ECME-Algorithm
        [~,~,~,num_steps_ecme,time_ecme,~]=iterate_studentT(X,w,'ECME',anz_steps,1,1);
        
        % Estimating using the GMMF with DAAREM
        [~,~,~,num_steps_daarem_gmmf,time_daarem_gmmf,~]=daarem_studentT(X,w,'GMMF',anz_steps,1,1);
        
        % Estimating using the MMF with DAAREM
        [~,~,~,num_steps_daarem_mmf,time_daarem_mmf,~]=daarem_studentT(X,w,'MMF',anz_steps,1,1);
        
        % Estimating using the EM-Algorithm with DAAREM
        [~,~,~,num_steps_daarem_em,time_daarem_em,~]=daarem_studentT(X,w,'EM',anz_steps,1,1);
        
        % Estimating using the accelerated EM-Algorithm with DAAREM
        [~,~,~,num_steps_daarem_aem,time_daarem_aem,~]=daarem_studentT(X,w,'aEM',anz_steps,1,1);
        
        % Estimating using the ECME-Algorithm with DAAREM
        [~,~,~,num_steps_daarem_ecme,time_daarem_ecme,~]=daarem_studentT(X,w,'ECME',anz_steps,1,1);
        
        % Estimating using the GMMF with SQUAREM
        [~,~,~,num_steps_squarem_gmmf,time_squarem_gmmf,~]=squarem_studentT(X,w,'GMMF',anz_steps,1,1);
        
        % Estimating using the MMF with SQUAREM
        [~,~,~,num_steps_squarem_mmf,time_squarem_mmf,~]=squarem_studentT(X,w,'MMF',anz_steps,1,1);
        
        % Estimating using the EM algorithm with SQUAREM
        [~,~,~,num_steps_squarem_em,time_squarem_em,~]=squarem_studentT(X,w,'EM',anz_steps,1,1);
        
        % Estimating using the accelerated EM algorithm with SQUAREM
        [~,~,~,num_steps_squarem_aem,time_squarem_aem,~]=squarem_studentT(X,w,'aEM',anz_steps,1,1);
        
        % Estimating using the ECME algorithm with SQUAREM
        [~,~,~,num_steps_squarem_ecme,time_squarem_ecme,~]=squarem_studentT(X,w,'ECME',anz_steps,1,1);
        
        times_gmmf(iter)=time_gmmf;
        steps_gmmf(iter)=num_steps_gmmf;
        times_mmf(iter)=time_mmf;
        steps_mmf(iter)=num_steps_mmf;
        times_em(iter)=time_em;
        steps_em(iter)=num_steps_em;
        times_aem(iter)=time_aem;
        steps_aem(iter)=num_steps_aem;
        times_ecme(iter)=time_ecme;
        steps_ecme(iter)=num_steps_ecme;
        times_daarem_em(iter)=time_daarem_em;
        steps_daarem_em(iter)=num_steps_daarem_em;
        times_daarem_ecme(iter)=time_daarem_ecme;
        steps_daarem_ecme(iter)=num_steps_daarem_ecme;
        times_daarem_aem(iter)=time_daarem_aem;
        steps_daarem_aem(iter)=num_steps_daarem_aem;
        times_daarem_mmf(iter)=time_daarem_mmf;
        steps_daarem_mmf(iter)=num_steps_daarem_mmf;
        times_daarem_gmmf(iter)=time_daarem_gmmf;
        steps_daarem_gmmf(iter)=num_steps_daarem_gmmf;
        times_squarem_gmmf(iter)=time_squarem_gmmf;
        steps_squarem_gmmf(iter)=num_steps_squarem_gmmf;
        times_squarem_mmf(iter)=time_squarem_mmf;
        steps_squarem_mmf(iter)=num_steps_squarem_mmf;
        times_squarem_em(iter)=time_squarem_em;
        steps_squarem_em(iter)=num_steps_squarem_em;
        times_squarem_aem(iter)=time_squarem_aem;
        steps_squarem_aem(iter)=num_steps_squarem_aem;
        times_squarem_ecme(iter)=time_squarem_ecme;
        steps_squarem_ecme(iter)=num_steps_squarem_ecme;
        
        
        if mod(iter,1)==0
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
    
    av_time_ecme=mean(times_ecme);
    std_time_ecme=(var(times_ecme))^.5;
    av_steps_ecme=mean(steps_ecme);
    std_steps_ecme=(var(steps_ecme))^.5;
    
    av_time_daarem_em=mean(times_daarem_em);
    std_time_daarem_em=(var(times_daarem_em))^.5;
    av_steps_daarem_em=mean(steps_daarem_em);
    std_steps_daarem_em=(var(steps_daarem_em))^.5;
    
    av_time_daarem_ecme=mean(times_daarem_ecme);
    std_time_daarem_ecme=(var(times_daarem_ecme))^.5;
    av_steps_daarem_ecme=mean(steps_daarem_ecme);
    std_steps_daarem_ecme=(var(steps_daarem_ecme))^.5;
    
    av_time_daarem_aem=mean(times_daarem_aem);
    std_time_daarem_aem=(var(times_daarem_aem))^.5;
    av_steps_daarem_aem=mean(steps_daarem_aem);
    std_steps_daarem_aem=(var(steps_daarem_aem))^.5;
    
    av_time_daarem_mmf=mean(times_daarem_mmf);
    std_time_daarem_mmf=(var(times_daarem_mmf))^.5;
    av_steps_daarem_mmf=mean(steps_daarem_mmf);
    std_steps_daarem_mmf=(var(steps_daarem_mmf))^.5;
    
    av_time_daarem_gmmf=mean(times_daarem_gmmf);
    std_time_daarem_gmmf=(var(times_daarem_gmmf))^.5;
    av_steps_daarem_gmmf=mean(steps_daarem_gmmf);
    std_steps_daarem_gmmf=(var(steps_daarem_gmmf))^.5;
    
    av_time_squarem_gmmf=mean(times_squarem_gmmf);
    std_time_squarem_gmmf=(var(times_squarem_gmmf))^.5;
    av_steps_squarem_gmmf=mean(steps_squarem_gmmf);
    std_steps_squarem_gmmf=(var(steps_squarem_gmmf))^.5;
    
    av_time_squarem_mmf=mean(times_squarem_mmf);
    std_time_squarem_mmf=(var(times_squarem_mmf))^.5;
    av_steps_squarem_mmf=mean(steps_squarem_mmf);
    std_steps_squarem_mmf=(var(steps_squarem_mmf))^.5;
    
    av_time_squarem_em=mean(times_squarem_em);
    std_time_squarem_em=(var(times_squarem_em))^.5;
    av_steps_squarem_em=mean(steps_squarem_em);
    std_steps_squarem_em=(var(steps_squarem_em))^.5;
    
    av_time_squarem_aem=mean(times_squarem_aem);
    std_time_squarem_aem=(var(times_squarem_aem))^.5;
    av_steps_squarem_aem=mean(steps_squarem_aem);
    std_steps_squarem_aem=(var(steps_squarem_aem))^.5;
    
    av_time_squarem_ecme=mean(times_squarem_ecme);
    std_time_squarem_ecme=(var(times_squarem_ecme))^.5;
    av_steps_squarem_ecme=mean(steps_squarem_ecme);
    std_steps_squarem_ecme=(var(steps_squarem_ecme))^.5;
    
    avg_times(:,i)=[av_time_em;av_time_aem;av_time_mmf;av_time_gmmf;av_time_ecme;...
        av_time_daarem_em;av_time_daarem_aem;av_time_daarem_mmf;av_time_daarem_gmmf;av_time_daarem_ecme;...
        av_time_squarem_em;av_time_squarem_aem;av_time_squarem_mmf;av_time_squarem_gmmf;av_time_squarem_ecme];
    std_times(:,i)=[std_time_em;std_time_aem;std_time_mmf;std_time_gmmf;std_time_ecme;...
        std_time_daarem_em;std_time_daarem_aem;std_time_daarem_mmf;std_time_daarem_gmmf;std_time_daarem_ecme;...
        std_time_squarem_em;std_time_squarem_aem;std_time_squarem_mmf;std_time_squarem_gmmf;std_time_squarem_ecme];
    avg_steps(:,i)=[av_steps_em;av_steps_aem;av_steps_mmf;av_steps_gmmf;av_steps_ecme;...
        av_steps_daarem_em;av_steps_daarem_aem;av_steps_daarem_mmf;av_steps_daarem_gmmf;av_steps_daarem_ecme;...
        av_steps_squarem_em;av_steps_squarem_aem;av_steps_squarem_mmf;av_steps_squarem_gmmf;av_steps_squarem_ecme];
    std_steps(:,i)=[std_steps_em;std_steps_aem;std_steps_mmf;std_steps_gmmf;std_steps_ecme;...
        std_steps_daarem_em;std_steps_daarem_aem;std_steps_daarem_mmf;std_steps_daarem_gmmf;std_steps_daarem_ecme;...
        std_steps_squarem_em;std_steps_squarem_aem;std_steps_squarem_mmf;std_steps_squarem_gmmf;std_steps_squarem_ecme];    
    
    fileID=fopen('results_with_accelerations_absolute','at');
    fprintf(fileID,'nu=%d sigma=[%d,%d;%d,%d]\n',nu,sigma(1,1),sigma(1,2),sigma(2,1),sigma(2,2));
    fprintf(fileID,'GMMF: Time: %d+-%d Steps: %d+-%d\n', av_time_gmmf,std_time_gmmf,av_steps_gmmf,std_steps_gmmf);
    fprintf(fileID,' MMF: Time: %d+-%d Steps: %d+-%d\n', av_time_mmf,std_time_mmf,av_steps_mmf,std_steps_mmf);
    fprintf(fileID,'  EM: Time: %d+-%d Steps: %d+-%d\n', av_time_em,std_time_em,av_steps_em,std_steps_em);
    fprintf(fileID,' aEM: Time: %d+-%d Steps: %d+-%d\n', av_time_aem,std_time_aem,av_steps_aem,std_steps_aem);
    fprintf(fileID,'ECME: Time: %d+-%d Steps: %d+-%d\n\n', av_time_ecme,std_time_ecme,av_steps_ecme,std_steps_ecme);
    fprintf(fileID,'With DAAREM acceleration:\n');
    fprintf(fileID,'GMMF: Time: %d+-%d Steps: %d+-%d\n', av_time_daarem_gmmf,std_time_daarem_gmmf,av_steps_daarem_gmmf,std_steps_daarem_gmmf);
    fprintf(fileID,' MMF: Time: %d+-%d Steps: %d+-%d\n', av_time_daarem_mmf,std_time_daarem_mmf,av_steps_daarem_mmf,std_steps_daarem_mmf);
    fprintf(fileID,'  EM: Time: %d+-%d Steps: %d+-%d\n', av_time_daarem_em,std_time_daarem_em,av_steps_daarem_em,std_steps_daarem_em);
    fprintf(fileID,' aEM: Time: %d+-%d Steps: %d+-%d\n', av_time_daarem_aem,std_time_daarem_aem,av_steps_daarem_aem,std_steps_daarem_aem);
    fprintf(fileID,'ECME: Time: %d+-%d Steps: %d+-%d\n\n', av_time_daarem_ecme,std_time_daarem_ecme,av_steps_daarem_ecme,std_steps_daarem_ecme);
    fprintf(fileID,'With SQUAREM acceleration:\n');
    fprintf(fileID,'GMMF: Time: %d+-%d Steps: %d+-%d\n', av_time_squarem_gmmf,std_time_squarem_gmmf,av_steps_squarem_gmmf,std_steps_squarem_gmmf);
    fprintf(fileID,' MMF: Time: %d+-%d Steps: %d+-%d\n', av_time_squarem_mmf,std_time_squarem_mmf,av_steps_squarem_mmf,std_steps_squarem_mmf);
    fprintf(fileID,'  EM: Time: %d+-%d Steps: %d+-%d\n', av_time_squarem_em,std_time_squarem_em,av_steps_squarem_em,std_steps_squarem_em);
    fprintf(fileID,' aEM: Time: %d+-%d Steps: %d+-%d\n', av_time_squarem_aem,std_time_squarem_aem,av_steps_squarem_aem,std_steps_squarem_aem);
    fprintf(fileID,'ECME: Time: %d+-%d Steps: %d+-%d\n\n', av_time_squarem_ecme,std_time_squarem_ecme,av_steps_squarem_ecme,std_steps_squarem_ecme);
    fclose(fileID);
end
if size(nus,2)==5
    fileID=fopen('results_with_accelerations_absolute','at');
    fprintf(fileID,'Latex table steps:\n\n');
    fprintf(fileID,'EM&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$\\\\\n',avg_steps(1,1),std_steps(1,1),avg_steps(1,2),std_steps(1,2),avg_steps(1,3),std_steps(1,3),avg_steps(1,4),std_steps(1,4),avg_steps(1,5),std_steps(1,5));
    fprintf(fileID,'aEM&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$\\\\\n',avg_steps(2,1),std_steps(2,1),avg_steps(2,2),std_steps(2,2),avg_steps(2,3),std_steps(2,3),avg_steps(2,4),std_steps(2,4),avg_steps(2,5),std_steps(2,5));
    fprintf(fileID,'MMF&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$\\\\\n',avg_steps(3,1),std_steps(3,1),avg_steps(3,2),std_steps(3,2),avg_steps(3,3),std_steps(3,3),avg_steps(3,4),std_steps(3,4),avg_steps(3,5),std_steps(3,5));
    fprintf(fileID,'GMMF&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$\\\\\n',avg_steps(4,1),std_steps(4,1),avg_steps(4,2),std_steps(4,2),avg_steps(4,3),std_steps(4,3),avg_steps(4,4),std_steps(4,4),avg_steps(4,5),std_steps(4,5));
    fprintf(fileID,'ECME&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$\\\\\n',avg_steps(5,1),std_steps(5,1),avg_steps(5,2),std_steps(5,2),avg_steps(5,3),std_steps(5,3),avg_steps(5,4),std_steps(5,4),avg_steps(5,5),std_steps(5,5));
    fprintf(fileID,'DAAREM-EM&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$\\\\\n',avg_steps(6,1),std_steps(6,1),avg_steps(6,2),std_steps(6,2),avg_steps(6,3),std_steps(6,3),avg_steps(6,4),std_steps(6,4),avg_steps(6,5),std_steps(6,5));
    fprintf(fileID,'DAAREM-aEM&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$\\\\\n',avg_steps(7,1),std_steps(7,1),avg_steps(7,2),std_steps(7,2),avg_steps(7,3),std_steps(7,3),avg_steps(7,4),std_steps(7,4),avg_steps(7,5),std_steps(7,5));
    fprintf(fileID,'DAAREM-MMF&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$\\\\\n',avg_steps(8,1),std_steps(8,1),avg_steps(8,2),std_steps(8,2),avg_steps(8,3),std_steps(8,3),avg_steps(8,4),std_steps(8,4),avg_steps(8,5),std_steps(8,5));
    fprintf(fileID,'DAAREM-GMMF&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$\\\\\n',avg_steps(9,1),std_steps(9,1),avg_steps(9,2),std_steps(9,2),avg_steps(9,3),std_steps(9,3),avg_steps(9,4),std_steps(9,4),avg_steps(9,5),std_steps(9,5));
    fprintf(fileID,'DAAREM-ECME&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$\\\\\n',avg_steps(10,1),std_steps(10,1),avg_steps(10,2),std_steps(10,2),avg_steps(10,3),std_steps(10,3),avg_steps(10,4),std_steps(10,4),avg_steps(10,5),std_steps(10,5));
    fprintf(fileID,'SQUAREM-EM&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$\\\\\n',avg_steps(11,1),std_steps(11,1),avg_steps(11,2),std_steps(11,2),avg_steps(11,3),std_steps(11,3),avg_steps(11,4),std_steps(11,4),avg_steps(11,5),std_steps(11,5));
    fprintf(fileID,'SQUAREM-aEM&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$\\\\\n',avg_steps(12,1),std_steps(12,1),avg_steps(12,2),std_steps(12,2),avg_steps(12,3),std_steps(12,3),avg_steps(12,4),std_steps(12,4),avg_steps(12,5),std_steps(12,5));
    fprintf(fileID,'SQUAREM-MMF&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$\\\\\n',avg_steps(13,1),std_steps(13,1),avg_steps(13,2),std_steps(13,2),avg_steps(13,3),std_steps(13,3),avg_steps(13,4),std_steps(13,4),avg_steps(13,5),std_steps(13,5));
    fprintf(fileID,'SQUAREM-GMMF&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$\\\\\n',avg_steps(14,1),std_steps(14,1),avg_steps(14,2),std_steps(14,2),avg_steps(14,3),std_steps(14,3),avg_steps(14,4),std_steps(14,4),avg_steps(14,5),std_steps(14,5));
    fprintf(fileID,'SQUAREM-ECME&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$&$%4.2f\\pm%4.2f$\n\n',avg_steps(15,1),std_steps(15,1),avg_steps(15,2),std_steps(15,2),avg_steps(15,3),std_steps(15,3),avg_steps(15,4),std_steps(15,4),avg_steps(15,5),std_steps(15,5));
    fprintf(fileID,'Latex table times:\n\n');
fprintf(fileID,'EM&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$\\\\\n',avg_times(1,1),std_times(1,1),avg_times(1,2),std_times(1,2),avg_times(1,3),std_times(1,3),avg_times(1,4),std_times(1,4),avg_times(1,5),std_times(1,5));
    fprintf(fileID,'aEM&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$\\\\\n',avg_times(2,1),std_times(2,1),avg_times(2,2),std_times(2,2),avg_times(2,3),std_times(2,3),avg_times(2,4),std_times(2,4),avg_times(2,5),std_times(2,5));
    fprintf(fileID,'MMF&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$\\\\\n',avg_times(3,1),std_times(3,1),avg_times(3,2),std_times(3,2),avg_times(3,3),std_times(3,3),avg_times(3,4),std_times(3,4),avg_times(3,5),std_times(3,5));
    fprintf(fileID,'GMMF&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$\\\\\n',avg_times(4,1),std_times(4,1),avg_times(4,2),std_times(4,2),avg_times(4,3),std_times(4,3),avg_times(4,4),std_times(4,4),avg_times(4,5),std_times(4,5));
    fprintf(fileID,'ECME&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$\\\\\n',avg_times(5,1),std_times(5,1),avg_times(5,2),std_times(5,2),avg_times(5,3),std_times(5,3),avg_times(5,4),std_times(5,4),avg_times(5,5),std_times(5,5));
    fprintf(fileID,'DAAREM-EM&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$\\\\\n',avg_times(6,1),std_times(6,1),avg_times(6,2),std_times(6,2),avg_times(6,3),std_times(6,3),avg_times(6,4),std_times(6,4),avg_times(6,5),std_times(6,5));
    fprintf(fileID,'DAAREM-aEM&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$\\\\\n',avg_times(7,1),std_times(7,1),avg_times(7,2),std_times(7,2),avg_times(7,3),std_times(7,3),avg_times(7,4),std_times(7,4),avg_times(7,5),std_times(7,5));
    fprintf(fileID,'DAAREM-MMF&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$\\\\\n',avg_times(8,1),std_times(8,1),avg_times(8,2),std_times(8,2),avg_times(8,3),std_times(8,3),avg_times(8,4),std_times(8,4),avg_times(8,5),std_times(8,5));
    fprintf(fileID,'DAAREM-GMMF&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$\\\\\n',avg_times(9,1),std_times(9,1),avg_times(9,2),std_times(9,2),avg_times(9,3),std_times(9,3),avg_times(9,4),std_times(9,4),avg_times(9,5),std_times(9,5));
    fprintf(fileID,'DAAREM-ECME&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$\\\\\n',avg_times(10,1),std_times(10,1),avg_times(10,2),std_times(10,2),avg_times(10,3),std_times(10,3),avg_times(10,4),std_times(10,4),avg_times(10,5),std_times(10,5));
    fprintf(fileID,'SQUAREM-EM&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$\\\\\n',avg_times(11,1),std_times(11,1),avg_times(11,2),std_times(11,2),avg_times(11,3),std_times(11,3),avg_times(11,4),std_times(11,4),avg_times(11,5),std_times(11,5));
    fprintf(fileID,'SQUAREM-aEM&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$\\\\\n',avg_times(12,1),std_times(12,1),avg_times(12,2),std_times(12,2),avg_times(12,3),std_times(12,3),avg_times(12,4),std_times(12,4),avg_times(12,5),std_times(12,5));
    fprintf(fileID,'SQUAREM-MMF&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$\\\\\n',avg_times(13,1),std_times(13,1),avg_times(13,2),std_times(13,2),avg_times(13,3),std_times(13,3),avg_times(13,4),std_times(13,4),avg_times(13,5),std_times(13,5));
    fprintf(fileID,'SQUAREM-GMMF&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$\\\\\n',avg_times(14,1),std_times(14,1),avg_times(14,2),std_times(14,2),avg_times(14,3),std_times(14,3),avg_times(14,4),std_times(14,4),avg_times(14,5),std_times(14,5));
    fprintf(fileID,'SQUAREM-ECME&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$&$%1.5f\\pm%1.5f$\n\n',avg_times(15,1),std_times(15,1),avg_times(15,2),std_times(15,2),avg_times(15,3),std_times(15,3),avg_times(15,4),std_times(15,4),avg_times(15,5),std_times(15,5));
end
