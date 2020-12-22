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
% This function implements the algorithms
% - MMF
% - GMMF
% - EM algorithm
% - accelerated EM algorithm
% - ECME algorithm
% for Student-t distributions with DAAREM acceleration.
%
% INPUTS:
%     X             - d x n array containing the samples
%     w             - 1 x n array containing the weights of the samples
%     step_algorithm- Choice of the algorithm, which determines the nu, mu
%                     and sigma step. step_algorithm has to be of type
%                     string.
%                     For step_algorithm='MMF' the function performs the
%                     MMF.
%                     For step_algorithm='GMMF' the function performs the
%                     GMMF.
%                     For step_algorithm='aEM' the function performs the
%                     accelerated EM algorith,
%                     For step_algorithm='ECME' the function perfors the
%                     ECME algorithm.
%                     Alternatively the step_algorithm can be declared by a
%                     function handle.
%     anz_steps     - Maximum number of iterations
%     stop          - stop==1 means we apply the stopping criteria, 
%                     stop==0 means we perform exactly anz_steps spteps.
%     abs_criteria  - Replaces the relative stopping criteria by an
%                     absolute stopping criteria
%     regularize    - we add in each step regularize*eye(d) to sigma to
%                     prevent that sigma becomes singular.
%                     If X contains at least n affine independent samples,
%                     this parameter should be set to 0.
%     save_obj      - save_obj==1 means we save in each step the negative
%                     log-likelihood value. Otherwise save_obj should be
%                     set to 0. Note that save_obj==1 slows down the
%                     performance of the algorithm.
%
% OUTPUTS:
%     mu,nu,sigma   - estimated parameters
%     num_steps     - if stop==1 then num_steps contains the number of
%                     steps, until the stopping criteria is reached.
%                     Otherwise NaN is returned.
%     time          - execution time of the algorithm
%     objective     - if save_obj==1, objective(r) contains the negative
%                     log-likelihood value of nu_r,mu_r,Sigma_r.
%                     Otherwise objective(r)=0.
function [mu,nu,sigma,num_steps,time,objective]=daarem_studentT(X,w,step_algorithm,anz_steps,stop,abs_criteria,regularize,save_obj)
if nargin<8
    save_obj=0;
end
if nargin<7
    regularize=0;
end
if nargin<6
    abs_criteria=0;
end
if nargin<5
    stop=0;
end
if ~isa(step_algorithm,'function_handle')
    if isa(step_algorithm,'char')
        switch step_algorithm
            case 'GMMF'
                step_algorithm=@(X,w,nu_r,mu_r,sigma_r,regularize,delta_r)studentT_step(X,w,nu_r,mu_r,sigma_r,regularize,@nu_step_gmmf,delta_r);
            case 'MMF'
                step_algorithm=@(X,w,nu_r,mu_r,sigma_r,regularize,delta_r)studentT_step(X,w,nu_r,mu_r,sigma_r,regularize,@nu_step_mmf,delta_r);
            case 'EM'
                step_algorithm=@(X,w,nu_r,mu_r,sigma_r,regularize,delta_r)studentT_step_em(X,w,nu_r,mu_r,sigma_r,regularize,@nu_step_em,delta_r);
            case 'aEM'
                step_algorithm=@(X,w,nu_r,mu_r,sigma_r,regularize,delta_r)studentT_step(X,w,nu_r,mu_r,sigma_r,regularize,@nu_step_em,delta_r);
            case 'ECME'
                step_algorithm=@(X,w,nu_r,mu_r,sigma_r,regularize,delta_r)studentT_step_em(X,w,nu_r,mu_r,sigma_r,regularize,@nu_step_gmmf,delta_r);
            otherwise
                error('Unknown algorithm.')
        end
    else
        error('Unknown algorithm. step_algorithm has to be a function handle or a string.')
    end
end
tol=1e-5;
if abs_criteria==1
    tol=1e-8;
end
[d,n]=size(X);
L=@(nu,mu,sigma,w,delta)-2*gammaln((d+nu)/2)+2*gammaln(nu/2)-nu.*log(nu)+(d+nu).*sum(w.*log(nu+delta))+log(det(sigma));
objective=zeros(1,anz_steps+1);
delta=@(samples,mu,sigma)sum(((sigma^(-1))*(samples-repmat(mu,1,n))).*(samples-repmat(mu,1,n)),1);
mu_r=sum(X,2)/n;
sigma_r=(X-repmat(mu_r,1,n))*(X-repmat(mu_r,1,n))'/n+regularize*eye(d);
nu_r=3;
tic;
weiter=1;
r=0;

% DAAREM Parameters
epsilon=0.01;
epsilon_c=0;
alpha=1.2;
m=3;
kappa=25;
D=2*kappa;

% DAAREM initialization
c_k=1;
s_k=0;
delta_r=delta(X,mu_r,sigma_r);
l_star=L(nu_r,mu_r,sigma_r,w,delta_r);
Fk=[];
Xk=[];
lambda_k=100000;
r_k=0;
norm_step=0;
daarem_step=0;
while weiter==1
    r=r+1;
    if sum(delta_r<0)>0 || isnan(nu_r)
        disp('hm')
        error('Numerical errors')
    end
    if r==1
        l_xk=L(nu_r,mu_r,sigma_r,w,delta_r);
    end
    if save_obj==1
        objective(r)=l_xk;
    end
    
    % fixed point algorithm step
    [nu_r_plus_one,mu_r_plus_one,sigma_r_plus_one,delta_r_plus_one]=step_algorithm(X,w,nu_r,mu_r,sigma_r,regularize,delta_r);
    
    % DAAREM step
    m_k=min(m,c_k);
    if r~=1
        Xk=Xk(:,(end-m_k+1):end);
    end
    x_k_alt=[nu_r;mu_r;sigma_r(:)];
    fk_vec_neu=[nu_r_plus_one-nu_r;mu_r_plus_one-mu_r;sigma_r_plus_one(:)-sigma_r(:)];
    if r==1
        fk_vec=fk_vec_neu;
        x_k=x_k_alt+fk_vec;
        l_xk=l(x_k);
    else
        Fk=[Fk,fk_vec_neu-fk_vec];
        Fk=Fk(:,(end-m_k+1):end);
        s_k_old=s_k;
        c_k_old=c_k;
        l_star_old=l_star;
        fk_vec=fk_vec_neu;
        daarem_delta_k=1/(1+alpha^(kappa-s_k));
        
        [Uk,Dk,Vk]=svd(Fk,'econ');
        uf=Uk'*fk_vec;
        if m_k==1
            d_vec=Dk(1);
        else
            d_vec=diag(Dk);
        end
        d_vec_sq=d_vec.^2;
        [lambda_k,r_k]=findDAAREMsolution(lambda_k,r_k,sqrt(daarem_delta_k),uf,d_vec);
        
        dd=(d_vec.*uf)./(d_vec_sq+lambda_k);
        gamma_k=Vk'*dd;
        t_kp1=x_k_alt+fk_vec-(Xk+Fk)*gamma_k;
        if any(isnan(t_kp1))||any(~isreal(t_kp1))
            disp('baeh')
        end
        l_tkp1=l(t_kp1);
        if l_tkp1<=l_xk+epsilon
            x_k=t_kp1;
            s_new=s_k+1;
            daarem_step=daarem_step+1;
            l_xk=l_tkp1;
        else
            x_k=x_k_alt+fk_vec;
            s_new=s_k;
            norm_step=norm_step+1;
            l_xk=l(x_k);
        end
        if mod(r-1,m)==0 && l_xk<=l_star+epsilon_c
            c_k=1;
            l_star=l_xk;
        elseif mod(r-1,m)==0 && l_xk>l_star+epsilon_c 
            s_new=max(s_new-m,-D);
            c_k=1;
            l_star=l_xk;
        else
            c_k=c_k+1;
        end
        s_k=s_new;
    end
    nu_r_plus_one_x=x_k(1);
    mu_r_plus_one_x=x_k(2:d+1);
    sigma_r_plus_one_x=reshape(x_k(2+d:end),d,d);
    delta_r_plus_one_x=nan;
    if nu_r_plus_one_x<=0 || min(eig(sigma_r_plus_one_x))<1e-5 || sum(sum(abs(sigma_r_plus_one-sigma_r_plus_one)))>1e-5
        nu_r_plus_one_x=nu_r_plus_one;
        mu_r_plus_one_x=mu_r_plus_one;
        sigma_r_plus_one_x=sigma_r_plus_one;
        delta_r_plus_one_x=delta_r_plus_one;
        if save_obj==1
            l_xk=L(nu_r_plus_one_x,mu_r_plus_one_x,sigma_r_plus_one_x,w,delta(X,mu_r_plus_one_x,sigma_r_plus_one_x));
        end
        x_k=[nu_r_plus_one;mu_r_plus_one;sigma_r_plus_one(:)];
        s_k=s_k_old;
        c_k=c_k_old;
        l_star=l_star_old;
    end
    Xk=[Xk,x_k-x_k_alt];
    
    if stop==1 && abs_criteria==0
        eps=((sum((mu_r-mu_r_plus_one).^2)+sum(sum((sigma_r-sigma_r_plus_one).^2)))/(sum(mu_r.^2)+sum(sum(sigma_r.^2))))^.5+((log(nu_r_plus_one)-log(nu_r))^2/log(nu_r)^2)^.5;
        if r==10000
            weiter=0;
        end
        if eps<tol
            weiter=0;
        end
    end
    if stop==1 && abs_criteria==1
        eps=(sum((mu_r-mu_r_plus_one).^2)+sum(sum((sigma_r-sigma_r_plus_one).^2))+(log(nu_r_plus_one)-log(nu_r))^2)^.5;
        if r==10000
            weiter=0;
        end
        if eps<tol
            weiter=0;
        end
    end
    mu_r=mu_r_plus_one_x;
    sigma_r=sigma_r_plus_one_x;
    nu_r=nu_r_plus_one_x;
    delta_r=delta_r_plus_one_x;
    if stop==0
        if r==anz_steps
            weiter=0;
        end
    end
end
if save_obj==1
    objective(r+1)=L(nu_r,mu_r,sigma_r,w,delta(X,mu_r,sigma_r));
end
time=toc;
if stop==1
    num_steps=r;
else
    num_steps=nan;
end
mu=mu_r; nu=nu_r;sigma=sigma_r;

    function out=l(x)
        try
            nu_x=x(1);
            mu_x=x(2:1+d);
            sigma_x=reshape(x(2+d:end),d,d);
            if nu_x<0 || sum(sum(abs(sigma_x-sigma_x')))>1e-5  || min(eig(sigma_x))<1e-5
                out=inf;
                return
            end
            out=L(nu_x,mu_x,sigma_x,w,delta(X,mu_x,sigma_x));
        catch
            out=inf;
        end
    end

    function [lambda_new,r_new]=findDAAREMsolution(lambda_old,r_old,sqrt_daarem_delta,uf,d_vec)
        mk2=min(m_k,size(d_vec,1));
        beta_LS_tilde=uf./d_vec;
        d_vec_times_uf=d_vec.*uf;
        norm_beta_LS_tilde=norm(beta_LS_tilde);
        v_k=sqrt_daarem_delta*norm_beta_LS_tilde;
        lambda_new=lambda_old-r_old/v_k;
        phi=@(lambda)norm((d_vec_times_uf)./(d_vec_sq+lambda))-v_k;
        phi_der=@(lambda)-1/norm((d_vec_times_uf)./(d_vec_sq+lambda)) * sum((d_vec(1:mk2).*uf(1:mk2)).^2./(d_vec_sq(1:mk2)+lambda).^3);
        L_=-phi(0)/phi_der(0);
        U_=norm(Fk'*fk_vec)/v_k;
        l_stop=(1+alpha^(kappa-s_k+.5))^(-.5)*norm_beta_LS_tilde;
        u_stop=(1+alpha^(kappa-s_k-.5))^(-.5)*norm_beta_LS_tilde;
        
        for iter=1:10
            if lambda_new<= L_ || lambda_new >= U_
                lambda_new=max(0.001*U_,sqrt(L_*U_));
            end
            phi_val=phi(lambda_new);
            phi_der_val=phi_der(lambda_new);
            norm_s=norm((d_vec.*uf)./(d_vec.^2+lambda_new));
            r_new=norm_s*phi_val/phi_der_val;
            if norm_s >= l_stop && norm_s <= u_stop
                return
            end
            if phi_val<0
                U_=lambda_new;
            end
            L_=max(L_,lambda_new-phi_val/phi_der_val);
            lambda_new=lambda_new-(norm_s/v_k)*(phi_val/phi_der_val);
        end
    end
end

