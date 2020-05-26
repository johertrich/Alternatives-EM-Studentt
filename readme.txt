This code belongs to the paper [1]. Please cite [1], if you use this code.
It contains an implementation of the accelerations of the EM algorithm from [1]. Further,
the SQUAREM [4] and DAAREM [2] acceleration are applied to the algorithms and an application for noise
level estimation is provided. For more details see [1].
For the noise level estimation, we use functions to detect homogenous areas. The code of these functions is part
of the implementation of [3]. The original code for the homogenous area detection is availible at 
https://github.com/csutour/RNLF. Note that the functions might be modified. The code for the homogenous area
detection can be found in the directory hom_area_detection and is protected by the CeCILL Licence, which is also
located in the directory hom_area_detection.
For questions and bug reports, please contact Johannes Hertrich (j.hertrich(at)math.tu-berlin.de).

CONTENTS:

1. REQUIREMENTS
2. USAGE
3. CLASSES AND FUNCTIONS
4. EXAMPLES
5. REFERENCES

--------------- 1. REQUIREMENTS ------------------------

This code was written and testet with Matlab R2019a. For some of the examples the Matlab coder and gcc are required.

--------------- 2. USAGE -------------------------------

To compute the ML-estimator of samples X of size d in a dxn matrix, weights w in a 1xn matrix and a method
out of 'EM', 'aEM', 'MMF', 'GMMF' run the command

iterate_studentT(X,w,method,anz_steps);

Here anz_steps is the maximal number of steps. To incorporate the SQUAREM acceleration use

squarem_studentT(X,w,method,anz_steps);

To incorporate the DAAREM acceleration use

daarem_studentT(X,w,method,anz_steps);

A detailed description of all functions and their corresponding input and output parameters can be found in 
Section 3.

--------------- 3. FUNCTIONS ---------------------------

In the following we describe all functions in the directory algorithms and specify the corresponding input
and output parameters. The same description can be found in the headers of the functions in the code.

- daarem_studentT(X,w,step_algorithm,anz_steps,stop,abs_criteria,regularize,save_obj):
This function implements the algorithms
- MMF
- GMMF
- EM algorithm
- accelerated EM algorithm
- ECME algorithm
for Student-t distributions with DAAREM acceleration.
Inputs:
    X             - d x n array containing the samples
    w             - 1 x n array containing the weights of the samples
    step_algorithm- Choice of the algorithm, which determines the nu, mu
                    and sigma step. step_algorithm has to be of type
                    string.
                    For step_algorithm='MMF' the function performs the
                    MMF.
                    For step_algorithm='GMMF' the function performs the
                    GMMF.
                    For step_algorithm='aEM' the function performs the
                    accelerated EM algorith,
                    For step_algorithm='ECME' the function perfors the
                    ECME algorithm.
                    Alternatively the step_algorithm can be declared by a
                    function handle.
    anz_steps     - Maximum number of iterations
    stop          - stop==1 means we apply the stopping criteria, 
                    stop==0 means we perform exactly anz_steps spteps.
                    If nargin<5, stop is set to 0.
    abs_criteria  - Replaces the relative stopping criteria by an
                    absolute stopping criteria
                    If nargin<6, abs_criteria is set to 0.
    regularize    - we add in each step regularize*eye(d) to sigma to
                    prevent that sigma becomes singular.
                    If X contains at least n affine independent samples,
                    this parameter should be set to 0.
                    If nargin<7, regularize is set to 0.
    save_obj      - save_obj==1 means we save in each step the negative
                    log-likelihood value. Otherwise save_obj should be
                    set to 0. Note that save_obj==1 slows down the
                    performance of the algorithm.
                    If nargin<8, stop is set to 0.
Outputs [mu,nu,sigma,num_steps,time,objective]:
    mu,nu,sigma   - estimated parameters
    num_steps     - if stop==1 then num_steps contains the number of
                    steps, until the stopping criteria is reached.
                    Otherwise NaN is returned.
    time          - execution time of the algorithm
    objective     - if save_obj==1, objective(r) contains the negative
                    log-likelihood value of nu_r,mu_r,Sigma_r.
                    Otherwise objective(r)=0.



- iterate_studentT(X,w,step_algorithm,anz_steps,stop,abs_criteria,regularize,save_obj):
This function implements the algorithms
- MMF
- GMMF
- EM algorithm
- accelerated EM algorithm
- ECME algorithm
for Student-t distributions.
Inputs:
    X             - d x n array containing the samples
    w             - 1 x n array containing the weights of the samples
    step_algorithm- Choice of the algorithm, which determines the nu, mu
                    and sigma step. step_algorithm has to be of type
                    string.
                    For step_algorithm='MMF' the function performs the
                    MMF.
                    For step_algorithm='GMMF' the function performs the
                    GMMF.
                    For step_algorithm='aEM' the function performs the
                    accelerated EM algorith,
                    For step_algorithm='ECME' the function perfors the
                    ECME algorithm.
                    Alternatively the step_algorithm can be declared by a
                    function handle.
    anz_steps     - Maximum number of iterations
    stop          - stop==1 means we apply the stopping criteria, 
                    stop==0 means we perform exactly anz_steps spteps.
                    If nargin<5, stop is set to 0.
    abs_criteria  - Replaces the relative stopping criteria by an
                    absolute stopping criteria
                    If nargin<6, abs_criteria is set to 0.
    regularize    - we add in each step regularize*eye(d) to sigma to
                    prevent that sigma becomes singular.
                    If X contains at least n affine independent samples,
                    this parameter should be set to 0.
                    If nargin<7, regularize is set to 0.
    save_obj      - save_obj==1 means we save in each step the negative
                    log-likelihood value. Otherwise save_obj should be
                    set to 0. Note that save_obj==1 slows down the
                    performance of the algorithm.
                    If nargin<8, stop is set to 0.
Outputs [mu,nu,sigma,num_steps,time,objective]:
    mu,nu,sigma   - estimated parameters
    num_steps     - if stop==1 then num_steps contains the number of
                    steps, until the stopping criteria is reached.
                    Otherwise NaN is returned.
    time          - execution time of the algorithm
    objective     - if save_obj==1, objective(r) contains the negative
                    log-likelihood value of nu_r,mu_r,Sigma_r.
                    Otherwise objective(r)=0.



- newton(start,f,der_f,max_steps,tol)
This function implements newtons method for determining the zeros in the nu-steps.
Inputs:
    start       - initial value
    f           - function handle of the function
    der_f       - function handle of the derivative
    max_steps   - maximum number of steps
                  If nargin<4, max_steps is set to 1000
    tol         - stopping criteria. Stop if abs(f(val))<tol.
                  If nargin<5 tol is set to 1e-5
Outputs [zero,evals]:
    zero        - root of f
    evals       - number of evaluations of f and its derivative.



- nu_step_em(nu_r,delta_r,d,n):
This function implements the nu-step of the EM algorithm and the
accelerated EM algorithm
Inputs:
      nu_r    - current nu
      delta_r - 1 x n array, where delta_r(i)=(x_i-mu_r)'*sigma_r^(-1)*(x_i-mu_r)
      d       - dimension of the data
      n       - number of samples
Outputs [nu_r_plus_one,evals]:
      nu_r_plus_one - next nu
      evals         - number of evaluations of the subfunctions A and B



- nu_step_gmmf(nu_r,delta_r,d,n):
This function implements the nu-step of the GMMF
Inputs:
      nu_r    - current nu
      delta_r - 1 x n array, where delta_r(i)=(x_i-mu_r)'*sigma_r^(-1)*(x_i-mu_r)
      d       - dimension of the data
      n       - number of samples
Outputs [nu_r_plus_one,evals]:
      nu_r_plus_one - next nu
      evals         - number of evaluations of the subfunctions A and B



- nu_step_mmf(nu_r,delta_r,d,n):
This function implements the nu-step of the MMF
Inputs:
      nu_r    - current nu
      delta_r - 1 x n array, where delta_r(i)=(x_i-mu_r)'*sigma_r^(-1)*(x_i-mu_r)
      d       - dimension of the data
      n       - number of samples
Outputs [nu_r_plus_one,evals]:
      nu_r_plus_one - next nu
      evals         - number of evaluations of the subfunctions A and B



- sample(n,mu,nu,sigma):
This function generates samples of a multivariate Student-t distribution.
Inputs:
    n           - int. Number of samples generated.
    mu          - dx1 matrix. Location of the samples.
    nu          - double. Degree of freedom of the samples.
    sigma       - dxd matrix. Scatter matrix
Outputs:
    X           - dxn matrix. Samples.



- squarem_studentT(X,w,step_algorithm,anz_steps,stop,abs_criteria,regularize,save_obj):
This function implements the algorithms
- MMF
- GMMF
- EM algorithm
- accelerated EM algorithm
- ECME algorithm
for Student-t distributions with SQUAREM acceleration.
Inputs:
    X             - d x n array containing the samples
    w             - 1 x n array containing the weights of the samples
    step_algorithm- Choice of the algorithm, which determines the nu, mu
                    and sigma step. step_algorithm has to be of type
                    string.
                    For step_algorithm='MMF' the function performs the
                    MMF.
                    For step_algorithm='GMMF' the function performs the
                    GMMF.
                    For step_algorithm='aEM' the function performs the
                    accelerated EM algorith,
                    For step_algorithm='ECME' the function perfors the
                    ECME algorithm.
                    Alternatively the step_algorithm can be declared by a
                    function handle.
    anz_steps     - Maximum number of iterations
    stop          - stop==1 means we apply the stopping criteria, 
                    stop==0 means we perform exactly anz_steps spteps.
                    If nargin<5, stop is set to 0.
    abs_criteria  - Replaces the relative stopping criteria by an
                    absolute stopping criteria
                    If nargin<6, abs_criteria is set to 0.
    regularize    - we add in each step regularize*eye(d) to sigma to
                    prevent that sigma becomes singular.
                    If X contains at least n affine independent samples,
                    this parameter should be set to 0.
                    If nargin<7, regularize is set to 0.
    save_obj      - save_obj==1 means we save in each step the negative
                    log-likelihood value. Otherwise save_obj should be
                    set to 0. Note that save_obj==1 slows down the
                    performance of the algorithm.
                    If nargin<8, stop is set to 0.
Outputs [mu,nu,sigma,num_steps,time,objective]:
    mu,nu,sigma   - estimated parameters
    num_steps     - if stop==1 then num_steps contains the number of
                    steps, until the stopping criteria is reached.
                    Otherwise NaN is returned.
    time          - execution time of the algorithm
    objective     - if save_obj==1, objective(r) contains the negative
                    log-likelihood value of nu_r,mu_r,Sigma_r.
                    Otherwise objective(r)=0.



- studentT_step(X,w,nu_r,mu_r,sigma_r,regularize,nu_step,delta_r)
This function performs one step of the MMF, GMMF or accelerated EM algorithm for
estimating the parameters of a Student-t distribution.
Inputs:
    X                 - d x n array containing the samples
    w                 - 1 x n array containing the weights of the samples
    nu_r,mu_r,sigma_r - current estimates of nu, mu and sigma
    regularize        - we add in each step regularize*eye(d) to sigma to
                        prevent that sigma becomes singular.
                        If X contains at least n affine independent samples,
                        this parameter should be set to 0.
    nu_step           - declares the nu-step of the algorithm.
                        Choose nu_step=@nu_step_em for the accelerated EM
                        algorithm.
                        Choose nu_step=@nu_step_gmmf for the GMMF.
                        Choose nu_step=@nu_step_mmf for the MMF.
    delta_r           - 1 x n array containing the values of delta_r. Set
                        delta_r=nan to compute the values within this
                        function.
                        If nargin<7 delta_r is computed wighin this function
Outputs [nu_r_plus_one,mu_r_plus_one,sigma_r_plus_one,delta_r_plus_one]:
    nu_r_plus_one     - updated estimate of nu.
    mu_r_plus_one     - updated estimate of mu.
    sigma_r_plus_one  - updated estimate of sigma.
    delta_r_plus_one  - updated values of delta.


- studentT_step_em(X,w,nu_r,mu_r,sigma_r,regularize,nu_step,delta_r)
This function performs one step of the EM algortihm or the ECME algorithm for
estimating the parameters of a Student-t distribution.
Inputs:
    X                 - d x n array containing the samples
    w                 - 1 x n array containing the weights of the samples
    nu_r,mu_r,sigma_r - current estimates of nu, mu and sigma
    regularize        - we add in each step regularize*eye(d) to sigma to
                        prevent that sigma becomes singular.
                        If X contains at least n affine independent samples,
                        this parameter should be set to 0.
    nu_step           - declares the nu-step of the algorithm.
                        Choose nu_step=@nu_step_em for the EM algorithm.
                        Choose nu_step=@nu_step_gmmf for the ECME algorithm.
    delta_r           - 1 x n array containing the values of delta_r. Set
                        delta_r=nan to compute the values within this
                        function.
                        If nargin<7 delta_r is computed wighin this function
Outputs [nu_r_plus_one,mu_r_plus_one,sigma_r_plus_one,delta_r_plus_one]:
    nu_r_plus_one     - updated estimate of nu.
    mu_r_plus_one     - updated estimate of mu.
    sigma_r_plus_one  - updated estimate of sigma.
    delta_r_plus_one  - returns NaN.


--------------- 4. EXAMPLES ----------------------------

In the following we give a short description of the scripts in the main directory. 
Note that you can find a detailed description in [1].


- historgrams_nu:
This script draws n=1000 samples from a student t distribution for a fixed nu out of [1,2,5,10,100,200],
with identity scatter matrix and zero mean in two dimensions. Then it uses the GMMF to compute the ML-Estimator. 
This procedure is repeated 1000 times and the histograms of the outcoming nu estimates is plotted.
The resulting histograms are given in Figure 3 in [1].


- Likelihood_plot:
The script draws n=1000 samples from a student t distribution with fixed nu out of [1,2,5,10,100,200] and
computes the ML-estimator using the EM, aEM, MMF and GMMF algorithm. Then it generates plots of 
the Likelihood value versus the number of steps of the corresponding algorithms.
Some exemplary resulting plots are given in Figure 2 in [1].


- nu_estimation_camera:
This script implements the noise level estimation as described in Section 6.3 in [1].
First the script detects homogenous areas in the noisy input image using kendalls tau test. For this
we use the code of [3], see the header of the readme. On the homogenous areas we estimate the 
scatter and degree of freedom of the noise on the image. Then histograms of results and the homogenous areas
are plotted.
The results are given in Figure 4 in [1].


- nu_estimation_muehle:
Same as nu_estimation_camera with another image.
The results are given in Figure 5 in [1].


- simulation_study:
Runtime comparison of the EM, aEM, MMF and GMMF algorithm using different choices of nu and sigma.
The script generates for fixed nu out of [1,2,5,10,100] and Sigma out of [0.1 Id, Id,10 Id, [2,-1;-1,2]] 
n=1000 samples- Then it runs the algorithm and records execution time and number of steps of the algorithms using
a relative stopping criteria, see [1]. This procedure is repeated 10000 times and the mean and standard
deviation of the execution times and number of steps are computed. The results are printed into the file 
'results'. 


- simulation_study_accelerations:
Runtime comparison of the EM, aEM, MMF, GMMF and ECME algorithms with no acceleration, SQUAREM acceleration and
DAAREM acceleration.
The script generates for nu out of [1,2,5,10,100] n=1000 samples of a student-t distribution with zero mean and 
0.1 Id scatter. Then it runs the algorithms and records execution time and number of steps of the algorithms 
using a relative stopping criteria, see [1]. This procedure is repeated 1000 times and the mean and standard
deviation of the execution times and number of steps are computed. The results are printed into the file 
'results_with_accelerations'.
Exemplary results are given in Table 3 and 4 in [1].


- simulation_study_accelerations_absolute:
Similar as simulations_study_acceleration with the only difference, that an absolute stopping criterion is used
instead of a relative one. The results are printed into the file  'results_with_accelerations_absolute'.
Exemplary results are given in Table 5 and 6 in [1].


--------------- 5. REFERENCES --------------------------

[1] M. Hasannasab, J. Hertrich, F. Laus, and G. Steidl. 
Alternatives to the EM algorithm for ML-estimation of location, scatter matrix and degree of freedom 
of the student-t distribution.
ArXiv preprint arXiv:1910.06623, 2019.

[2] N. C. Henderson and R. Varadhan. 
Damped Anderson acceleration with restarts and monotonicity control for accelerating EM and EM-like algorithms.
Journal of Computational and Graphical Statistics, 28(4):834–846, 2019

[3] C. Sutour, C.A. Deledalle and J.F. Aujol.
Estimation of the noise level function based on a nonparametric detection of homogeneous image regions. 
SIAM Journal on Imaging Sciences, 8(4):2622–2661, 2015.

[4] R. Varadhan and C. Roland. 
Simple and globally convergent methods for accelerating the convergence of any EM algorithm. 
Scandinavian Journal of Statistics. Theory and Applications, 35(2):335–353, 2008.

