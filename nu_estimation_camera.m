% This code belongs to the paper
%
% M. Hasannasab, J. Hertrich, F. Laus, and G. Steidl. 
% Alternatives to the EM algorithm for ML-estimation of location, scatter
% matrix and degree of freedom of the student-t distribution.
% ArXiv preprint arXiv:1910.06623, 2019.
%
% If you use this code, please cite the paper.
%
% This script generates plots and diagrams as in Figure 4 in the Paper
close all;
addpath('hom_area_detection')
addpath('algorithms')

% Load image
img = double(imread('images/cameraman.pgm'));
for nu=[1,5]
    rng(6)
    noise=10*randn(size(img))./sqrt(chi2rnd(nu,size(img))/nu);
    
    img_nse = img + noise;
    
    param.W = 16;
    [noise, noise_info] = noise_estimation(img_nse, param);
    M=[];
    S=[];
    N = [];
    
    
    for l_block_start=1:noise_info.W:size(img_nse,1)-noise_info.W+1
        for k_block_start=1:noise_info.W:size(img_nse,2)-noise_info.W+1
            bloc = noise_info.hom(l_block_start:l_block_start-1+noise_info.W, k_block_start:k_block_start-1+noise_info.W);
            if sum(abs(bloc(:))) == noise_info.W*noise_info.W
                A = reshape(img_nse(l_block_start:l_block_start-1+noise_info.W, k_block_start:k_block_start-1+noise_info.W), [1,1, length(bloc(:))]);
                X=squeeze(A)';
                [mu_r,Sigma_r,nu_r]=iterate_studentT(X,size(A,2),'GMMF',300,1);
                S = [S; Sigma_r] ;
                N = [N; nu_r];
                M = [M; mu_r];
            end
            
        end
    end
    
    
    % Show results
    figure
    plothomogeneous(img_nse, noise_info.W, noise_info.hom, img)
    title(['Homogeneous areas. nu=' num2str(nu)])
    figure
    hist(N(:)), axis square
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',[.5 .5 .5],'EdgeColor','k')
    title(['Histogram nu. nu=' num2str(nu)])
    figure
    hist(S(:)), axis square
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',[.5 .5 .5],'EdgeColor','k')
    title(['Histogram sigma. nu=' num2str(nu)])
    
    disp(char(['Estimated nu = ', num2str(mean(N))]))
    disp(char(['geometric nu = ', num2str(geomean(N))]))
    disp(char(['Estimated sigma = ', num2str(mean(S))]))
    disp(char(['geometric sigma = ', num2str(geomean(S))]))
    
end
