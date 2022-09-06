addpath(genpath(cd))
clear
clc
%% parameter seting
opts.mu = 1e-6;
opts.rho = 1.1;
opts.max_iter = 500;
opts.DEBUG = 1;

%% Dataset setting
img_name = 'Indian_pines';
T = importdata(['./data-HSI/' img_name '_noised.mat']);
T_real = importdata(['./data-HSI/' img_name '_corrected.mat']);
gt = importdata(['./data-HSI/' img_name '_gt.mat']);
CTrain = [6 144 84 24 50 75 3 49 2 97 247 62 22 130 38 10];
lambda = 5e-3;

% img_name = 'Salinas';
% T = importdata(['.\data-HSI\' img_name '_corrected.mat']);
% T_real = importdata(['.\data-HSI\' img_name '_corrected.mat']);
% gt = importdata(['.\data-HSI\' img_name '_gt.mat']);
% CTrain = [200 373 198 140 268 396 358 1127 620 328 107 193 92 107 727 180];
% lambda = 5e-3;

% img_name = 'SalinasA';
% T = importdata(['./data-HSI/' img_name '_noised.mat']);
% T_real = importdata(['./data-HSI/' img_name '_corrected.mat']);
% gt = importdata(['./data-HSI/' img_name '_gt.mat']);
% CTrain = [39 0 0 0 0 0 0 0 0 134 62 153 67 80 ];
% lambda = 1e-2;

% img_name = 'BotswanaA';
% % T = importdata(['.\data-HSI\' img_name '_noised.mat']);
% T = importdata(['.\data-HSI\' img_name '_noised.mat']);
% T_real = importdata(['.\data-HSI\' img_name '_corrected.mat']);
% gt = importdata(['.\data-HSI\' img_name '_gt.mat']);
% CTrain = [9 10 0 0 0 0 26 0 31 25 31 18 27 5];
% lambda = 5e-1;

% img_name = 'PaviaUA';
% T = importdata(['.\data-HSI\' img_name '_noised.mat']);
% T_real = importdata(['.\data-HSI\' img_name '.mat']);
% gt = importdata(['.\data-HSI\' img_name '_gt.mat']);
% CTrain = [228 267 210 108 70 503 133 180 86];

%lambda = 1/sqrt(n3*max(n1,n2));
%lambda = 1/sqrt(max(n3,n1*n2));

[n1,n2,n3] = size(T);

%% algorithms setting
opts.fun = 'oddfunc' ; opts.gamma = 0.01;
% opts.fun = 'etp' ;  opts.gamma = 10;
% opts.fun = 'laplace' ; opts.gamma = 0.1;
% opts.fun = 'exponential' ; opts.gamma = 0.05 ;
% opts.fun = 'lp' ;        opts.gamma = 0.5 ;
% opts.fun = 'scad' ;      opts.gamma = 100 ;
% opts.fun = 'cappedl1' ; opts.gamma = 1000 ;
% opts.fun = 'mcp' ; opts.gamma = 10 ;
% opts.fun = 'geman' ;  opts.gamma = 10 ;
% opts.fun = 'logarithm' ; opts.gamma = 10 ;

tic
[TRec,  S] = trpca_explr_plus(T_real, lambda,opts);
time = toc;

%% classification
[OA_SVM1, OA_SVM2, AA_SVM1, AA_SVM2, Kappa_SVM1, Kappa_SVM2,ave_TPR_SVM1,ave_TPR_SVM2] = Classification_V2(T_real,TRec,gt,20,CTrain);

% %% save result
% file = fopen('classification_result.txt','w');
% fprintf(file,'  %s  %s  %s  %s  %s  %s  %s \n','time',  'OA_SVM_org',  'OA_SVM_rev',  'AA_SVM_org',  'AA_SVM_rev', 'Kappa_SVM_org', 'Kappa_SVM_rev');
% fprintf(file,'  %f  %f  %f  %f  %f  %f  %f \n',[time  OA_SVM1  OA_SVM2  AA_SVM1  AA_SVM2  Kappa_SVM1  Kappa_SVM2]);
% fclose(file);
% %% 
% if n3>3
%     subplot(121);
%     imagesc(T(:,:,100));
%     subplot(122);
%     imagesc(TRec(:,:,100));
% end;
