
clc
clear all
%close all
addpath('ProgramFiles/')
addpath('ProblemGenerator/')
addpath('Plotting/')

randn('seed',1);
rand('seed',1);

N_anneal=5;

T = 1e4; % the length of time-series
n = 1e1; % dimension of the problem

[X, K_orig, gamma] = generate_clustering_illia(T,n);

K = 3;

% map data to [-1,1]
X=X-(ones(T,1)*mean(X'))';
X=X./max(max(abs(X)));

% compute "standard" K-means
mydisp('- kmeans')
tic
out_kmeans = KMeans_Classify(X,K,N_anneal);
out_kmeans.time = toc;

% compute LSD
if false
    mydisp('- LSD')
    tic
    out_lsd = LSD_Classify(X,K);
    out_lsd.time = toc;
end

% compute SPAM
if true
    mydisp('- SPAM')
    tic
    out_spam = SPAM_Classify(X,K,N_anneal);
    out_spam.time = toc;
end

% disp info
if exist('out_kmeans','var')
    disp(['- KMEANS:  ' num2str(out_kmeans.time) 's,   L = ' num2str(out_kmeans.L)]);
end

if exist('out_lsd','var')
    disp(['- LSD:     ' num2str(out_lsd.time) 's,   L = ' num2str(out_lsd.L)]);
end

if exist('out_spam','var')
    disp(['- SPAM:    ' num2str(out_spam.time) 's,   L = ' num2str(out_spam.L)]);
end

% plot solution
if true
    what_to_plot = [1,2];
    subplots = 3;
    
    figure
    subplot_idx = 1;
    if true
        subplot(1,subplots,subplot_idx)
        title('original')
        plot_classification(X,gamma,what_to_plot);
        subplot_idx = subplot_idx + 1;
    end
    
    if exist('out_kmeans','var') && true
        subplot(1,subplots,subplot_idx)
        title('k-means')
        plot_classification(X,out_kmeans.gamma,what_to_plot);
        subplot_idx = subplot_idx + 1;
    end
    
    if exist('out_lsd','var') && true
        subplot(1,subplots,subplot_idx)
        title('LSD')
        plot_classification_soft(X,out_lsd.gamma,what_to_plot);
        subplot_idx = subplot_idx + 1;
    end
    
    if exist('out_spam','var') && true
        subplot(1,subplots,subplot_idx)
        title('SPAM')
        plot_classification_soft(X,out_spam.gamma,what_to_plot);
        subplot_idx = subplot_idx + 1;
    end
    
end
