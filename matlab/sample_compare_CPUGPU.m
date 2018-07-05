clc
clear all
path(pathdef);

%close all
addpath('ProgramFiles/')
addpath('ProgramFiles/Kmeans')
addpath('ProgramFiles/KmeansGPU')
addpath('ProblemGenerator/')
addpath('Plotting/')

% prepare GPU
g = gpuDevice(1);
reset(g);
wait(g);

K=3;

% set number of annealing steps
N_anneal=10;

% the length of time-series
Ts = 2.^[10:18];
%Ts = 2.^[15:17];

% the dimension of the problem
ns = [4,16];%[4,16,64,256,1024];

% number of random generation of the problem
N_rand = 1;

% stopping criteria for kmeans
myeps = 1e-6;

compute_kmeansMatlab = true;
compute_kmeans = true;
compute_kmeansgpu = true;

timer_all = tic;
for n_rand=1:N_rand
    
    for nidx = 1:length(ns)
        for Tidx = 1:length(Ts)
            
            
            disp([' - n = ' num2str(ns(nidx)) ' (' num2str(nidx) '/' num2str(length(ns)) '), ' ...
                'T = ' num2str(Ts(Tidx)) ' (' num2str(Tidx) '/' num2str(length(Ts)) '), ' ...
                'n_rand = ' num2str(n_rand) '/' num2str(N_rand)])
            
            %%% Create artificial point-clustering "k-means" problem
            T = Ts(Tidx); % the length of time-series
            n = ns(nidx); % dimension of the problem
            
            % prepare random seed
            randn('seed',n_rand);
            rand('seed',n_rand);
            
            X = generate_clustering_illia(T,n);
            
            % map data to [-1,1]
            X=X-(ones(T,1)*mean(X'))';
            X=X./max(max(abs(X)));

            % compute "standard" K-means
            if compute_kmeansMatlab
                disp('    - kmeansMatlab')
                
                mytimer_kmeans = tic;
                out_kmeans= Kmeans_Classify(X,K,N_anneal,myeps);
                out_kmeans.time = toc(mytimer_kmeans);
                out_kmeans.gamma = []; % throw away solution
                out_kmeans.C = []; % throw away solution
                out_kmeans.computed = true;
            else
                out_kmeans.L = Inf;
                out_kmeans.time = -1;
                out_kmeans.time_its = -1;
                out_kmeans.it_all = 0;
                out_kmeans.computed = false;
            end
            
            % compute "standard" K-means
            if compute_kmeans
                disp('    - kmeans')
                
                mytimer_kmeans = tic;
                out_kmeans= myKmeans_Classify(X,K,N_anneal,myeps);
                out_kmeans.time = toc(mytimer_kmeans);
                out_kmeans.gamma = []; % throw away solution
                out_kmeans.C = []; % throw away solution
                out_kmeans.computed = true;
            else
                out_kmeans.L = Inf;
                out_kmeans.time = -1;
                out_kmeans.time_its = -1;
                out_kmeans.it_all = 0;
                out_kmeans.computed = false;
            end
            
            % compute "standard" K-means
            if compute_kmeansgpu
                disp('    - kmeansGPU')
                mytimer_kmeansgpu = tic;
                out_kmeansgpu = myKmeansGPU_Classify(X,K,N_anneal,myeps);
                out_kmeansgpu.time = toc(mytimer_kmeansgpu);
                out_kmeansgpu.gamma = []; % throw away solution
                out_kmeansgpu.C = []; % throw away solution
                out_kmeansgpu.computed = true;
            else
                out_kmeansgpu.L = Inf;
                out_kmeansgpu.time = -1;
                out_kmeansgpu.time_its = -1;
                out_kmeansgpu.it_all = 0;
                out_kmeansgpu.computed = false;
            end
            
            X = [];
        end
    end
end

save("results/compare_CPUGPU.mat")

disp('-------------------------------------------')
disp(['total time : ' num2str(toc(timer_all)) ' s'])
disp('-------------------------------------------')
