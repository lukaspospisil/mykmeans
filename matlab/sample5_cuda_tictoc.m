clear all
path(pathdef);

%close all
addpath('ProblemGenerator/')
addpath('Plotting/')

% set random generator
randn('seed',13);
rand('seed',13);

K=3;

% problem parameters
Ts = [2000:2000:20000,5e4,1e5,5e5,1e6, 5e6,1e7]; % number of points
%Ts = [2000:2000:20000]; % number of points
N = 2; % dimension of the data

% prepare device
gpu = gpuDevice(1);
reset(gpu);

time_Gamma0 = zeros(length(Ts),1);
time_all = zeros(length(Ts),1);
time_Theta = zeros(length(Ts),1);
time_Gamma = zeros(length(Ts),1);
time_L = zeros(length(Ts),1);

for Tidx = 1:length(Ts)
    T = Ts(Tidx);
    
    disp(['T = ' num2str(T)])
    
    % generate problem
    Xcpu = generate_clustering(T,N);
    X = gpuArray(Xcpu);

    % run K-means algorithm
    % with one annealing step

    % generate random feasible gamma
    wait(gpu);
    timer_Gamma0 = tic;
    
    Gamma = gpuArray.zeros(K,T);
    randidx = gpuArray.randi([1,K],1,T);
    for k=1:K
        Gamma(k,randidx==k) = 1;
    end
    
    wait(gpu);
    time_Gamma0(Tidx) = toc(timer_Gamma0);
    
    it = 0;
    L = Inf;
    Theta = gpuArray.zeros(N,K);

    timer_all = tic;
    while it < 1000
        
        % solve Theta-problem
        wait(gpu);
        timer_Theta = tic;
        
        for k=1:K
            sumgamma = sum(Gamma(k,:));
            if sumgamma ~= 0
                for n = 1:N
                    Theta(n,k) = dot(X(n,:),Gamma(k,:))/sumgamma;
                end
            end
        end
        
        wait(gpu);
        time_Theta(Tidx) = time_Theta(Tidx) + toc(timer_Theta);
        
        % solve gamma problem
        wait(gpu);
        timer_Gamma = tic;
        
        g = gpuArray.zeros(size(Gamma));
        for k=1:K
            g(k,:) = dot(X - Theta(:,k),X - Theta(:,k));
        end
        [~,idx] = min(g,[],1);
        Gamma = gpuArray.zeros(K,T);
        for k=1:K
            Gamma(k,idx==k) = 1;
        end
        
        wait(gpu);
        time_Gamma(Tidx) = time_Gamma(Tidx) + toc(timer_Gamma);
        
        % compute function value
        L_old = L;

        wait(gpu);
        timer_L = tic;

        L = sum(sum(bsxfun(@times,Gamma,g)))/(T*N);
        
        wait(gpu);
        time_L(Tidx) = time_L(Tidx) + toc(timer_L);
        
        % check stopping criteria
        Ldelta = L_old - L;
        if Ldelta < 1e-6
            break;
        end
        
        it = it + 1;
    end

    wait(gpu);
    time_all(Tidx) = toc(timer_all)/it;
    time_Theta(Tidx) = time_Theta(Tidx)/it;
    time_Gamma(Tidx) = time_Gamma(Tidx)/it;
    time_L(Tidx) = time_L(Tidx)/it;
end

clear Gamma X Theta g
save('results/kmeans5.mat')

