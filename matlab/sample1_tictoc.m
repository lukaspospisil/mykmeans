clc
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
Ts = [2000:2000:20000]; % number of points
N = 2; % dimension of the data

time_gamma0 = zeros(length(Ts),1);
time_C = zeros(length(Ts),1);
time_gamma = zeros(length(Ts),1);
time_L = zeros(length(Ts),1);

for Tidx = 1:length(Ts)
    T = Ts(Tidx);
    
    disp(['T = ' num2str(T)])
    
    % generate problem
    X = generate_clustering(T,N);

    % run K-means algorithm
    % with one annealing step

    % generate random feasible gamma
    tic
    
    Gamma = zeros(K,T);
    for t=1:T
        randidx = randi([1,K],1,1);
        Gamma(randidx,t) = 1;
    end
    
    time_gamma0(Tidx) = toc;
    
    it = 0;
    L = Inf;
    Theta = zeros(N,K);
    
    while it < 1000
        
        % solve Theta-problem
        tic
        
        for k=1:K
            sumgamma = sum(Gamma(k,:));
            if sumgamma ~= 0
                for n = 1:N
                    Theta(n,k) = dot(X(n,:),Gamma(k,:))/sumgamma;
                end
            end
        end
        
        time_C(Tidx) = time_C(Tidx) + toc;
        
        % solve gamma problem
        tic
        
        for t=1:T
            g = zeros(K,1);
            for k=1:K
                g(k) = dot(X(:,t) - Theta(:,k),X(:,t) - Theta(:,k));
            end
            
            Gamma(:,t) = zeros(K,1);
            [~,minidx] = min(g);
            Gamma(minidx,t) = 1;
        end
        
        time_gamma(Tidx) = time_gamma(Tidx) + toc;
        
        % compute function value
        tic
        
        L_old = L;
        
        
        L = 0;
        for k=1:K
            for t=1:T
                L = L + ...
                  Gamma(k,t)*dot(X(:,t) - Theta(:,k),X(:,t) - Theta(:,k));
            end
        end
        L = L/(T*N);
        
        time_L(Tidx) = time_L(Tidx) + toc;
        
        % check stopping criteria
        Ldelta = L_old - L;
        if Ldelta < 1e-6
            break;
        end
        
        it = it + 1;
    end
    
    time_C(Tidx) = time_C(Tidx)/it;
    time_gamma(Tidx) = time_gamma(Tidx)/it;
    time_L(Tidx) = time_L(Tidx)/it;
end

%gamma = [];
%X = [];
%C = [];
%save('results/kmeans1.mat')

