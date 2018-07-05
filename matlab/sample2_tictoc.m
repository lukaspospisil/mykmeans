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
Ts = [2000:2000:20000, 5e4, 1e5, 5e5, 1e6]; % number of points
n = 2; % dimension of the data


time_gamma0 = zeros(length(Ts),1);
time_C = zeros(length(Ts),1);
time_gamma = zeros(length(Ts),1);
time_L = zeros(length(Ts),1);

for Tidx = 1:length(Ts)
    T = Ts(Tidx);
    
    disp(['T = ' num2str(T)])
    
    % generate problem
    X = generate_clustering(T,n);
    
    % generate random feasible gamma
    tic
    gamma = rand(K,T);
    gamma_sum = sum(gamma,1);
    for k=1:K
        gamma(k,:) = gamma(k,:)./gamma_sum;
    end
    time_gamma0(Tidx) = toc;
    
    it = 0;
    L = Inf;
    C = zeros(n,K);
    
    while it < 1000
        
        % solve C-problem
        tic
        
        for k=1:K
            sumgamma = sum(gamma(k,:));
            if sumgamma ~= 0
                for nn = 1:n
                    C(nn,k) = dot(X(nn,:),gamma(k,:))/sumgamma;
                end
            end
        end
        
        time_C(Tidx) = time_C(Tidx) + toc;
        
        % solve gamma problem
        tic
        g = zeros(size(gamma));
        for k=1:K
            g(k,:) = dot(X - C(:,k),X - C(:,k));
        end
        [~,idx] = min(g,[],1);
        gamma = zeros(K,T);
        for k=1:K
            gamma(k,idx==k) = 1;
        end

        time_gamma(Tidx) = time_gamma(Tidx) + toc;
        
        % compute function value
        tic
        L_old = L;
        L = sum(sum(bsxfun(@times,gamma,g)))/(T*n);
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

gamma = [];
X = [];
C = [];
save('results/kmeans2.mat')

