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
T = 1000; % number of points
n = 2; % dimension of the data

% generate problem
X = generate_clustering(T,n);


% run K-means algorithm

% generate random feasible gamma
gamma = rand(K,T);
for t=1:T
    gamma(:,t) = gamma(:,t)/sum(gamma(:,t));
end

it = 0;
L = Inf;
C = zeros(n,K);

while it < 1000
    
    % solve C-problem
    for k=1:K
        sumgamma = sum(gamma(k,:));
        if sumgamma ~= 0
            for nn = 1:n
                C(nn,k) = dot(X(nn,:),gamma(k,:))/sumgamma;
            end
        end
    end
    
    % solve gamma problem
    for t=1:T
        g = zeros(K,1);
        for k=1:K
            g(k) = dot(X(:,t) - C(:,k),X(:,t) - C(:,k));
        end
        
        gamma(:,t) = zeros(K,1);
        [~,minidx] = min(g);
        gamma(minidx,t) = 1;
    end
    
    % compute function value
    L_old = L;
    L = 0;
    for k=1:K
        for t=1:T
            L = L + gamma(k,t)*dot(X(:,t) - C(:,k),X(:,t) - C(:,k));
        end
    end
    L = L/(T*n);
    
    % check stopping criteria
    Ldelta = L_old - L;
    if Ldelta < 1e-6
        break;
    end
    
    it = it + 1;
end


plot_classification( X,gamma, 1:2);
