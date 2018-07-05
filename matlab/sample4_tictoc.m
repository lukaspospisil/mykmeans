clc
%clear all
path(pathdef);

%close all
addpath('ProblemGenerator/')
addpath('Plotting/')

% set random generator
randn('seed',13);
rand('seed',13);

K=3;

% problem parameters
%Ts = [2000:2000:20000]; % number of points
Ts = [2000:2000:20000, 5e4, 1e5, 5e5, 1e6, 5e6, 1e7]; % number of points
n = 2; % dimension of the data


time_gamma0 = zeros(length(Ts),1);
time_C = zeros(length(Ts),1);
time_gamma = zeros(length(Ts),1);
time_L = zeros(length(Ts),1);

% initialize parallel enviroment
numOfWorkers = feature('numcores'); % number of parallel processes
if ~exist('poolobj')
	poolobj = parpool(numOfWorkers); % start parallel environment
end

for Tidx = 1:length(Ts)
    T = Ts(Tidx);
    
    disp(['T = ' num2str(T)])
    
    % generate problem
    X = generate_clustering(T,n);
    X = distributed(X); % distribute data
    
    % generate random feasible gamma
    tic

    spmd
        X_local = getLocalPart(X);
        T_local = size(X_local,2);
        gamma_local = rand(K,T_local);
        gamma_sum = sum(gamma_local,1);
        for k=1:K
            gamma_local(k,:) = gamma_local(k,:)./gamma_sum;
        end
        
        L = Inf;
    end

    time_gamma0(Tidx) = toc;
    
    it = 0;
    
    while it < 1000
        
        % solve C-problem
        tic

        spmd
            sumgamma_local = zeros(K,1);
            mydot_local = zeros(n,K);
            for k=1:K
                sumgamma_local(k) = sum(gamma_local(k,:));
                for nn = 1:n
                    mydot_local(nn,k) = dot(X_local(nn,:),gamma_local(k,:));
                end
            end
            
            % sum local contributions
            sumgamma = gop(@plus,sumgamma_local);
            mydot = gop(@plus,mydot_local);
            
            % compute new C
            C = zeros(n,K);
            for k=1:K
                if sumgamma(k) ~= 0
                    for nn = 1:n
                        C(nn,k) = mydot(nn,k)/sumgamma(k);
                    end
                end
            end
        end
        
        time_C(Tidx) = time_C(Tidx) + toc;
        
        % solve gamma problem
        tic

        spmd
            g_local = zeros(size(gamma_local));
            for k=1:K
                g_local(k,:) = dot(X_local - C(:,k),X_local - C(:,k));
            end
            [~,idx] = min(g_local,[],1);
            gamma_local = zeros(K,T_local);
            for k=1:K
                gamma_local(k,idx==k) = 1;
            end
        end
            
        time_gamma(Tidx) = time_gamma(Tidx) + toc;
        
        % compute function value
        tic
        L_old = L;
        spmd
            L_local = sum(sum(bsxfun(@times,gamma_local,g_local)));
            L = gop(@plus,L_local)/(T*n);
        end
        
        time_L(Tidx) = time_L(Tidx) + toc;
        
        % check stopping criteria
        Ldelta = L_old{1} - L{1};
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
save('results/kmeans3.mat')

