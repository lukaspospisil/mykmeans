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
%Ts = [2000:2000:20000]; % number of points
Ts = [2000:2000:20000, 5e4, 1e5, 5e5, 1e6]; % number of points
n = 2; % dimension of the data


time_gamma0 = zeros(length(Ts),1);
time_C = zeros(length(Ts),1);
time_gamma = zeros(length(Ts),1);
time_L = zeros(length(Ts),1);

% initialize parallel enviroment
numOfWorkers = feature('numcores'); % number of parallel processes
if ~exist('poolobj')
    %	poolobj = parpool(numOfWorkers); % start parallel environment
end

for Tidx = 1:length(Ts)
    T = Ts(Tidx);
    
    disp(['T = ' num2str(T)])
    
    % generate problem
    X = generate_clustering(T,n);
    
    % prepare local indexes
    opt_nmb = floor(T/numOfWorkers);
    rest_nmb = mod(T,numOfWorkers);
    
    local_Tidx = cell(numOfWorkers+1,1);
    X_local = cell(numOfWorkers,1);
    T_local = cell(numOfWorkers,1);
    g_local = cell(numOfWorkers,1);
    C_local = cell(numOfWorkers,1);    

    local_Tidx{1} = 1;
    for proc_id=1:numOfWorkers
        local_Tidx{proc_id+1} = local_Tidx{proc_id} + opt_nmb;
        if proc_id <= rest_nmb
            local_Tidx{proc_id+1} = local_Tidx{proc_id+1} + 1;
        end
        
        X_local{proc_id} = X(:,local_Tidx{proc_id} : (local_Tidx{proc_id+1}-1));
        T_local{proc_id} = local_Tidx{proc_id+1} - local_Tidx{proc_id};
        L_local{proc_id} = 1;
    end
    
    % generate random feasible gamma
    tic
    
    for proc_id = 1:numOfWorkers
        gamma_local{proc_id} = rand(K,T_local{proc_id});
        gamma_sum{proc_id} = sum(gamma_local{proc_id},1);
        for k=1:K
            gamma_local{proc_id}(k,:) = gamma_local{proc_id}(k,:)./gamma_sum{proc_id};
        end
    end

    L = Inf;
    
    time_gamma0(Tidx) = toc;
    
    it = 0;
    
    while it < 1000
        
        % solve C-problem
        tic
        
        parfor proc_id = 1:numOfWorkers
            sumgamma_local{proc_id} = zeros(K,1);
            mydot_local{proc_id} = zeros(n,K);
            for k=1:K
                sumgamma_local{proc_id}(k) = sum(gamma_local{proc_id}(k,:));
                for nn = 1:n
                    mydot_local{proc_id}(nn,k) = dot(X_local{proc_id}(nn,:),gamma_local{proc_id}(k,:));
                end
            end
        end
        
        % sum local contributions
        sumgamma = plus(sumgamma_local{:});
        mydot = plus(mydot_local{:});
        
        % compute new C
        C = zeros(n,K);
        for k=1:K
            if sumgamma(k) ~= 0
                for nn = 1:n
                    C(nn,k) = mydot(nn,k)/sumgamma(k);
                end
            end
        end
        time_C(Tidx) = time_C(Tidx) + toc;
        
        % distribute C
        for proc_id = 1:numOfWorkers
            C_local{proc_id} = C;
        end        
        
        % solve gamma problem
        tic
        
        for proc_id = 1:numOfWorkers
            g_local{proc_id} = zeros(size(gamma_local{proc_id}));
            for k=1:K
                g_local{proc_id}(k,:) = dot(X_local{proc_id} - C_local{proc_id}(:,k),X_local{proc_id} - C_local{proc_id}(:,k));
            end
            [~,idx] = min(g_local{proc_id},[],1);
            gamma_local{proc_id} = zeros(K,T_local{proc_id});
            for k=1:K
                gamma_local{proc_id}(k,idx==k) = 1;
            end
        end
        
        time_gamma(Tidx) = time_gamma(Tidx) + toc;
        
        % compute function value
        tic
        L_old = L;
        parfor proc_id = 1:numOfWorkers
            L_local{proc_id} = sum(sum(bsxfun(@times,gamma_local{proc_id},g_local{proc_id})));
        end
        L = plus(L_local{:});
 
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
save('results/kmeans3.mat')

