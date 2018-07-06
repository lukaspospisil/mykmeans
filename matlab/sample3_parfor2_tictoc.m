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
Ts = [2000:2000:20000,5e4,1e5,5e5,1e6,5e6,1e7]; % number of points
%Ts = [2000:2000:20000]; % number of points
N = 2; % dimension of the data

time_Gamma0 = zeros(length(Ts),1);
time_all = zeros(length(Ts),1);
time_Theta = zeros(length(Ts),1);
time_Gamma = zeros(length(Ts),1);
time_L = zeros(length(Ts),1);

for Tidx = 1:length(Ts)
    T = Ts(Tidx);
    
    disp(['T = ' num2str(T)])
    
    % get parallel pool or start new one
    p = gcp('nocreate');
    
    % prepare local indexes
    opt_nmb = floor(T/p.NumWorkers);
    rest_nmb = mod(T,p.NumWorkers);
    local_Tidx = cell(1,p.NumWorkers);
    for proc_id=1:p.NumWorkers
        if proc_id == 1
            my_begin = 1;
        else
            my_begin = my_end+1;
        end
        
        my_end = my_begin + opt_nmb - 1;
        if proc_id <= rest_nmb
            my_end = my_end + 1;
        end
        
        local_Tidx{proc_id} = my_begin:my_end;
        
    end
    
    % generate problem
    X = generate_clustering(T,N);
    
    % run K-means algorithm
    % with one annealing step
    
    % generate random feasible gamma
    timer_Gamma0 = tic;
    
    Gamma = zeros(K,T);
    for t=1:T
        randidx = randi([1,K],1,1);
        Gamma(randidx,t) = 1;
    end
    
    time_Gamma0(Tidx) = toc(timer_Gamma0);
    
    it = 0;
    L = Inf;
    Theta = zeros(N,K);
    
    timer_all = tic;
    while it < 1000
        
        % solve Theta-problem
        timer_Theta = tic;
        
        for k=1:K
            sumgamma = sum(Gamma(k,:));
            if sumgamma ~= 0
                for n = 1:N
                    Theta(n,k) = dot(X(n,:),Gamma(k,:))/sumgamma;
                end
            end
        end
        
        time_Theta(Tidx) = time_Theta(Tidx) + toc(timer_Theta);
        
        % solve gamma problem
        
        Theta_d = cell(p.NumWorkers,1);
        Gamma_d = cell(1,p.NumWorkers);
        X_d = cell(1,p.NumWorkers);
        for proc_id=1:p.NumWorkers
            Theta_d{proc_id} = Theta;
            Gamma_d{proc_id} = Gamma(:,local_Tidx{proc_id});
            X_d{proc_id} = X(:,local_Tidx{proc_id});
        end
        
        g_d = cell(p.NumWorkers,1);
        
        timer_Gamma = tic;
        
        parfor proc_id=1:p.NumWorkers
            g_d{proc_id} = zeros(size(Gamma_d{proc_id}));
            for k=1:K
                g_d{proc_id}(k,:) = dot(X_d{proc_id} - Theta_d{proc_id}(:,k),...
                    X_d{proc_id} - Theta_d{proc_id}(:,k));
            end
            [~,idx] = min(g_d{proc_id},[],1);
            Gamma_d{proc_id} = zeros(K,length(local_Tidx{proc_id}));
            for k=1:K
                Gamma_d{proc_id}(k,idx==k) = 1;
            end
        end
        
        time_Gamma(Tidx) = time_Gamma(Tidx) + toc(timer_Gamma);
        
        for proc_id=1:p.NumWorkers
            Gamma(:,local_Tidx{proc_id}) = Gamma_d{proc_id};
        end
        
        % compute function value
        timer_L = tic;
        
        L_old = L;
        
        parfor proc_id=1:p.NumWorkers
            L_d{proc_id} = sum(sum(bsxfun(@times,Gamma_d{proc_id},g_d{proc_id})));
        end
        L = sum([L_d{:}])/(T*N);
        
        time_L(Tidx) = time_L(Tidx) + toc(timer_L);
        
        % check stopping criteria
        Ldelta = L_old - L;
        if Ldelta < 1e-6
            break;
        end
        
        it = it + 1;
    end
    
    time_all(Tidx) = toc(timer_all)/it;
    time_Theta(Tidx) = time_Theta(Tidx)/it;
    time_Gamma(Tidx) = time_Gamma(Tidx)/it;
    time_L(Tidx) = time_L(Tidx)/it;
end

clear Gamma X Theta g
save('results/kmeans3a.mat')

