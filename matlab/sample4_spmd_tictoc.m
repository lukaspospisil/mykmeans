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
    
    % generate problem
    
    X_orig = generate_clustering(T,N);
    X_dist = distributed(X_orig);
    
    spmd
        X_local = getLocalPart(X_dist);
        T_local = size(X_local,2);
        
        % generate random feasible gamma
        timer_Gamma0 = tic;
        
        Gamma_local = zeros(K,T_local);
        for t=1:T_local
            randidx = randi([1,K],1,1);
            Gamma_local(randidx,t) = 1;
        end
        time_Gamma0(Tidx) = toc(timer_Gamma0);
        
        it = 0;
        L = Inf;
        Theta = zeros(N,K);

        timer_all = tic;        
        stop_it = false;
        while it < 1000 && ~stop_it
        
            % solve Theta-problem
            timer_Theta = tic;
            for k=1:K
                sumgamma_local = sum(Gamma_local(k,:));
                sumgamma = gop(@plus,sumgamma_local);
                
                for n = 1:N
                    mydot_local = dot(X_local(n,:),Gamma_local(k,:));
                    mydot = gop(@plus,mydot_local);
                
                    if sumgamma ~= 0
                        Theta(n,k) = mydot/sumgamma;
                    end
                end
            end
            time_Theta(Tidx) = time_Theta(Tidx) + toc(timer_Theta);            
        
            % solve Gamma problem
            timer_Gamma = tic;
            g_local = zeros(size(Gamma_local));
            for k=1:K
                g_local(k,:) = dot(X_local - Theta(:,k),X_local - Theta(:,k));
            end
            [~,idx] = min(g_local,[],1);

            Gamma_local = zeros(K,T_local);
            for k=1:K
                Gamma_local(k,idx==k) = 1;
            end
            time_Gamma(Tidx) = time_Gamma(Tidx) + toc(timer_Gamma);        
            
            
            % compute function value
            L_old = L;
            
            timer_L = tic;
            L_local = sum(sum(bsxfun(@times,Gamma_local,g_local)));
            L = gop(@plus, L_local)/(T*N);
            time_L(Tidx) = time_L(Tidx) + toc(timer_L);
            
            % check stopping criteria
            Ldelta = L_old - L;
            if Ldelta < 1e-6
                stop_it = true;
            end
        
            it = it + 1;
        end
        
        time_all(Tidx) = toc(timer_all)/it;
        time_Theta(Tidx) = time_Theta(Tidx)/it;
        time_Gamma(Tidx) = time_Gamma(Tidx)/it;
        time_L(Tidx) = time_L(Tidx)/it;
        
    end

end

p = gcp('nocreate');
time_all = sum([time_all{:}],2)/p.NumWorkers;
time_Theta = sum([time_Theta{:}],2)/p.NumWorkers;
time_Gamma = sum([time_Gamma{:}],2)/p.NumWorkers;
time_L = sum([time_L{:}],2)/p.NumWorkers;

clear Gamma X Theta g timer_L timer_all timer_Gamma timer_Gamma0 g_local timer_Theta X_dist X_local X_orig time_Gamma0 p randidx Gamma_local idx it k L L_local L_old Ldelta mydot mydot_local n stop_it sumgamma sumgamma_local t T_local
save('results/kmeans4a.mat')

