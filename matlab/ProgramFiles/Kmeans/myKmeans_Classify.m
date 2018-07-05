function [out_kmeans]= myKmeans_Classify(X, K, N_anneal, myeps)

maxit = 1e4;
[C,gamma,it_all] = mykmeans(X,K,N_anneal,myeps,maxit);

out_kmeans.it_all = it_all;
out_kmeans.gamma=gamma;
out_kmeans.C = C;

out_kmeans.L = get_L( X,C,gamma );

end

function [C_best,gamma_best,it_sum] = mykmeans(X,K,N_anneal,myeps,maxit)
[d,T]=size(X);

C_best = zeros(d,K);
gamma_best = zeros(K,T);
L_best = Inf;
it_sum = 0;

for i=1:N_anneal
    % initial gamma for this annealing step
    gamma=rand(K,T);
    gamma_sum = sum(gamma,1);
    for k=1:K
        gamma(k,:)=gamma(k,:)./gamma_sum;
    end
    
    C = C_best;
    
    L = Inf;
    
    for it=1:maxit
        
        % find new centroids
        for k=1:K
            sumgamma = sum(gamma(k,:));
            if sumgamma ~= 0
                C(:,k) = sum(bsxfun(@times,X,gamma(k,:)),2)/sumgamma;
            end
        end
        
        % find new affiliation vector
        for k=1:K
            gamma(k,:) = dot(X - C(:,k),X - C(:,k)); % reuse this array
        end
        [~,idx] = min(gamma,[],1);
        gamma = zeros(K,T);
        gamma(sub2ind(size(gamma),idx,1:T)) = 1;
        
        % compute objective function
        L_old = L;
        L = norm(X - C*gamma,'fro')^2/(T*d);
        
        if abs(L - L_old) < myeps
            break;
        end
    end
    
    if L < L_best
        L_best = L;
        C_best = C;
        gamma_best = gamma;
    end
    
    it_sum = it_sum + it;
end

end
