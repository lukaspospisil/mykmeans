function [out_kmeans]= myKmeansGPU_Classify(X, K, N_anneal, myeps)

maxit = 1e4;
[C,gamma,it_all,time_its] = mykmeansGPU(X,K,N_anneal,myeps,maxit);

out_kmeans.it_all = it_all;
out_kmeans.gamma=gamma;
out_kmeans.C = C;
out_kmeans.time_its = time_its;

out_kmeans.L = get_L( X,C,gamma );

end

function [C_best,gamma_best,it_sum,time_its] = mykmeansGPU(X,K,N_anneal,myeps,maxit)
[d,T]=size(X);

Xgpu = gpuArray(X);
Cgpu_best = gpuArray(zeros(d,K,'double'));
gammagpu_best = gpuArray(zeros(K,T,'double'));

L_best = Inf;
it_sum = 0;
time_its = 0;

for i=1:N_anneal
    % initial gamma for this annealing step
    gammagpu=gpuArray(rand(K,T));
    gammagpu_sum = sum(gammagpu,1);
    for k=1:K
        gammagpu(k,:)=gammagpu(k,:)./gammagpu_sum;
    end
    
    Cgpu = Cgpu_best;
    
    L = Inf;

    timer_its = tic;
    for it=1:maxit
        
        % find new centroids
        for k=1:K
            sumgammagpu = sum(gammagpu(k,:));
            if sumgammagpu ~= 0
                Cgpu(:,k) = sum(bsxfun(@times,Xgpu,gammagpu(k,:)),2)/sumgammagpu;
            end
        end
        
        % find new affiliation vector
        for k=1:K
            gammagpu(k,:) = dot(Xgpu - Cgpu(:,k),Xgpu - Cgpu(:,k));
        end
        [~,idx] = min(gammagpu,[],1);
        gammagpu = gpuArray(zeros(K,T,'double'));
        gammagpu(sub2ind(size(gammagpu),idx,1:T)) = 1;
        
        % compute objective function
        L_old = L;
        L = norm(Xgpu - Cgpu*gammagpu,'fro')^2/(T*d);
        
        if abs(L - L_old) < myeps
            break;
        end
    end
    time_its = time_its + toc(timer_its);
    
    if L < L_best
        L_best = L;
        Cgpu_best = Cgpu;
        gammagpu_best = gammagpu;
    end
    
    it_sum = it_sum + it;
end

C_best = gather(Cgpu_best);
gamma_best = gather(gammagpu_best);

end


