function [out_kmeans]= KMeans_Classify(X, K, N_anneal)   
    [d,T]=size(X);
    
    [idx_fin,out_kmeans.C,out_kmeans.L]=kmeans(X',K,'Replicates',N_anneal,'MaxIter',1000);
    
    out_kmeans.gamma=zeros(K,T);
    for ttt=1:T
       out_kmeans.gamma(idx_fin(ttt),ttt) = 1;
    end
    out_kmeans.C = out_kmeans.C';

    out_kmeans.L = norm(X - out_kmeans.C*out_kmeans.gamma,'fro')^2/(T*d);

end
