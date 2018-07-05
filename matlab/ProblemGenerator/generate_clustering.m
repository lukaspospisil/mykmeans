function [ X,K, gamma, mu ] = generate_clustering( T,n )
%GENERATE_CLUSTERING Generate artificial point-clustering "k-means" problem
%   X - generated data, each point is stored in column

K = 4;

% fixed problem parameters
sigma{1}=0.2/(n-2)*eye(n);
sigma{1}(1:2,1:2)=[0.1 0.05;0.05 0.1];

sigma{2}=0.2/(n-2)*eye(n);
sigma{2}(1:2,1:2)=[0.1 -0.05;-0.05 0.1];

sigma{3}=0.2/(n-2)*eye(n);
sigma{3}(1:2,1:2)=eye(2);

sigma{4}=0.2/(n-2)*eye(n);
sigma{4}(1:2,1:2)=eye(2);

mu1=[0;0;zeros(n-2,1)];
mu2=[0.8;1.6;zeros(n-2,1)];
mu3=[1.6;0;zeros(n-2,1)];
mu4=[0.8;0.8;zeros(n-2,1)];

mu = [mu1,mu2,mu3,mu4];

gamma = zeros(K,T);
X = zeros(n,T);
for k=1:K
    Tsize = T/K;
    t_begin = floor((k-1)*Tsize)+1;
    t_end = min(floor(k*Tsize),T);
    
    X(:,t_begin:t_end) = mvnrnd(mu(:,k),sigma{k},t_end-t_begin+1)'; 
    
    gamma(k,t_begin:t_end) = ones(1,t_end-t_begin+1);
end


end

