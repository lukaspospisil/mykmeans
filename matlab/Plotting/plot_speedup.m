% prepare data
timess = zeros(length(ns),length(Ts),6);

time_all = 0;
for nidx = 1:length(ns)
    for Tidx = 1:length(Ts)
        for alg = [5,6]
            timess_n_rand = [];
            it_n_rand = [];
            
            % compute averages
            for n_randidx = 1:length(Ns_rand)
                if computed(nidx,Tidx,n_randidx,alg) == 1
%                    timess_n_rand = [timess_n_rand,aa];
                    timess_n_rand = [timess_n_rand,out{nidx,Tidx,n_randidx,alg}.time_its];
%                    timess_n_rand = [timess_n_rand,out{nidx,Tidx,n_randidx,alg}.time_its];
                    it_n_rand = [it_n_rand,out{nidx,Tidx,n_randidx,alg}.it_all];
                end
            end
            
            timess(nidx,Tidx,alg) = mean(timess_n_rand./it_n_rand); %./it_n_rand
        end
    end
end

% plot 3D
if and(length(ns) > 1, length(Ts) > 1)
    [ngrid,Tgrid] = meshgrid(ns,Ts);
    
    figure
    hold on
    
    title('speed up')

    surf(ngrid,Tgrid,timess(:,:,1)'./timess(:,:,2)',...
            'EdgeColor',[1.0,0.0,0.0], ...
            'FaceColor',[1.0,0.5,0.5], ...
            'EdgeAlpha',0.6, ...
            'FaceAlpha',0.8);
    
    
    surf(ngrid,Tgrid,timess(:,:,5)'./timess(:,:,6)',...
            'EdgeColor',[0.0,0.0,1.0], ...
            'FaceColor',[0.2,0.2,1.0], ...
            'EdgeAlpha',0.6, ...
            'FaceAlpha',0.8);
    
    legend('kmeans CPU/kmeans GPU','spam CPU/spam GPU');
    
    xlabel('$n$','Interpreter','latex')
    ylabel('$T$','Interpreter','latex')
    zlabel('speed up','Interpreter','latex')
    
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
%    set(gca, 'ZScale', 'log')
    
    view(30,30)
    grid on
    
    hold off
end

