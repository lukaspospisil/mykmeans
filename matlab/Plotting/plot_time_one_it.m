% prepare data
timess = zeros(length(ns),length(Ts),6);

for nidx = 1:length(ns)
    for Tidx = 1:length(Ts)
        for alg = 1:6
            timess_n_rand = [];
            it_n_rand = [];
            % compute averages
            for n_randidx = 1:length(Ns_rand)
                if computed(nidx,Tidx,n_randidx,alg) == 1
                    timess_n_rand = [timess_n_rand,out{nidx,Tidx,n_randidx,alg}.time];
                    it_n_rand = [it_n_rand,out{nidx,Tidx,n_randidx,alg}.it_all];
                end
            end
            
            timess(nidx,Tidx,alg) = mean(timess_n_rand./it_n_rand);
        end
    end
end

to_plot = [];

to_plot{1}.EdgeColor = [1.0,0.0,0.0];
to_plot{1}.FaceColor = [1.0,0.5,0.5];
to_plot{1}.alg = 1;
to_plot{1}.name = 'Kmeans';

to_plot{2}.EdgeColor = [0.8,0.4,0.4];
to_plot{2}.FaceColor = [1.0,0.6,0.6];
to_plot{2}.alg = 2;
to_plot{2}.name = 'Kmeans GPU';

to_plot{3}.EdgeColor = [0.0,0.0,1.0];
to_plot{3}.FaceColor = [0.2,0.2,1.0];
to_plot{3}.alg = 5;
to_plot{3}.name = 'SPAM';

to_plot{4}.EdgeColor = [0.4,0.4,0.8];
to_plot{4}.FaceColor = [0.6,0.6,1.0];
to_plot{4}.alg = 6;
to_plot{4}.name = 'SPAM GPU';

% plot 3D
if and(length(ns) > 1, length(Ts) > 1)
    [ngrid,Tgrid] = meshgrid(ns,Ts);
    
    figure
    hold on
    
    title('one iteration time')
    
    mylegend = cell(numel(to_plot),1);
    for alg = 1:numel(to_plot)
        surf(ngrid,Tgrid,timess(:,:,to_plot{alg}.alg)',...
            'EdgeColor',to_plot{alg}.EdgeColor, ...
            'FaceColor',to_plot{alg}.FaceColor, ...
            'EdgeAlpha',0.6, ...
            'FaceAlpha',0.8);
        mylegend{alg} = to_plot{alg}.name;
    end
    
    legend(mylegend);
    
    xlabel('$n$','Interpreter','latex')
    ylabel('$T$','Interpreter','latex')
    zlabel('time of one iteration $[s]$','Interpreter','latex')
    
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    set(gca, 'ZScale', 'log')
    
    view(30,30)
    grid on
    
    hold off
end

