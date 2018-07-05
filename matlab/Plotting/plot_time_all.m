% prepare data
timess = zeros(length(ns),length(Ts),6);
timess_S.all = zeros(length(ns),length(Ts),6);
timess_S.assemble = zeros(length(ns),length(Ts),6);
timess_S.solve = zeros(length(ns),length(Ts),6);
timess_Gamma.all = zeros(length(ns),length(Ts),6);
timess_Gamma.assemble = zeros(length(ns),length(Ts),6);
timess_Gamma.solve = zeros(length(ns),length(Ts),6);
timess_L = zeros(length(ns),length(Ts),6);
times_test = zeros(length(ns),length(Ts),6);

time_all = 0;
for nidx = 1:length(ns)
    for Tidx = 1:length(Ts)
        for alg = [5,6]
            timess_n_rand = [];
            timess_n_rand_S.all = [];
            timess_n_rand_S.assemble = [];
            timess_n_rand_S.solve = [];
            timess_n_rand_Gamma.all = [];
            timess_n_rand_Gamma.assemble = [];
            timess_n_rand_Gamma.solve = [];
            timess_n_rand_L = [];
            timess_n_rand_test = [];

            it_n_rand = [];
            
            % compute averages
            for n_randidx = 1:length(Ns_rand)
                if computed(nidx,Tidx,n_randidx,alg) == 1

%                    timess_n_rand = [timess_n_rand,aa];
                    
                    timess_n_rand = [timess_n_rand,out{nidx,Tidx,n_randidx,alg}.time_its];

                    timess_n_rand_S.all = [timess_n_rand_S.all,out{nidx,Tidx,n_randidx,alg}.time_sub.time_S.all];
                    timess_n_rand_S.assemble = [timess_n_rand_S.assemble,out{nidx,Tidx,n_randidx,alg}.time_sub.time_S.assemble];
                    timess_n_rand_S.solve = [timess_n_rand_S.solve,out{nidx,Tidx,n_randidx,alg}.time_sub.time_S.solve];

                    timess_n_rand_Gamma.all = [timess_n_rand_Gamma.all,out{nidx,Tidx,n_randidx,alg}.time_sub.time_Gamma.all];
                    timess_n_rand_Gamma.assemble = [timess_n_rand_Gamma.assemble,out{nidx,Tidx,n_randidx,alg}.time_sub.time_Gamma.assemble];
                    timess_n_rand_Gamma.solve = [timess_n_rand_Gamma.solve,out{nidx,Tidx,n_randidx,alg}.time_sub.time_Gamma.solve];

                    timess_n_rand_test = [timess_n_rand_test,gather(out{nidx,Tidx,n_randidx,alg}.it_spgqp_all)];
                    
                    timess_n_rand_L = [timess_n_rand_L,out{nidx,Tidx,n_randidx,alg}.time_sub.time_L];
%                    timess_n_rand = [timess_n_rand,out{nidx,Tidx,n_randidx,alg}.time_its];
                    it_n_rand = [it_n_rand,out{nidx,Tidx,n_randidx,alg}.it_all];
                end
            end
            
            timess(nidx,Tidx,alg) = mean(timess_n_rand./it_n_rand); %./it_n_rand

            timess_S.all(nidx,Tidx,alg) = mean((timess_n_rand_S.all)./it_n_rand);
            timess_S.assemble(nidx,Tidx,alg) = mean((timess_n_rand_S.assemble)./it_n_rand);
            timess_S.solve(nidx,Tidx,alg) = mean((timess_n_rand_S.solve)./it_n_rand);

            timess_Gamma.all(nidx,Tidx,alg) = mean((timess_n_rand_Gamma.all)./it_n_rand);
            timess_Gamma.assemble(nidx,Tidx,alg) = mean((timess_n_rand_Gamma.assemble)./it_n_rand);
            timess_Gamma.solve(nidx,Tidx,alg) = mean((timess_n_rand_Gamma.solve)./it_n_rand);
            
            timess_L(nidx,Tidx,alg) = mean(timess_n_rand_L./it_n_rand); %./it_n_rand
             
        end
    end
end

% plot 3D
if and(length(ns) > 1, length(Ts) > 1)
    [ngrid,Tgrid] = meshgrid(ns,Ts);
    
    edgecolor1 = [0.0,0.0,1.0];
    facecolor1 = [0.2,0.2,1.0];

    edgecolor2 = [0.4,0.4,0.8];
    facecolor2 = [0.6,0.6,1.0];

    figure
    
    subplot(3,3,1)
    hold on
    title('total time')
    plot_all_one(timess(:,:,5)',ngrid,Tgrid,edgecolor1, facecolor1);
    plot_all_one(timess(:,:,6)',ngrid,Tgrid,edgecolor2, facecolor2);
    legend('CPU','GPU')
    zlabel('time [s]','Interpreter','latex')
    hold off
    
    subplot(3,3,2)
    hold on
    title('L computation')
    plot_all_one(timess_L(:,:,5)',ngrid,Tgrid,edgecolor1, facecolor1);
    plot_all_one(timess_L(:,:,6)',ngrid,Tgrid,edgecolor2, facecolor2);
    zlabel('time [s]','Interpreter','latex')
    hold off

    
    subplot(3,3,4)
    hold on
    title('Gamma all')
    plot_all_one(timess_Gamma.all(:,:,5)',ngrid,Tgrid,edgecolor1, facecolor1);
    plot_all_one(timess_Gamma.all(:,:,6)',ngrid,Tgrid,edgecolor2, facecolor2);
    zlabel('time [s]','Interpreter','latex')
    hold off

    subplot(3,3,5)
    hold on
    title('Gamma assemble')
    plot_all_one(timess_Gamma.assemble(:,:,5)',ngrid,Tgrid,edgecolor1, facecolor1);
    plot_all_one(timess_Gamma.assemble(:,:,6)',ngrid,Tgrid,edgecolor2, facecolor2);
    zlabel('time [s]','Interpreter','latex')
    hold off

    subplot(3,3,6)
    hold on
    title('Gamma solve')
    plot_all_one(timess_Gamma.solve(:,:,5)',ngrid,Tgrid,edgecolor1, facecolor1);
    plot_all_one(timess_Gamma.solve(:,:,6)',ngrid,Tgrid,edgecolor2, facecolor2);
    zlabel('time [s]','Interpreter','latex')
    hold off
    
    
    subplot(3,3,7)
    hold on
    title('S all')
    plot_all_one(timess_S.all(:,:,5)',ngrid,Tgrid,edgecolor1, facecolor1);
    plot_all_one(timess_S.all(:,:,6)',ngrid,Tgrid,edgecolor2, facecolor2);
    zlabel('time [s]','Interpreter','latex')
    hold off

    subplot(3,3,8)
    hold on
    title('S assemble')
    plot_all_one(timess_S.assemble(:,:,5)',ngrid,Tgrid,edgecolor1, facecolor1);
    plot_all_one(timess_S.assemble(:,:,6)',ngrid,Tgrid,edgecolor2, facecolor2);
    zlabel('time [s]','Interpreter','latex')
    hold off

    subplot(3,3,9)
    hold on
    title('S solve')
    plot_all_one(timess_S.solve(:,:,5)',ngrid,Tgrid,edgecolor1, facecolor1);
    plot_all_one(timess_S.solve(:,:,6)',ngrid,Tgrid,edgecolor2, facecolor2);
    zlabel('time [s]','Interpreter','latex')
    hold off

    
    
end

