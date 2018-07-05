clear all

kmeans1 = load('results/kmeans1.mat');
kmeans2 = load('results/kmeans2.mat');
kmeans3 = load('results/kmeans3.mat');

mylegend{1} = '1st implementation';
mylegend{2} = 'vectorization';
mylegend{3} = 'spmd';

plot1 = true;
plot2 = true;
plot3 = true;

setlog = true;

figure
hold on
title('$\Gamma$ computation','interpreter','latex')
if plot1
    plot(kmeans1.Ts,kmeans1.time_gamma,'bs-')
end
if plot2
    plot(kmeans2.Ts,kmeans2.time_gamma,'rs-')
end
if plot3
    plot(kmeans3.Ts,kmeans3.time_gamma,'gs-','Color',[0,0.6,0])
end
legend(mylegend)
xlabel('$T$','Interpreter','latex')
ylabel('time $[s]$','Interpreter','latex')
if setlog
    set(gca,'Xscale','log')
    set(gca,'Yscale','log')
end
hold off

figure
hold on
title('$C$ computation','interpreter','latex')
if plot1
    plot(kmeans1.Ts,kmeans1.time_C,'bs-')
end
if plot2
    plot(kmeans2.Ts,kmeans2.time_C,'rs-')
end
if plot3
    plot(kmeans3.Ts,kmeans3.time_C,'gs-','Color',[0,0.6,0])
end
legend(mylegend)
xlabel('$T$','Interpreter','latex')
ylabel('time $[s]$','Interpreter','latex')
if setlog
    set(gca,'Xscale','log')
    set(gca,'Yscale','log')
end
hold off

figure
hold on
title('$L$ computation','interpreter','latex')
if plot1
    plot(kmeans1.Ts,kmeans1.time_L,'bs-')
end
if plot2
    plot(kmeans2.Ts,kmeans2.time_L,'rs-')
end
if plot3
    plot(kmeans3.Ts,kmeans3.time_L,'gs-','Color',[0,0.6,0])
end
legend(mylegend)
xlabel('$T$','Interpreter','latex')
ylabel('time $[s]$','Interpreter','latex')
if setlog
    set(gca,'Xscale','log')
    set(gca,'Yscale','log')
end
hold off

