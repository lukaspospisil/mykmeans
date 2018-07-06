clear all

plot1 = true;
plot2 = true;
plot3a = true;
plot3b = true;

mylegend{1} = '1st implementation';
mylegend{2} = 'vectorization';
mylegend{3} = 'parfor, 2 workers';
mylegend{4} = 'parfor, 4 workers';

if plot1
    kmeans1 = load('results/kmeans1.mat');
end

if plot2
    kmeans2 = load('results/kmeans2.mat');
end

if plot3a
    kmeans3a = load('results/kmeans3a.mat');
end

if plot3b
    kmeans3b = load('results/kmeans3b.mat');
end


setlog = true;

figure
hold on
title('$\Gamma$ computation','interpreter','latex')
if plot1
    plot(kmeans1.Ts,kmeans1.time_Gamma,'bs-')
end
if plot2
    plot(kmeans2.Ts,kmeans2.time_Gamma,'rs-')
end
if plot3a
    plot(kmeans3a.Ts,kmeans3a.time_Gamma,'gs-','Color',[0,0.6,0])
end
if plot3b
    plot(kmeans3b.Ts,kmeans3b.time_Gamma,'ms-','Color',[0.5,0.3,1.0])
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
title('$\Theta$ computation','interpreter','latex')
if plot1
    plot(kmeans1.Ts,kmeans1.time_Theta,'bs-')
end
if plot2
    plot(kmeans2.Ts,kmeans2.time_Theta,'rs-')
end
if plot3a
    plot(kmeans3a.Ts,kmeans3a.time_Theta,'gs-','Color',[0,0.6,0])
end
if plot3b
    plot(kmeans3b.Ts,kmeans3b.time_Theta,'ms-','Color',[0.5,0.3,1.0])
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
if plot3a
    plot(kmeans3a.Ts,kmeans3a.time_L,'gs-','Color',[0,0.6,0])
end
if plot3b
    plot(kmeans3b.Ts,kmeans3b.time_L,'ms-','Color',[0.5,0.3,1.0])
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
title('one iteration','interpreter','latex')
if plot1
    plot(kmeans1.Ts,kmeans1.time_all,'bs-')
end
if plot2
    plot(kmeans2.Ts,kmeans2.time_all,'rs-')
end
if plot3a
    plot(kmeans3a.Ts,kmeans3a.time_all,'gs-','Color',[0,0.6,0])
end
if plot3b
    plot(kmeans3b.Ts,kmeans3b.time_all,'ms-','Color',[0.5,0.3,1.0])
end

legend(mylegend)
xlabel('$T$','Interpreter','latex')
ylabel('time $[s]$','Interpreter','latex')
if setlog
    set(gca,'Xscale','log')
    set(gca,'Yscale','log')
end
hold off



