function plot_classification_soft( X,pi, plot_dimensions)

% this version works with continuous pi, but only for 3 clusters = R,G,B
% plot_dimension defines which components of X will be plotted - should
% work with 2 or 3 dimensions

% plot the solution (clusters = colors)
K = size(pi,1);
T = size(pi,2);

%figure;
hold on
for t=1:T
    mycolor = [0,0,0];
    if K >= 1
       mycolor(1) = pi(1,t); 
    end
    if K >= 2
       mycolor(2) = pi(2,t); 
    end
    if K >= 3
       mycolor(3) = pi(3,t); 
    end
    
    if length(plot_dimensions) == 2
        plot(X(plot_dimensions(1),t),X(plot_dimensions(2),t),'*','Color',mycolor);
    end
    if length(plot_dimensions) == 3
        plot3(X(plot_dimensions(1),t),X(plot_dimensions(2),t),X(plot_dimensions(3),t),'s','Color',mycolor);
    end    
end
xlabel(['$x_{' num2str(plot_dimensions(1)) '}$'],'Interpreter','latex');
ylabel(['$x_{' num2str(plot_dimensions(2)) '}$'],'Interpreter','latex');
if length(plot_dimensions) == 3
    zlabel(['$x_{' num2str(plot_dimensions(3)) '}$'],'Interpreter','latex');
end
set(gca,'FontSize',12,'LineWidth',1);

end

