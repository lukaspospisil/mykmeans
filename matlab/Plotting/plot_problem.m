function plot_problem( X, plot_dimensions)

% plot the solution (clusters = colors)

%figure;
hold on
if length(plot_dimensions) == 2
    plot(X(plot_dimensions(1),:),X(plot_dimensions(2),:),'*','Color',[0.4,0.4,0.4]);
end
if length(plot_dimensions) == 3
    plot3(X(plot_dimensions(1),:),X(plot_dimensions(2),:),X(plot_dimensions(3),ii),'s','Color',[0.4,0.4,0.4]);
end
xlabel(['$x_{' num2str(plot_dimensions(1)) '}$'],'Interpreter','latex');
ylabel(['$x_{' num2str(plot_dimensions(2)) '}$'],'Interpreter','latex');
if length(plot_dimensions) == 3
    zlabel(['$x_{' num2str(plot_dimensions(3)) '}$'],'Interpreter','latex');
end
set(gca,'FontSize',12,'LineWidth',1);

end

