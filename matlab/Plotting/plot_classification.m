function plot_classification( X,pi, plot_dimensions)

% plot the solution (clusters = colors)
K = size(pi,1);
mymap = lines(K); % define own pallete of colors


mylegend = cell(K,1);
%figure;
hold on
for k=1:K
    ii=find(pi(k,:)==1);
    if length(plot_dimensions) == 2
        plot(X(plot_dimensions(1),ii),X(plot_dimensions(2),ii),'*','Color',mymap(k,:));
    end
    if length(plot_dimensions) == 3
        plot3(X(plot_dimensions(1),ii),X(plot_dimensions(2),ii),X(plot_dimensions(3),ii),'s','Color',mymap(k,:));
    end
    
    mylegend{k} = ['Label ' num2str(k-1)];
end
legend(mylegend);
xlabel(['$x_{' num2str(plot_dimensions(1)) '}$'],'Interpreter','latex');
ylabel(['$x_{' num2str(plot_dimensions(2)) '}$'],'Interpreter','latex');
if length(plot_dimensions) == 3
    zlabel(['$x_{' num2str(plot_dimensions(3)) '}$'],'Interpreter','latex');
end
set(gca,'FontSize',12,'LineWidth',1);

end

