function best_plot_ever(fig)
    font = 22;
    a = 36; % set this parameter and keep it forever
    b = 0.55; % feel free to play with this ratio
    set(gca,'FontSize',font)
    set(findall(fig,'-property','Box'),'Box','off') % optional
    set(findall(fig,'-property','Interpreter'),'Interpreter','latex')
    set(findall(fig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
    set(fig,'Units','centimeters','Position',[3 3 a b*a])
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
    set(gca, 'XDir', 'normal', 'YDir', 'normal');
    grid on
end