
%%%%%%%%%%%%%%%% cas (iv)/cas(v) %%%%%%%%%%%%%%%%%%%%%%
epsilon = [2,0.5,0.1,0.05];
epsinv = 1./epsilon;

%relerrH1 = [  0.3195  ,  0.1986  ,  0.2043  ,  0.1661]; % cas (iv) avec f = 1 
relerrH1 = [1.0546  ,  0.2658  ,  0.2389  ,  0.1705]; % cas (v) avec f = 2*pi^2*sin(pi*x)*sin(pi*y)
relerrL2 = [ 1.0579  ,  0.1102  ,  0.0252   , 0.0118  ]; % cas (v) avec f = 2*pi^2*sin(pi*x)*sin(pi*y)
%relerrL2 = [  0.2730  ,  0.0693  ,  0.0159  ,  0.0074]; % cas (iv) avec f = 1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot norme L2
    fig1 = figure();
    loglog(1./epsilon,relerrL2,'-.+r','MarkerSize', 16, 'linewidth' , 2)
    xlabel('$1/\varepsilon$','Interpreter','latex','FontSize',22);
    ylabel('L2 relative error','FontSize',22);
    set(gca, 'FontSize', 18)
    %legend('err ini-ex','err ini-hom')
    ylim([1e-3 1.2])
    grid on;
    
    text(epsinv(4),relerrL2(4),['[' num2str(epsinv(4)) ',' num2str(relerrL2(4)) ']'],'FontSize',12)
    text(epsinv(3),relerrL2(3),['[' num2str(epsinv(3)) ',' num2str(relerrL2(3)) ']'],'FontSize',12)
    
% plot semi-norme H1
    fig2 = figure();
    loglog(1./epsilon,relerrH1,'-.+b','MarkerSize', 16, 'linewidth' , 2)
    xlabel('$1/\varepsilon$','Interpreter','latex','FontSize',22);
    ylabel('semi-H1 relative error','FontSize',22);
    set(gca, 'FontSize', 18)
    ylim([1e-3 1.2])
    %legend('err ini-ex','err ini-hom')
    grid on;
    
    %text(epsinv(4),relerrH1(4),[' [' num2str(epsinv(4)) ',' num2str(relerrH1(4)) ']'],'FontSize',12)
    %text(epsinv(3),relerrH1(3),[' [' num2str(epsinv(3)) ',' num2str(relerrH1(3)) ']'],'FontSize',12)
    
