
%%%%%%%%%%%%%%%% cas (iv) %%%%%%%%%%%%%%%%%%%%%%
epsilon1 = [2,0.5,0.1,0.05];
epsilon2 = [4,1,0.4,0.2,0.08,0.03];
epsilon = [epsilon1 ,epsilon2];
epsinv = 1./epsilon;

relerrH11 = [  0.2979   , 0.1842   , 0.1856  ,  0.1447];
relerrL21 = [  0.2928   , 0.0950   , 0.0415   , 0.0215];


relerrH12 = [0.0637  ,  0.1724  ,  0.1827  ,  0.1964  ,  0.1717  ,  0.0860];
relerrL22 = [  0.0442  ,  0.1082   , 0.0742 ,   0.0605  ,  0.0285  ,  0.0095];

relerrH1 = [ relerrH11, relerrH12];
relerrL2 = [relerrL21 , relerrL22];

epsilon = [2,0.5,0.2,0.1,0.05,0.03];
epsinv = 1./epsilon;
relerrH1 =[0.2979   , 0.1842, 0.1964,0.1856,0.1447, 0.0860];
relerrL2 = [0.2928   , 0.0950, 0.0605,0.0415,0.0215,0.0095 ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot norme L2
    fig1 = figure();
    loglog(1./epsilon,relerrL2,'-.+r','MarkerSize', 16)
    xlabel('1/eps','FontSize',22);
    ylabel('L2 relative error','FontSize',22);
    set(gca, 'FontSize', 18)
    %legend('err ini-ex','err ini-hom')
    ylim([1e-3 1])
    grid on;
    
    text(epsinv(4),relerrL2(4),[' [' num2str(epsinv(4)) ',' num2str(relerrL2(4)) ']'],'FontSize',12)
    text(epsinv(5),relerrL2(5),[' [' num2str(epsinv(5)) ',' num2str(relerrL2(5)) ']'],'FontSize',12)
    
% plot semi-norme H1
    fig2 = figure();
    loglog(1./epsilon,relerrH1,'-.+b','MarkerSize', 16)
    xlabel('1/eps','FontSize',22);
    ylabel('semi-H1 relative error','FontSize',22);
    set(gca, 'FontSize', 18)
    ylim([1e-3 1])
    %legend('err ini-ex','err ini-hom')
    grid on;
    
    text(epsinv(4),relerrH1(4),[' [' num2str(epsinv(4)) ',' num2str(relerrH1(4)) ']'],'FontSize',12)
    text(epsinv(5),relerrH1(5),[' [' num2str(epsinv(5)) ',' num2str(relerrH1(5)) ']'],'FontSize',12)
    