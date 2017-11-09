%% Begin routine: FEM routine with P1 elements for the heat equation homogenization problem
clear all;
close all;
clc;

% change directory
cd('/home/philippe/Documents/MATLAB/X01/TP3');

% load the meshes
addpath('/home/philippe/Documents/MATLAB/X01/TP3/Mesh')
mesh_name = [{'geomCarre02.msh'},{'geomCarre01.msh'},{'geomCarre005.msh'},{'geomCarre002.msh'},{'geomCarre001.msh'}];
h = [0.2, 0.1, 0.05, 0.02, 0.01]; % size of the meshes with decreasing order 
%Nb_mesh = length(h)-2; % run the routine for the 3 first meshes
num_mesh = 3; % run the main program for mesh size 0.05
h =  h(num_mesh);
% -----------------------------------------------------------------------------------------------------------------------

% initialize the value(s) for epsilon
% Parametre de la taille de la microstructure périodique du problème
epsilon = [10*h];
Nb_epsilon = length(epsilon);
% h should be small compared to epsilon
% -----------------------------------------------------------------------------------------------------------------------

% le choix du terme source se trouve dans la routine A.m
global Aeff; Aeff = zeros(2,2);
Anum = 3; % si on change Anum, on doit changer f.m aussi dans typeA=0
% -----------------------------------------------------------------------------------------------------------------------

% initialize the errors
%hinv = zeros(1,Nb_epsilon); % parametric study wrt the mesh size
epsinv = zeros(1,Nb_epsilon);
relerrL2_ini = zeros(1,Nb_epsilon);
relerrL2_hom = zeros(1,Nb_epsilon);
relerrH1_ini = zeros(1,Nb_epsilon);
relerrH1_hom = zeros(1,Nb_epsilon);
% -----------------------------------------------------------------------------------------------------------------------

%% Problem solving
for i = 1:Nb_epsilon % boucle sur epsilon
    
disp(['mesh size = ' num2str(h) ]);
disp(['epsilon = ' num2str(epsilon(i))]);
tic;

% inversion et resolution des problemes
% ----------

% 1) probleme initial, condition de dirichlet
typeBC = 'dirichlet'; cellule = 'no'; typeA = Anum;
[KK_ini,~,PP_0,LL_ini,Coorneu,Numtri] = assemble_fem(epsilon(i),typeA,typeBC,cellule,mesh_name{1,num_mesh});
AA = KK_ini;
LL_0 = PP_0*LL_ini; 
AA_0 = PP_0*AA*PP_0'; 
UU_0 = AA_0 \ LL_0;

UU_ini = PP_0'*UU_0; % solution of the initial Poisson problem

% 2) problemes de cellules, conditions periodiques
typeBC = 'periodique'; cellule = 'yes'; typeA = Anum;
[KK,MM,PP_per,LL_cel,~,~] = assemble_fem(epsilon(i),typeA,typeBC,cellule,mesh_name{1,num_mesh}); 

eta = 0.1; % constante de pénalisation
AA_cel = KK + eta*MM;
AA_cel_per = PP_per*AA_cel*PP_per';

LL_cel1 = LL_cel(:,1); LL_cel2 = LL_cel(:,2);
LL_cel1_per = PP_per*LL_cel1;
LL_cel2_per = PP_per*LL_cel2;
UU_cel1_per = AA_cel_per \ LL_cel1_per;
UU_cel2_per = AA_cel_per \ LL_cel2_per;
    
% solution des 2 problèmes de cellule
UU_cel1 = PP_per'*UU_cel1_per;
UU_cel2 = PP_per'*UU_cel2_per;
UU_cel = [UU_cel1,UU_cel2];

% construction du tenseur homogénéisé
Aeff = (Coorneu+UU_cel)'*KK*(Coorneu+UU_cel); Aeff = Aeff/4;

% 3) Probleme homogeneisé, conditions dirichlet
typeBC = 'dirichlet'; cellule = 'no'; typeA = 0;
[KK_hom,~,PP_0,LL_hom,~,~] = assemble_fem(epsilon(i),typeA,typeBC,cellule,mesh_name{1,num_mesh});
AA = KK_hom;
LL_0 = PP_0*LL_hom; 
AA_0 = PP_0*AA*PP_0'; 
UU_0 = AA_0 \ LL_0;

% solution du problème homogénéisé
UU_hom = PP_0'*UU_0;
time = toc;
disp(['elapsed time = ' num2str(time) ' sec']);

% 4) exact solution
UU_exact = sin(pi*Coorneu(:,1)).*sin(pi*Coorneu(:,2));


% ------------------------------------------------------------

% 5) calcul des erreurs
    %hinv(i) = 1/h(i);
    epsinv(i) = 1/epsilon(i);
    err_iniex = UU_ini - UU_exact;
    err_inihom = UU_ini - UU_hom;
    
    % calcul des normes L2
    normL2ex = sqrt(UU_exact'*MM*UU_exact);
    normL2hom = sqrt(UU_hom'*MM*UU_hom);
    
    %calcul des normes H1
    normH1ex = sqrt(UU_exact'*KK_ini*UU_exact);
    normH1hom = sqrt(UU_hom'*KK_ini*UU_hom);
    
% Calcul de l erreur L2
    % erreur du probleme initial
    errL2_ini = sqrt((err_iniex)'*MM*err_iniex);
    relerrL2_ini(i) = errL2_ini/normL2ex;
    % erreur pour le probleme homogénéisé
    errL2_hom = sqrt((err_inihom)'*MM*err_inihom);
    relerrL2_hom(i) = errL2_hom/normL2hom;   
    
% Calcul de la semi norme H1
    % erreur du probleme initial
    errH1_ini = sqrt((err_iniex)'*KK_ini*err_iniex);
    relerrH1_ini(i) = errH1_ini/normH1ex;
    % erreur pour le probleme homogénéisé
    errH1_hom = sqrt((err_inihom)'*KK_ini*err_inihom);
    relerrH1_hom(i) = errH1_hom/normH1hom;

end


%% Plot figures

% visualisation
    affiche(UU_ini, Numtri, Coorneu, sprintf('Probleme ini'));
    affiche(UU_cel1, Numtri, Coorneu, sprintf('Cellule 1'));
    affiche(UU_cel2, Numtri, Coorneu, sprintf('Cellule 2'));
    affiche(UU_hom, Numtri, Coorneu, sprintf('Probleme homogeneise'));
    
%plot norme L2
    fig2 = figure();
    loglog(epsinv,relerrL2_ini,'+-.r','MarkerSize', 16)
    hold on;
    loglog(epsinv,relerrL2_hom,'+-.k','MarkerSize', 16)
    xlabel('1/eps','FontSize',22);
    ylabel('L2 relative error','FontSize',18);
    set(gca, 'FontSize', 18)
    legend('err ini-ex','err ini-hom')
    %xlim([1 300])
    grid on;
    
% plot semi-norme H1
    fig3 = figure();
    loglog(epsinv,relerrH1_ini,'b+-.','MarkerSize', 16)
    hold on;
    loglog(epsinv,relerrH1_hom,'k+-.','MarkerSize', 16)
    xlabel('1/eps','FontSize',22);
    ylabel('semi-H1 relative error','FontSize',18);
    set(gca, 'FontSize', 18)
    %xlim([1 300])
    legend('err ini-ex','err ini-hom')
    grid on;
    
%% save figures

 %cd('/home/philippe/Documents/MATLAB/X01/TP2/FigTP3')
 %print(figure(1),'case5','-dpng')
 %print(fig2,'L2_case5','-depsc')
 %print(fig3,'semiH1_case5','-depsc')
