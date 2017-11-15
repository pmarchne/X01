%% Begin routine: FEM routine with P1 elements for the heat equation homogenization problem
clear all;
close all;
clc;

% change directory
cd('/home/philippe/Documents/MATLAB/X01/TP3');

% load the mesh
addpath('/home/philippe/Documents/MATLAB/X01/TP3/Mesh')
mesh_name = [{'geomCarre02.msh'},{'geomCarre01.msh'},{'geomCarre005.msh'},{'geomCarre002.msh'},{'geomCarre001.msh'}];
h = [0.2, 0.1, 0.05, 0.02, 0.01]; % size of the meshes with decreasing order 
num_mesh = 4; % choose the mesh to use
h =  h(num_mesh);
% -----------------------------------------------------------------------------------------------------------------------

% initialize the value(s) for epsilon
% Parametre de la taille de la microstructure du problème: On fait varier
% le nombre d'élements par période en changeant epsilon: la période doit
% etre choisie assez grande par rapport au pas de maillage
%epsilon = [2,0.5,0.1,0.05];
epsilon = [4,1,0.4,0.2,0.08,0.03];
Nb_epsilon = length(epsilon);
% -----------------------------------------------------------------------------------------------------------------------

% le choix du terme source se trouve dans la routine A.m
global Aeff; Aeff = zeros(2,2);
Anum = 4; % choisir Anum et changer f.m pour typeA == 0
% -----------------------------------------------------------------------------------------------------------------------

% initialize the errors
epsinv = zeros(1,Nb_epsilon);
relerrL2 = zeros(1,Nb_epsilon);
relerrH1 = zeros(1,Nb_epsilon);
% -----------------------------------------------------------------------------------------------------------------------

%% Problem solving
for i = 1:Nb_epsilon % boucle sur epsilon
    
disp(['mesh size = ' num2str(h) ]);
disp(['epsilon = ' num2str(epsilon(i))]);


% inversion et resolution des problemes
% ----------

% 1) probleme initial, condition de dirichlet
tic;
typeBC = 'dirichlet'; cellule = 'no';
[KK_ini,~,PP_0,LL_ini,Coorneu,Numtri] = assemble_fem(epsilon(i),Anum,typeBC,cellule,mesh_name{1,num_mesh});
AA = KK_ini;
LL_0 = PP_0*LL_ini; 
AA_0 = PP_0*AA*PP_0'; 
UU_0 = AA_0 \ LL_0;

UU_ini = PP_0'*UU_0; % solution of the initial Poisson problem
time1 = toc;
disp(['solving time initial problem = ' num2str(time1) ' sec']);


% 2) problemes de cellules, conditions periodiques
tic;
typeBC = 'periodique'; cellule = 'yes';
[KK_cel,MM_cel,PP_per,LL_cel,Coorneu_cel,Numtri_cel] = assemble_fem(epsilon(i),Anum,typeBC,cellule,'geomCarre_per_002.msh'); 
% on utilise le maillage periodique pour resoudre le probleme de cellule

eta = 5; % constante de pénalisation
AA_cel = KK_cel + eta*MM_cel;
AA_cel_per = PP_per*AA_cel*PP_per';

% resolution des problèmes de cellule
LL_cel_per = PP_per*LL_cel;
UU_cel_per = AA_cel_per \ LL_cel_per;
UU_cel = PP_per'*UU_cel_per;

% construction du tenseur homogénéisé
Aeff = (Coorneu_cel+UU_cel)'*KK_cel*(Coorneu_cel+UU_cel);

time2 = toc;
disp(['solving time cell problems = ' num2str(time2) ' sec']);

% 3) Probleme homogeneisé, conditions dirichlet
tic;
typeBC = 'dirichlet'; cellule = 'no';
[KK_hom,MM_hom,~,LL_hom,~,~] = assemble_fem(epsilon(i),0,typeBC,cellule,mesh_name{1,num_mesh});

AA_hom = KK_hom;
LL_0 = PP_0*LL_hom; 
AA_0 = PP_0*AA_hom*PP_0'; 
UU_0 = AA_0 \ LL_0;

% solution du problème homogénéisé
UU_hom = PP_0'*UU_0;

time3 = toc;
disp(['solving time homogenized problem = ' num2str(time3) ' sec']);

% ------------------------------------------------------------
% 4) calcul des erreurs
    tic;
    %hinv(i) = 1/h(i);
    epsinv(i) = 1/epsilon(i);
    err = abs(UU_ini - UU_hom);
    %U1 = (UU_cel(:,1) + UU_0)*KK + (UU_cel(:,2) + UU_0)*KK;
    %err_new = abs(UU_ini - UU_hom - epsilon(i)*U1); erreur en racine de
    %epsilon
    %[KK_err,MM_err,~,~,~,~] = assemble_fem(epsilon(i),1,'dirichlet','no',mesh_name{1,num_mesh});
    
    % Calcul de erreur L2
    normL2hom = sqrt(UU_hom'*MM_hom*UU_hom);
    errL2 = sqrt((err)'*MM_hom*err);
    relerrL2(i) = errL2/normL2hom; 
    
    %calcul norme H1
    normH1hom = sqrt(UU_hom'*KK_hom*UU_hom);
    errH1 = sqrt((err)'*KK_hom*err);
    relerrH1(i) = errH1/normH1hom;
    
    time4 = toc;
    disp(['solving time errors = ' num2str(time4) ' sec']);
    disp(['total time = ' num2str(time1+time2+time3+time4) ' sec']);
    fprintf('\n');
    
    % visualisation
    affiche(UU_ini, Numtri, Coorneu, sprintf('ini, eps = %g', epsilon(i)));
    affiche(UU_cel(:,1), Numtri_cel, Coorneu_cel, sprintf('Cel 1, eps = %g', epsilon(i)));
    affiche(UU_cel(:,2), Numtri_cel, Coorneu_cel, sprintf('Cel 2, eps = %g', epsilon(i)));
    affiche(UU_hom, Numtri, Coorneu, sprintf('hom, eps = %g', epsilon(i)));
end
%% Plot figures
    
%plot norme L2
    fig2 = figure();
    loglog(epsinv,relerrL2,'+-.k','MarkerSize', 16)
    xlabel('1/eps','FontSize',22);
    ylabel('L2 relative error','FontSize',18);
    set(gca, 'FontSize', 18)
    %legend('err ini-ex','err ini-hom')
    %xlim([1 300])
    grid on;
    
% plot semi-norme H1
    fig3 = figure();
    loglog(epsinv,relerrH1,'k+-.','MarkerSize', 16)
    xlabel('1/eps','FontSize',22);
    ylabel('semi-H1 relative error','FontSize',18);
    set(gca, 'FontSize', 18)
    %xlim([1 300])
    %legend('err ini-ex','err ini-hom')
    grid on;
    
%% save figures

 %cd('/home/philippe/Documents/MATLAB/X01/TP2/FigTP3')
 %print(figure(1),'case5','-dpng')
 %print(fig2,'L2_case5','-depsc')
 %print(fig3,'semiH1_case5','-depsc')
 
%x0=50;
%y0=50;
%width=700;
%height=600;
%set(gcf,'units','points','position',[x0,y0,width,height])
