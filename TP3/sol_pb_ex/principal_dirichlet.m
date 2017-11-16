clear all;
close all;
clc;
addpath('/home/philippe/Documents/MATLAB/X01/TP3/sol_pb_ex/Mesh')
% lecture du maillage et affichage
% ---------------------------------
nom_maillages = [{'geomCarre02.msh'},{'geomCarre01.msh'},{'geomCarre005.msh'},{'geomCarre002.msh'},{'geomCarre001.msh'}];
h = [0.2, 0.1, 0.05, 0.02, 0.01]; H = length(h)-1;
epsilon = 1; % Parametre du terme source defini dans A.m
type = 4;

hinv = zeros(1,H); 
relerrL2 = zeros(1,H);
relerrH1 = zeros(1,H);

for k = 1:H % boucle sur les differents maillages
    tic;
    [Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillages{1,k});

% ----------------------
% calcul des matrices EF
% ----------------------

% declarations
% ------------
    KK = sparse(Nbpt,Nbpt); % matrice de rigidite
    MM = sparse(Nbpt,Nbpt); % matrice de rigidite
    LL = zeros(Nbpt,1);     % vecteur second membre
    
% boucle sur les triangles
% ------------------------
    for l=1:Nbtri
  % Coordonnees des sommets du triangles
    S1 = Coorneu(Numtri(l,1),:);
    S2 = Coorneu(Numtri(l,2),:);
    S3 = Coorneu(Numtri(l,3),:);
  % calcul des matrices elementaires du triangle l 
  
    [Kel]=matK_elem(S1, S2, S3, epsilon, type);
    Mel=matM_elem(S1, S2, S3);
    for i=1:3
          I = Numtri(l,i);
        for j=1:3
              J = Numtri(l,j);
            KK(I,J) = KK(I,J) + Kel(i,j);
            MM(I,J) = MM(I,J) + Mel(i,j);
        end % for j
    end %  for i

    end % for l

% Calcul du second membre L
% -------------------------

    FF = f(Coorneu(:,1),Coorneu(:,2),epsilon,type);
    LL = MM * FF;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction de la matrice de projection Dirichlet
% ---------------------------

    % extraction des noeuds avec reference nulle
    Ref0 = Refneu(Refneu(:)==0,:);    
    Nbzero = size(Ref0,1);
    % initialisation de la matrice de projection
    PP = zeros(Nbzero,Nbpt);
    PP(1:end,(Nbpt-Nbzero+1):end)=eye(Nbzero,Nbzero);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
    FF_0 = PP*FF;
    LL_0 = PP*LL;  
    
% inversion
% ----------
    AA = KK;
      
    AA_0 = PP*AA*PP'; 
    UU_0 = AA_0 \ LL_0;
    UU = PP'*UU_0;
    UU_exact = sin(pi*Coorneu(:,1)).*sin(pi*Coorneu(:,2));  
% visualisation
% -------------
    if(k==H)
    affiche(abs(UU-UU_exact), Numtri, Coorneu, sprintf('Dirichlet - %s', nom_maillages{1,k}));
    end
    validation = 'oui';
% validation
% ----------
    if strcmp(validation,'oui')
    
% Calcul de l erreur L2
    err = UU_exact - UU;
    
    hinv(k) = 1/h(k);
    errL2 = sqrt((err)'*MM*err);
    normL2 = sqrt(UU_exact'*MM*UU_exact);
    relerrL2(k) = errL2/normL2;
    
% Calcul de la semi norme H1
    
    errH1 = sqrt((err)'*KK*err);
    normH1 = sqrt(UU_exact'*KK*UU_exact);
    relerrH1(k) = errH1/normH1;
    end
    toc;
end
%plot norme L2
    fig2 = figure(2);
    loglog(hinv,relerrL2,'+-.r','MarkerSize', 16)
    xlabel('1/h','FontSize',22);
    ylabel('L2 relative error','FontSize',22);
    set(gca, 'FontSize', 18)
    xlim([1 300])
    grid on;
    hold on;
    
% plot semi-norme H1
    fig3 = figure(3);
    loglog(hinv,relerrH1,'b+-.','MarkerSize', 16)
    xlabel('1/h','FontSize',22);
    ylabel('semi-H1 relative error','FontSize',22);
    set(gca, 'FontSize', 18)
    xlim([1 300])
    grid on;
    hold on;
