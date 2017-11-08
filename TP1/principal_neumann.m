% =====================================================
%
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour l'equation de Laplace suivante, avec conditions de
% Neumann sur le maillage nom_maillage.msh
%
% | -Delta  u + u= f,   dans \Omega
% |         du/dn = 0,   sur le bord
%
% =====================================================
clear all;
close all;
%clc;
addpath('/home/philippe/Documents/MATLAB/X01/TP1/Mesh')
% lecture du maillage et affichage
% ---------------------------------
nom_maillages = [{'geomCarre02.msh'},{'geomCarre01.msh'},{'geomCarre005.msh'},{'geomCarre002.msh'},{'geomCarre001.msh'}];
h = [0.2, 0.1, 0.05, 0.02, 0.01]; H = length(h)-2;
psi = 1; % Parametre du terme source defini dans A.m

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
  
    Kel=matK_elem(S1, S2, S3, psi);
    Mel=matM_elem(S1, S2, S3);
    for i=1:3
          I = Numtri(l,i);
        for j=1:3
              J = Numtri(l,j);
            KK(I,J) = KK(I,J) + Kel(i,j);
            MM(I,J) = MM(I,J) + Mel(i,j);
        end % for j
    end %  for i
  % On fait l'assemmblage de la matrice globale et du second membre

    end % for l

% Calcul du second membre L
% -------------------------
    FF = f(Coorneu(:,1),Coorneu(:,2),psi);
    LL = MM * FF;
% inversion
% ----------
    AA = MM + KK; 
    UU = AA \ LL;
    UU_exact = cos(pi*Coorneu(:,1)).*cos(2*pi*Coorneu(:,2));
% visualisation
% -------------
    if(k==H)
    affiche(abs(UU), Numtri, Coorneu, sprintf('Neumann - %s', nom_maillages{1,k}));
    affiche(abs(UU_exact), Numtri, Coorneu, sprintf('Neumann - %s', nom_maillages{1,k}));
    end
    validation = 'oui';
% validation
% ----------
    if strcmp(validation,'oui')
    
% Calcul de l erreur L2
    err = UU_exact - UU;
    errL2 = sqrt((err)'*MM*err);
    normL2 = sqrt(UU_exact'*MM*UU_exact);
%plot norme L2
    fig2 = figure(2);
    loglog((1/h(k)),(errL2/normL2),'r+','MarkerSize', 14)
    xlabel('1/h','FontSize',20);
    ylabel('L2 Normalized error','FontSize',20);
    set(gca, 'FontSize', 16)
    xlim([1 300])
    grid on;
    hold on;
    
% Calcul de la semi norme H1
    errH1 = sqrt((err)'*KK*err);
    normH1 = sqrt(UU_exact'*KK*UU_exact);
% plot semi-norme H1
    fig3 = figure(3);
    loglog((1/h(k)),(errH1/normH1),'b+','MarkerSize', 14)
    xlabel('1/h','FontSize',20);
    ylabel('H1 Normalized error','FontSize',20);
    set(gca, 'FontSize', 16)
    xlim([1 300])
    grid on;
    hold on;
    
    end
    toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% % save figures
% cd('/home/philippe/Documents/MATLAB/X01/TP1/Fig')
% print(figure(1),'solution001_psi2','-dpng')
% print(fig2,'L2error_psi2','-depsc')
% print(fig3,'H1error_psi2','-depsc')
