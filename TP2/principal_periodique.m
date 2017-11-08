% =====================================================
%
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour l'equation de Laplace suivante, avec conditions periodiques
% sur le maillage nom_maillage.msh
%
% | -div(A grad u ) u + u= f,   dans \Omega
% |         u periodique,   sur le bord
%
% =====================================================
clear all;
close all;
%clc;
addpath('/home/philippe/Documents/MATLAB/X01/TP2/Mesh')
% lecture du maillage et affichage
% ---------------------------------
nom_maillages = [{'geomCarre02.msh'},{'geomCarre01.msh'},{'geomCarre005.msh'},{'geomCarre002.msh'},{'geomCarre001.msh'}];
h = [0.2, 0.1, 0.05, 0.02, 0.01]; H = length(h)-2;
psi = 0; % Parametre du terme source defini dans A.m

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
    FF = f(Coorneu(:,1),Coorneu(:,2),psi,'periodique');
    LL = MM * FF;
    
% Construction de la matrice de projection
% ---------------------------
    % noeuds du coin
    Refneu(1:4,:) = 5; % nouvelle reference pour differentier les noeuds du coin 
    RefCoin = Refneu(Refneu(:)==5,:);
    % noeuds en y = 0
    Ref1 = Refneu(Refneu(:)==1,:);
    % noeuds en x = 2
    Ref2 = Refneu(Refneu(:)==2,:);
    % noeuds en y = 2 
    Ref3 = Refneu(Refneu(:)==3,:);
    % noeuds en x = 0
    Ref4 = Refneu(Refneu(:)==4,:);
    % noeuds interieurs
    Ref0 = Refneu(Refneu(:)==0,:);
    Nbzero = size(Ref0,1);
    
    % prise en compte des noeuds du coin
    Ncoin=4;
    MatCoin = 1/4*ones(1,Ncoin); % vecteur pour traiter les noeuds du coin
    
    % prise en compte des noeuds du bord sans les coins
    Nbord = (Nbaretes-4)/2; % ne marche pas si le nombre d'aretes est impair
    MatBord = 1/2*eye(Nbord); 
    MatBord2 = [flip(1/2*eye(Nbord/2)),zeros(Nbord/2,Nbord/2); zeros(Nbord/2,Nbord/2), flip(1/2*eye(Nbord/2))];
    
    MatInt = eye(Nbzero); % matrice qui prend en compte les noeuds de reference 0 a l'interieur du domaine
    sum = Nbzero + Ncoin + (Nbaretes-4);
    
    % assemblage de la matrice de projection
    bloc1 = [MatCoin , zeros(1,sum-4)];
    bloc2 = [zeros(Nbord,4), MatBord , MatBord2, zeros(Nbord,Nbzero)];
    bloc3 = [zeros(Nbzero,sum-Nbzero) , MatInt];
    PP_per = [bloc1 ; bloc2 ; bloc3];
    
    FF_per = PP_per*FF;
    LL_per = PP_per*LL;  
% inversion
% ----------    
    AA = MM + KK;
    AA_per = PP_per*AA*PP_per'; 
    UU_per = AA_per \ LL_per;
    UU = PP_per'*UU_per;
% visualisation
% -------------
    if(k==H)
    affiche(abs(UU), Numtri, Coorneu, sprintf('Periodique - %s', nom_maillages{1,k}));
    end
    validation = 'oui';
% validation
% ----------
    if strcmp(validation,'oui')
    UU_exact = cos(pi*Coorneu(:,1)).*cos(2*pi*Coorneu(:,2));
% Calcul de l erreur L2
    err = UU_exact - UU;
    errL2 = sqrt((err)'*MM*err);
    normL2 = sqrt(UU_exact'*MM*UU_exact);
%plot norme L2
    fig2 = figure(2);
    loglog((1/h(k)),(errL2/normL2),'r+','MarkerSize', 14)
    xlabel('1/h','FontSize',20);
    ylabel('L2 relative error','FontSize',20);
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
    ylabel('semi-H1 relative error','FontSize',20);
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
%cd('/home/philippe/Documents/MATLAB/X01/TP2/Fig')
%print(figure(1),'solution_per002_psi16','-dpng')
% print(fig2,'L2error_per_psi0','-depsc')
% print(fig3,'H1error_per_psi0','-depsc')

