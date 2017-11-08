clear all;
close all;
clc;
addpath('/home/philippe/Documents/MATLAB/X01/TP2/Mesh')
% lecture du maillage et affichage
% ---------------------------------
nom_maillages = [{'geomCarre02.msh'},{'geomCarre01.msh'},{'geomCarre005.msh'},{'geomCarre002.msh'},{'geomCarre001.msh'}];
h = [0.2, 0.1, 0.05, 0.02, 0.01]; H = length(h)-2;
epsilon = 0.01; % Parametre du terme source defini dans A.m
type = 3;

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
    LL1 = zeros(Nbpt,1);
    LL2 = zeros(Nbpt,1);
    
% boucle sur les triangles
% ------------------------
    for l=1:Nbtri
  % Coordonnees des sommets du triangles
    S1 = Coorneu(Numtri(l,1),:);
    S2 = Coorneu(Numtri(l,2),:);
    S3 = Coorneu(Numtri(l,3),:);
  % calcul des matrices elementaires du triangle l 
  
    [Kel,fcell1,fcell2]=matK_elem(S1, S2, S3, epsilon, type);
    Mel=matM_elem(S1, S2, S3);
    for i=1:3
          I = Numtri(l,i);
        for j=1:3
              J = Numtri(l,j);
            KK(I,J) = KK(I,J) + Kel(i,j);
            MM(I,J) = MM(I,J) + Mel(i,j);
        end % for j
        LL1(I,1) = LL1(I,1) + fcell1(i);
        LL2(I,1) = LL2(I,1) + fcell2(i);
    end %  for i

    end % for l

% Calcul du second membre L
% -------------------------

    FF = f(Coorneu(:,1),Coorneu(:,2),epsilon,type,'dirichlet');
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
   
% Construction de la matrice de projection Periodique
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    LL1 = -LL1; LL2 = -LL2;
    LL_cel1_per = PP_per*LL1;
    LL_cel2_per = PP_per*LL2;
    
    FF_0 = PP*FF;
    LL_0 = PP*LL;  
    
% inversion
% ----------
    AA = KK;
    eta = 0.05;
    AA_cel = KK + eta*MM;
    
    
    AA_0 = PP*AA*PP'; 
    UU_0 = AA_0 \ LL_0;
    UU = PP'*UU_0;
    UU_exact = sin(pi*Coorneu(:,1)).*sin(pi*Coorneu(:,2));
    
    AA_cel_per = PP_per*AA_cel*PP_per'; 
    UU_cel1_per = AA_cel_per \ LL_cel1_per;
    UU_cel2_per = AA_cel_per \ LL_cel2_per;
    
    UU_cel1 = PP_per'*UU_cel1_per;
    UU_cel2 = PP_per'*UU_cel2_per;
    
    Aeff(1,1) = (Coorneu(:,1)+UU_cel1)'*KK*(Coorneu(:,1)+UU_cel1);
    Aeff(2,1) = (Coorneu(:,1)+UU_cel1)'*KK*(Coorneu(:,2)+UU_cel2);
    Aeff(1,2) = (Coorneu(:,2)+UU_cel2)'*KK*(Coorneu(:,1)+UU_cel1);
    Aeff(2,2) = (Coorneu(:,2)+UU_cel2)'*KK*(Coorneu(:,2)+UU_cel2);
    Aeff = Aeff/4;
    
    
% visualisation
% -------------
    if(k==H)
    affiche(abs(UU), Numtri, Coorneu, sprintf('Dirichlet - %s', nom_maillages{1,k}));
    affiche(abs(UU_cel1), Numtri, Coorneu, sprintf('Cellule 1 - %s', nom_maillages{1,k}));
    affiche(abs(UU_cel2), Numtri, Coorneu, sprintf('Cellule 2 - %s', nom_maillages{1,k}));
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
    
% % save figures
 %cd('/home/philippe/Documents/MATLAB/X01/TP2/FigTP3')
 %print(figure(1),'case5','-dpng')
 %print(fig2,'L2_case5','-depsc')
 %print(fig3,'semiH1_case5','-depsc')
