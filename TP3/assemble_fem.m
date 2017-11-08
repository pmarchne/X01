function [KK,MM,PP,LL,Coorneu,Numtri] = assemble_fem(epsilon,typeA,BC,cellule,nom_maillage)

[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);

%% calcul des matrices EF

% declarations
% ------------
    KK = sparse(Nbpt,Nbpt); % matrice de rigidite
    MM = sparse(Nbpt,Nbpt); % matrice de rigidite
    LL = zeros(Nbpt,1);     % vecteur second membre
    LLcel1 = zeros(Nbpt,1);    % vecteur second membre du probleme de cellule 1
    LLcel2 = zeros(Nbpt,1);    % vecteur second membre du probleme de cellule 2
    
% boucle sur les triangles
% ------------------------
    for l=1:Nbtri
  % Coordonnees des sommets du triangles
    S1 = Coorneu(Numtri(l,1),:);
    S2 = Coorneu(Numtri(l,2),:);
    S3 = Coorneu(Numtri(l,3),:);
  % calcul des matrices elementaires du triangle l 
  
    [Kel,fcell1,fcell2]=matK_elem(S1, S2, S3, epsilon, typeA);
    Mel=matM_elem(S1, S2, S3);
    for i=1:3
          I = Numtri(l,i);
        for j=1:3
              J = Numtri(l,j);
            KK(I,J) = KK(I,J) + Kel(i,j);
            MM(I,J) = MM(I,J) + Mel(i,j);
        end % for j
        % assemblage des seconds membres des problemes de cellule
        LLcel1(I,1) = LLcel1(I,1) + fcell1(i);
        LLcel2(I,1) = LLcel2(I,1) + fcell2(i);
    end %  for i

    end % for l

%% Calcul du second membre L du probleme initial
% -------------------------

    FF = f(Coorneu(:,1),Coorneu(:,2),epsilon,typeA,'dirichlet');
    LL = MM * FF;
    
    LLcel1 = -LLcel1; LLcel2 = -LLcel2;
    

    
%% Constructions des matrices pour la pseudo elimination


% Construction de la matrice de projection Dirichlet
if strcmp(BC,'dirichlet')
    % extraction des noeuds avec reference nulle
    Ref0 = Refneu(Refneu(:)==0,:);    
    Nbzero = size(Ref0,1);
    % initialisation de la matrice de projection
    PP = zeros(Nbzero,Nbpt);
    PP(1:end,(Nbpt-Nbzero+1):end)=eye(Nbzero,Nbzero);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
elseif strcmp(BC,'periodique') 
% Construction de la matrice de projection Periodique
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
    PP = [bloc1 ; bloc2 ; bloc3];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
end 

if strcmp(cellule,'yes')
    LL = [LLcel1,LLcel2];
end
    
    %FF_0 = PP*FF;
end
