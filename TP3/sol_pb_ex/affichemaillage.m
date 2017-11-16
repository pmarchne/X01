function affichemaillage(nom_maillage, titre);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% affichemaillage: 
% pour visualiser un maillage triangulaire 2D
%
% SYNOPSIS affichemaillage(nom_maillage, titre)
%          
% INPUT  * nom_maillage : racine du fichier de maillage .msh (string) 
%        * titre (optionel) un titre (string)
%
% OUTPUT une fenetre graphique
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% control on the input args
if (nargin<2), titre = ''; end;

% lecture du fichier nom_maillage.amdba
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri]=lecture_msh(nom_maillage);

%visualisation du maillage
figure;
hold on

% maillage
trimesh(Numtri,Coorneu(:,1),Coorneu(:,2),zeros(Nbpt,1),'Edgecolor','red');
view(2);
axis('equal');
xlabel('x','FontSize',24);
ylabel('y','FontSize',24);
grid on;
%xt = get(gca, 'XTick');
set(gca, 'FontSize', 24)
% ajouter eventuellement un titre
%title(titre);

hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
