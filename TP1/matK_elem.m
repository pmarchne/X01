function [Kel_new] = matK_elem(S1, S2, S3, psi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mat_elem :
% calcul la matrices de raideur elementaire en P1 lagrange
%
% SYNOPSIS [Kel] = mat_elem(S1, S2, S3)
%          
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle 
%                      (vecteurs reels 1x2)
%
% OUTPUT - Kel matrice de raideur elementaire (matrice 3x3)
%
% NOTE (1) le calcul est exacte (pas de condensation de masse)
%      (2) calcul direct a partir des formules donnees par 
%          les coordonnees barycentriques 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% Transformation F_l du triangle de reference vers le triangle physique
F_l = @(xref,yref) [(x3-x1)*xref + (x2-x1)*yref + x1 , (y3-y1)*xref + (y2-y1)*yref + y1];
% Jacobienne de F_l
Jacobian = [ x3-x1 , x2-x1 ; y3-y1 , y2-y1];

% coordonees des points pour la quadrature sur le triangle physique
X1q = F_l(1/3,1/3);
X2q = F_l(1/5,1/5);
X3q = F_l(1/5,3/5);
X4q = F_l(3/5,1/5);
omega1 = -9/32; omega234 = 25/96;

% quadrature du terme source
A_quad = omega1*A(X1q(1),X1q(2),psi) + omega234*(A(X2q(1),X2q(2),psi)+A(X3q(1),X3q(2),psi)+A(X4q(1),X4q(2),psi));

% les 3 normales a l'arete opposees (de la longueur de l'arete)
norm = zeros(3, 2);
norm(1, :) = [y2-y3, x3-x2];
norm(2, :) = [y3-y1, x1-x3];
norm(3, :) = [y1-y2, x2-x1];
% les 3 normales a l'arete opposees (de la longueur de l'arete) sur le triangle de reference
norm_ref = zeros(3, 2);
norm_ref = [-1,-1;0,1;1,0];

% D est, au signe pres, deux fois l'aire du triangle
D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
if (abs(D) <= eps) 
  error('l aire d un triangle est nulle!!!'); 
end;

% calcul de la matrice de raideur
% -------------------------------
Kel = zeros(3,3);
Kel_new = zeros(3,3);
for i=1:3
  for j=1:3
    Kel(i,j) = norm(i,1)*norm(j,1) + norm(i,2)*norm(j,2);
    % nouvelle matrice de rigidite elementaire
    Kel_new(i,j) = ( (A_quad*inv(Jacobian')*(norm_ref(i,:))')' * (inv(Jacobian')*norm_ref(j,:)') ) * abs(det(Jacobian));
  end; % j
end; % i

Kel = 1/(2*abs(D)) * Kel;
Kel_new = Kel_new;
% retourner Kel ou Kel_new selon le probleme sans ou avec terme source
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
