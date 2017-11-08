function [Kel,fel_cellule1,fel_cellule2] = matK_elem(S1, S2, S3, psi,type)

% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% Transformation F_l du triangle de reference vers le triangle physique
F_l = @(xref,yref) [(x2-x1)*xref + (x3-x1)*yref + x1 , (y2-y1)*xref + (y3-y1)*yref + y1];
% Jacobienne de F_l
Jacobian = [ x2-x1 , x3-x1 ; y2-y1 , y3-y1];

% coordonees des points pour la quadrature sur le triangle physique
X1q = F_l(1/3,1/3);
X2q = F_l(1/5,1/5);
X3q = F_l(1/5,3/5);
X4q = F_l(3/5,1/5);
omega1 = -9/32 ;
omega234 = 25/96;

% quadrature du terme source
A_quad = omega1*A(X1q(1),X1q(2),psi,type) + omega234*(A(X2q(1),X2q(2),psi,type)+A(X3q(1),X3q(2),psi,type)+A(X4q(1),X4q(2),psi,type));

% les 3 normales sur le triangle de reference
norm_ref = zeros(3, 2);
norm_ref = [-1,-1;1,0;0,1];

% calcul de la matrice de raideur
% -------------------------------
Kel = zeros(3,3);
for i=1:3
  for j=1:3
    Kel(i,j) = (     (  A_quad*inv(Jacobian')*(norm_ref(i,:))'  )' * (inv(Jacobian')*norm_ref(j,:)') ) * abs(det(Jacobian));
  end; % j
end; % i

fel_cellule1 = zeros(1,3);
fel_cellule2 = zeros(1,3);
for p=1:3
    fel_cellule1(p) = (     (  A_quad*(norm_ref(2,:))'  )' * (inv(Jacobian')*norm_ref(p,:)') ) * abs(det(Jacobian));
    fel_cellule2(p) = (     (  A_quad*(norm_ref(3,:))'  )' * (inv(Jacobian')*norm_ref(p,:)') ) * abs(det(Jacobian));
end

end
