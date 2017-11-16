function Mel = matM_elem(S1, S2, S3)

% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% D est, au signe pres, deux fois l'aire du triangle
D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
if (abs(D) <= eps) 
  error('l aire d un triangle est nulle!!!'); 
end;

% calcul de la matrice de masse
% -----------------------------
Mel = zeros(3,3);
for i=1:3
	for j=1:3
		Mel(i,j) = 1/12;
	end;
    Mel(i,i) = 1/6; 
end;

Mel = abs(D)/2 * Mel;
end
