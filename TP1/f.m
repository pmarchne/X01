function val = f(x,y,psi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f :
% Evaluation de la fonction second membre.
%
% SYNOPSIS val = f(x,y)
%          
% INPUT * x,y : les 2 coordonnees du point ou on veut evaluer la fonction.
%
% OUTPUT - val: valeur de la fonction sur ce point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Etape de validation
% Sans terme source, le calcul de f a partir de la donnee u donne:
%f0 = @(x,y) (cos(pi*x).*cos(2*pi*y)*(1+10*pi^2));
%val = f0(x,y);


%Avec le terme source A, et en fonction de psi, un calcul formel donne l'expression suivante:
f = @(x,y,psi) (pi^2/psi*sin(pi*x).*sin(pi/psi*y).*cos(2*pi*y).*cos(pi/psi*x) - 2*pi^2/psi*sin(2*pi*y).*sin(pi/psi*x).*cos(pi*x).*cos(pi/psi*y) + ...
    5*pi^2*(sin(pi/psi*x).*sin(pi/psi*y) + 2).*cos(pi*x).*cos(2*pi*y) + cos(pi*x).*cos(2*pi*y));

val = f(x,y,psi);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
