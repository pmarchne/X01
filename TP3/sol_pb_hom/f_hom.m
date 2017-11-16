function val = f_hom(x,y)

% Second membre pour Probleme homogénéisé
       val = 2*pi^2.*sin(pi*x).*sin(pi*y);
       %val = ones(length(x),1);
       %val = x.*y;
end
