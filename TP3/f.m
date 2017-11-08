function val = f(x,y,epsilon,type,probleme)

% Probleme de Dirichlet
% Sans terme source, le calcul de f a partir de la donnee u donne:
if strcmp(probleme,'dirichlet')
    if(type == 0)
        val  = (4*sqrt(3)+2*sqrt(15))*pi^2*sin(pi*x).*sin(pi*y);
    elseif (type == 1)
        val = 2*pi^2*sin(pi*x).*sin(pi*y);
    elseif (type == 2)
        val = 3*pi^2*sin(pi*x).*sin(pi*y);
    elseif (type == 3)
        val = pi^2*(-2/epsilon*cos(pi*x).*cos(2*pi/epsilon*x) + ...
                  sin(pi*x).*sin(2*pi/epsilon*x) + 6*sin(pi*x)).*sin(pi*y);
    elseif (type == 4)
        val = pi^2*(-2/epsilon*cos(pi*x).*cos(2*pi/epsilon*x) + ...
                2*sin(pi*x).*sin(2*pi/epsilon*x) + 6*sin(pi*x)).*sin(pi*y);
    elseif (type == 5)
        val = 4*pi^2* ...
          (-epsilon*(0.5*sin(2*pi*epsilon*x) + 1).*sin(pi*x).*cos(pi*y).*cos(2*pi*epsilon*y) - ...
           epsilon*(sin(2*pi*epsilon*y) + 4).*sin(pi*y).*cos(pi*x).*cos(2*pi*epsilon*x)*0.5 + ...
           (0.5*sin(2*pi*epsilon*x) + 1).*(sin(2*pi*epsilon*y) + 4).*sin(pi*x).*sin(pi*y));
    end

% Probleme periodique
% Sans terme source, le calcul de f a partir de la donnee u donne:
elseif strcmp(probleme,'periodique')
    f0 = @(x,y) (cos(pi*x).*cos(2*pi*y)*(1+10*pi^2));
    val = f0(x,y);

else disp('error: ecrire dirichlet ou periodique')
    

end
