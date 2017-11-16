function val = f(x,y,epsilon,type)

% Probleme de Dirichlet
% calcul de f a partir de la donnee u donne:

    if (type == 1)
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
          (-1/epsilon*(0.5*sin(2*pi/epsilon*x) + 1).*sin(pi*x).*cos(pi*y).*cos(2*pi/epsilon*y) - ...
           1/epsilon*(sin(2*pi/epsilon*y) + 4).*sin(pi*y).*cos(pi*x).*cos(2*pi/epsilon*x)*0.5 + ...
           (0.5*sin(2*pi/epsilon*x) + 1).*(sin(2*pi/epsilon*y) + 4).*sin(pi*x).*sin(pi*y));
    end
    

end
