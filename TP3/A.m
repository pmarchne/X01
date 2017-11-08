function val = A(x,y,epsilon,type)
global Aeff
if(type == 0)
    val = Aeff;
elseif (type == 1)
    val = eye(2);
elseif (type == 2)
    val = [1 0; 0 2];
elseif (type == 3)
    val = [2+sin(2*pi/epsilon*x) 0 ; 0 4];
elseif (type == 4)
    val = [2+sin(2*pi/epsilon*x) 0 ; 0 4+sin(2*pi/epsilon*x)];
elseif (type == 5)
    val = (sin(2*pi*epsilon*x) + 2)*(sin(2*pi*epsilon*y) + 4)*eye(2);
end

end
