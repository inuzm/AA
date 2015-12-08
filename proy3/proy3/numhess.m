% numhess.m calcula la matriz hessiana num?rica de una funci?n en un punto.
% Argumentos de entrada:
% f: funci?n de la que se obtendr? el gradiente.
% x: punto en el que se calcular? el gradiente de f.
% h: argumento opcional de entrada que establece con qu? precisi?n se har?
% la derivaci?n.
function[H] = numhess(f,x,h)

% Si no se da un valor para h se le asignar? el valor de 1.0e-4.
if( nargin < 3 )
    h = 1.0e-4;
end

% Recordando que las dimensiones de la matriz hessiana son del mismo tama?o
% que el vector x en el que se evalua se obtiene el tama?o y se usa la
% matriz identidad para poder hacer los c?lculos.
n = length(x);
I = eye(n);
H = zeros(n);

% Se calcula la matriz hessiana.
for i = 1:n
    for j = 1:n
        H(i,j) = (feval(f, x+h*I(:,i)+h*I(:,j))-...
            feval(f, x-h*I(:,i)+h*I(:,j))-feval(f, x+h*I(:,i)-h*I(:,j))+...
            feval(f, x-h*I(:,i)-h*I(:,j)))/(4*h^2);
    end
end

end