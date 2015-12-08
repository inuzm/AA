% numgrad.m calcula el gradiente num?rico de una funci?n en un punto.
% Argumentos de entrada:
% f: funci?n de la que se obtendr? el gradiente.
% x: punto en el que se calcular? el gradiente de f.
% h: argumento opcional de entrada que establece con qu? precisi?n se har?
% la derivaci?n.
function[grad] = numgrad(f,x,h)

% Si no se da un valor para h se le asignar? el valor de 1.0e-8.
if( nargin < 3 )
    h = 1.0e-8;
end

% Recordando que el gradiente es del mismo tama?o que el vector x en el
% que se evalua se obtiene el tama?o y se usa la matriz identidad para
% poder hacer los c?lculos.
n = length(x);
I = eye(n);
grad = zeros(n,1);

% Se calcula el gradiente.
for k = 1:n
    grad(k) = (feval(f, x+h*I(:,k))-feval(f, x-h*I(:,k)))/(2*h);
end

end