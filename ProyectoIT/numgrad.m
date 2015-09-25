% numgrad.m calcula el gradiente numérico de una función en un punto.
% Argumentos de entrada:
% f: función de la que se obtendrá el gradiente.
% x: punto en el que se calculará el gradiente de f.
% h: argumento opcional de entrada que establece con qué precisión se hará
% la derivación.
function[grad] = numgrad(f,x,h)

% Si no se da un valor para h se le asignará el valor de 1.0e-8.
if( nargin < 3 )
    h = 1.0e-8;
end

% Recordando que el gradiente es del mismo tamaño que el vector x en el
% que se evalua se obtiene el tamaño y se usa la matriz identidad para
% poder hacer los cálculos.
n = length(x);
I = eye(n);
grad = zeros(n,1);

% Se calcula el gradiente.
for k = 1:n
    grad(k) = (feval(f, x+h*I(:,k))-feval(f, x-h*I(:,k)))/(2*h);
end

end
