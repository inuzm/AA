% numhess.m calcula la matriz hessiana numérica de una función en un punto.
% Argumentos de entrada:
% f: función de la que se obtendrá el gradiente.
% x: punto en el que se calculará el gradiente de f.
% h: argumento opcional de entrada que establece con qué precisión se hará
% la derivación.
function[H] = numhess(f,x,h)

% Si no se da un valor para h se le asignará el valor de 1.0e-4.
if( nargin < 3 )
    h = 1.0e-4;
end

% Recordando que las dimensiones de la matriz hessiana son del mismo tamaño
% que el vector x en el que se evalua se obtiene el tamaño y se usa la
% matriz identidad para poder hacer los cálculos.
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
