function[x, k, t] = GC(A, b, x, tol, maxiter)
%
% La función gradiente conjugado obtiene la solución de un
% sistema de ecuaciones de la forma:
% A * x = b
% con el algoritmo del gradiente conjugado que se encuentra
% en el libro de Nocedal.
%
% Argumentos de entrada obligatorios:
% A: la matriz que representa las ecuaciones lineales.
% b: el vector que resultará de A * x.
% Argumentos de entrada opcionales:
% x: el punto inicial para el método iterativo.
% tol: la tolerancia que se usará para la condición de
% paro.
% maxiter: el máximo de iteraciones que llevará a cabo el
% método.
%
% Argumentos de salida:
% x: la aproximación al vector solución del problema.
% k: las iteraciones en las cuales se alcanzó la
% aproximación de salida.

% Se revisa que la matriz sea cuadrada. En caso de no serlo
% se dará un mensaje de error.
[m, n] = size(A);

if( m ~= n )
    error('La matriz debe ser una matriz cuadrada.')
    return;
end

% Si alguno de los argumentos opcionales de entrada no se
% dio se le dará un valor predeterminado a cada uno de los
% argumentos faltantes.

if( nargin < 5 )
    maxiter = 2 * n; 
    if( nargin < 4 )
        tol = 1.0e-8;
        if( nargin < 3 )
            x = zeros(n, 1);
        end
    end
end

% Se asignan los valores de r0, p0 y se modifica la
% tolerancia para que sea la tolerancia inicial multiplicada
% por la norma dos de r0.
tic;
r = A * x - b;
p = -r;
k = 0;

% Se lleva a cabo el método iterativo:

while( norm(r) > tol && k < maxiter )
    alpha = (r' * r) / (p' * A * p);
    x = x + alpha * p;
    rtemp = r;
    r = r + alpha * A * p;
    beta = (r' * r) / (rtemp' * rtemp);
    p = -r + beta * p;
    k = k + 1;
end
t = toc;    
end
