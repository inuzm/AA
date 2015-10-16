function[x, k] = GCPre(A, b, x, M, tol, maxiter)
%
% La función gradiente conjugado prec obtiene la solución
% de un sistema de ecuaciones de la forma:
% A * x = b
% con el algoritmo del gradiente conjugado precondicionado
% que se encuentra en el libro de Nocedal.
%
% Argumentos de entrada obligatorios:
% A: la matriz que representa las ecuaciones lineales.
% b: el vector que resultará de A * x.
% Argumentos de entrada opcionales:
% x: el punto inicial para el método iterativo.
% M: la matriz con la cual se precondicionará usando M' * M.
% tol: la tolerancia que se usará para la condición de
% paro.
% maxiter: el maximo de iteraciones que llevará a cabo el
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

if( nargin < 6 )
    maxiter = 2 * n; 
    if( nargin < 5 )
        tol = 1.0e-8;
        if( nargin < 4 )
	% Si no se diera matriz M entonces se asignará	
	% la matriz identidad, resultando en el algoritmo
	% usual del gradiente conjugado.
            M = eye(m);
            if( nargin < 3)
                x = zeros(n, 1);
            end
        end
    end
end

% Se asignan los valores de r0, p0, y0 y se modifica la
% tolerancia para que sea la tolerancia inicial multiplicada
% por la norma dos de r0.

r = A * x - b;
k = 0;
y = M' \ r;
y = M \ y;
p = -y;

% Se lleva a cabo el método iterativo:

while( norm(r) > tol && k < maxiter )
    alpha = (r' * y) / (p' * A * p);
    x = x + alpha * p;
    ytemp = y;
    rtemp = r;
    r = r + alpha * A * p;
    y = M' \ r;
    y = M \ y;
    beta = (r' * y) / (rtemp' * ytemp);
    p = -y + beta * p;
    k = k + 1;
end
    
end
