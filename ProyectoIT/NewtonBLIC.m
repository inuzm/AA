% NewtonBLIC.m resuelve el problema de la forma:
% 		min f(x)
% sin restricciones. Esto lo hace con el paso de Newton e interpolación
% cuadrática.
% Argumentos de entrada:
% func: la función f que se quiere minimizar.
% x: punto inicial para el método.
% tol: la tolerancia para la condición de paro; este argumento es opcional.
% maxiter: el máximo de iteraciones que se harán en el método; este
% argumento es opcional
% Argumentos de salida:
% x: punto en el que la norma del gradiente es menor que la tolerancia.
% k: número de iteraciones en el que se llegó al valor de x.
function[x, k] = NewtonBLIC(func, x, tol, maxiter)

% Si no se dan todos los argumentos entonces se revisa cuáles no se dieron.
% Si no se dio maxiter se le asignará el valor de 10000. Además, si no se
% da tol se le asignará el valor de 1.0e-8.
if( nargin < 4 )
	maxiter = 10000;
	if( nargin < 3 )
		tol = 1.0e-8;
	end
end

% Se calculan el gradiente y la matriz hessiana de la función en el punto
% inicial y se le da el valor de 0 a k.
gf = numgrad(func, x);
Hf = numhess(func, x);
k = 0;

% Mientras no se cumpla la condición de paro, que la norma del gradiente
% sea menor que cierta tolerancia o que se llegue al máximo de iteraciones,
% se realiza la minimización de f con el método de Newton.
while( norm(gf) > tol && k < maxiter )
    p = - Hf \gf;
	alpha = alphaIC(func, x, p);
    x = x + alpha * p;
	gf = numgrad(func, x);
	Hf = numhess(func, x);
	k = k + 1;
end

end


% alphaIC es una función auxiliar, para NewtonBLIC, que calcula el tamaño 
% de paso que se debe dar con interpolación cuadrática.
function[alpha] = alphaIC(func, x, p)

% Se asigna el valor de 1.0e-4 a c para usarlo en la condición de Armijo.
c = 1.0e-4;
phi0 = feval(func, x);
dphi0 = numgrad(func, x)' * p;
a1 = 1;
phi1 = feval(func, x + a1 * p);

% Mientras la condición de Armijo no se cumpla se seguirá haciendo
% interpolación cuadrática con phi(0), phi'(0) y phi(alpha_k).
while( phi1 > phi0 + c * a1 * dphi0 )
	a1 = - (dphi0 * a1^2) / (2 * (phi1 - phi0 - a1 * dphi0));
	phi1 = feval(func, x + a1 * p);
end
alpha = a1;

end 
