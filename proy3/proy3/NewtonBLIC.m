% NewtonBLIC.m resuelve el problema de la forma:
% 		min f(x)
% sin restricciones. Esto lo hace con el paso de Newton e interpolaci?n
% cuadr?tica.
% Argumentos de entrada:
% func: la funci?n f que se quiere minimizar.
% x: punto inicial para el m?todo.
% tol: la tolerancia para la condici?n de paro; este argumento es opcional.
% maxiter: el m?ximo de iteraciones que se har?n en el m?todo; este
% argumento es opcional
% Argumentos de salida:
% x: punto en el que la norma del gradiente es menor que la tolerancia.
% k: n?mero de iteraciones en el que se lleg? al valor de x.
function[x, k] = NewtonBLIC(func, x, tol, maxiter)

% Si no se dan todos los argumentos entonces se revisa cu?les no se dieron.
% Si no se dio maxiter se le asignar? el valor de 10000. Adem?s, si no se
% da tol se le asignar? el valor de 1.0e-8.
if( nargin < 4 )
	maxiter = 10000;
	if( nargin < 3 )
		tol = 1.0e-8;
	end
end

% Se calculan el gradiente y la matriz hessiana de la funci?n en el punto
% inicial y se le da el valor de 0 a k.
gf = numgrad(func, x);
Hf = numhess(func, x);
k = 0;

% Mientras no se cumpla la condici?n de paro, que la norma del gradiente
% sea menor que cierta tolerancia o que se llegue al m?ximo de iteraciones,
% se realiza la minimizaci?n de f con el m?todo de Newton.
while( norm(gf) > tol && k < maxiter )
    p = - Hf \gf;
	alpha = alphaIC(func, x, p);
    x = x + alpha * p;
	gf = numgrad(func, x);
	Hf = numhess(func, x);
	k = k + 1;
end

end


% alphaIC es una funci?n auxiliar, para NewtonBLIC, que calcula el tama?o 
% de paso que se debe dar con interpolaci?n cuadr?tica.
function[alpha] = alphaIC(func, x, p)

% Se asigna el valor de 1.0e-4 a c para usarlo en la condici?n de Armijo.
c = 1.0e-4;
phi0 = feval(func, x);
dphi0 = numgrad(func, x)' * p;
a1 = 1;
phi1 = feval(func, x + a1 * p);

% Mientras la condici?n de Armijo no se cumpla se seguir? haciendo
% interpolaci?n cuadr?tica con phi(0), phi'(0) y phi(alpha_k).
while( phi1 > phi0 + c * a1 * dphi0 )
	a1 = - (dphi0 * a1^2) / (2 * (phi1 - phi0 - a1 * dphi0));
	phi1 = feval(func, x + a1 * p);
end
alpha = a1;

end 