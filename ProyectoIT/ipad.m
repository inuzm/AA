% ipad.m modela las ventas del iPad con la función logística.
% Para realizar esto utiliza la información de las ventas durante
% 16 trimestres; esta información nos permite construir residuales
% que se usan para estimar ciertos parámetros del modelo logístico.
close all;

% El vector t representa a los primeros 16 trimestres mientras que
% el vector v tiene la información de la venta del iPad durante
% los respectivos trimestres.
t = [1:16]';
v = [3.27 4.19 7.33 4.69 9.25 11.12 15.30 11.80 17.00 14.00 22.90...
	 19.50 14.60 14.10 26.00 16.35]';

% La función logística con
% 		r = x(1)
% 		K = x(2)
% 		P0 = x(3)
P = @(x, t1) x(2) / (1 + (x(2) / x(3) - 1) * exp(-x(1) * t1));

% Con base en los datos obtenidos y la función se crea la
% función obejetivo que es un medio de la suma de los residuales.
fobj = @(x) 0;
for k = 1:16
	fobj = @(x) (P(x, t(k)) - v(k))^2 / 2 + fobj(x);
end

% Se establece el punto inicial para el método de búsqueda lineal
% con paso de Newton.
x0 = [0.1 30 3.5]';

% Se resuelve el problema con el paso de Newton y búsqueda lineal
% con interpolación cuadrática y se imprime el resultado.
[x_sol, iter] = NewtonBLIC(fobj, x0);
fprintf('Los parámetros son:\n');
fprintf('r  = %1.4e\n', x_sol(1));
fprintf('K  = %1.4e\n', x_sol(2));
fprintf('P0 = %1.4e\n', x_sol(3));

% Se le asignan los parámetros encontrados a la función logística
% y se calculan los valores aproximados en cada trimestre y se
% guardan en la variable Pt.
P_sol = @(t1) P(x_sol, t1);

Pt = zeros(16, 1);
for k = 1:16
	Pt(k) = P_sol(t(k));
end

% Se grafica la solución final junto con los datos que se usaron.
plot(t, Pt, 'or', 'linewidth', 3)
hold on;
plot(t, v, 'db', 'linewidth', 3)
xlabel('Trimestre')
ylabel('Ventas (en millones)')
legend('Ventas estimadas', 'Ventas observadas')
hold off;
