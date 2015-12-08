% ipad.m modela las ventas del iPad con la funci?n log?stica.
% Para realizar esto utiliza la informaci?n de las ventas durante
% 16 trimestres; esta informaci?n nos permite construir residuales
% que se usan para estimar ciertos par?metros del modelo log?stico.
close all;

% El vector t representa a los primeros 16 trimestres mientras que
% el vector v tiene la informaci?n de la venta del iPad durante
% los respectivos trimestres.
t = [1:16]';
v = [3.27 4.19 7.33 4.69 9.25 11.12 15.30 11.80 17.00 14.00 22.90...
	 19.50 14.60 14.10 26.00 16.35]';

% La funci?n log?stica con
% 		r = x(1)
% 		K = x(2)
% 		P0 = x(3)
P = @(x, t1) x(2) / (1 + (x(2) / x(3) - 1) * exp(-x(1) * t1));

% Con base en los datos obtenidos y la funci?n se crea la
% funci?n obejetivo que es un medio de la suma de los residuales.
robj = @(x)[];
for k = 1:16
	robj = @(x) [robj(x); P(x, t(k))-v(k)];
end

% Se establece el punto inicial para el m?todo de b?squeda lineal
% con paso de Newton.
x0 = [0.1 30 3.5]';

% Se resuelve el problema con el paso de Newton y b?squeda lineal
% con interpolaci?n cuadr?tica y se imprime el resultado.
% [x_sol, iter] = NewtonBLIC(fobj, x0)
%[x_sol, iter] = GaussNewton(robj, x0);
[x_sol, iter] = GaussNewtonRC(robj, x0);
fprintf('Los par?metros son:\n');
fprintf('r  = %1.4e\n', x_sol(1));
fprintf('K  = %1.4e\n', x_sol(2));
fprintf('P0 = %1.4e\n', x_sol(3));

% Se le asignan los par?metros encontrados a la funci?n log?stica
% y se calculan los valores aproximados en cada trimestre y se
% guardan en la variable Pt.
P_sol = @(t1) P(x_sol, t1);

Pt = zeros(16, 1);
for k = 1:16
	Pt(k) = P_sol(t(k));
end

% Se grafica la soluci?n final junto con los datos que se usaron.
plot(t, Pt, 'or', 'linewidth', 3)
hold on;
plot(t, v, 'db', 'linewidth', 3)
xlabel('Trimestre')
ylabel('Ventas (en millones)')
legend('Ventas estimadas', 'Ventas observadas')
hold off;