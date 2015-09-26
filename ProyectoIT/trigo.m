% trigo.m modela la producción de trigo por hectárea de un cultivo
% con la función logística.
% Para realizar esto se utiliza la información de la producción
% durante 24 años; esta información nos permite construir
% residuales que se usan para estimar ciertos parámetros del modelo
% logístico.
close all;

% El vector t representa a los primeros 24 años mientras que
% el vector prodt tiene la información de la producción de trigo
% durante los respectivos años.
t = [1:24]';
prodt = [11.72, 13.38, 14.10, 13.87, 14.80, 15.68, 14.36, 16.30,...
	 16.91, 18.16, 18.43, 18.70, 20.46, 19.16, 20.02, 22.41,...
	 21.21, 22.81, 23.97, 23.27, 23.80, 25.59, 24.93, 26.59]'; 

% La función logística con
% 		r = x(1)
% 		K = x(2)
% 		P0 = x(3)
P = @(x, t1) x(2) / (1 + (x(2) / x(3) - 1) * exp(-x(1) * t1));

% Con base en los datos obtenidos y la función se crea la
% función obejetivo que es un medio de la suma de los residuales.
fobj = @(x) 0;
for k = 1:24
	fobj = @(x) (P(x, t(k)) - prodt(k))^2 / 2 + fobj(x);
end

% Se establece el punto inicial para el método de búsqueda lineal
% con paso de Newton.
x0 = [0.05 1 30]';

% Se resuelve el problema con el paso de Newton y búsqueda lineal
% con interpolación cuadrática y se imprime el resultado.
[x_sol, iter] = NewtonBLIC(fobj, x0, 1.0e-6);
fprintf('Los parámetros son:\n');
fprintf('r  = %1.4e\n', x_sol(1));
fprintf('K  = %1.4e\n', x_sol(2));
fprintf('P0 = %1.4e\n', x_sol(3));

% Se le asignan los parámetros encontrados a la función logística
% y se calculan los valores aproximados en cada trimestre y se
% guardan en la variable Pt.
P_sol = @(t1) P(x_sol, t1);
Pt = zeros(16, 1);

for k = 1:24
	Pt(k) = P_sol(t(k));
end

% Se grafica la solución final junto con los datos que se usaron.
plot(t, Pt, 'or', 'linewidth', 3)
hold on;
plot(t, prodt, 'db', 'linewidth', 3)
xlabel('Año')
ylabel('Producción por hectárea')
legend('Producción estimada', 'Producción observada')
hold off;
