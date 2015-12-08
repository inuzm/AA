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
% función obejetivo que es un vector de residuales.
robj = @(x) [];
for k = 1:24
	robj = @(x) [robj(x); P(x, t(k))-prodt(k)];
end

% Se establecen los puntos iniciales para los métodos.
x0 = [0.005 1 30]'; % Gauss Newton
x1 = [0.005 13 4]'; % L-M

% Se resuelve el problema con ambos métodos:
[x_sol, iter] = MetodoGaussNewton(robj, x0);
[x_sol1, iter1] = LevenbergMarquadt(robj, x1);

fprintf('Los parámetros para GN son:\n');
fprintf('r  = %1.4e\n', x_sol(1));
fprintf('K  = %1.4e\n', x_sol(2));
fprintf('P0 = %1.4e\n', x_sol(3));

fprintf('Los parámetros para LM son:\n');
fprintf('r  = %1.4e\n', x_sol1(1));
fprintf('K  = %1.4e\n', x_sol1(2));
fprintf('P0 = %1.4e\n', x_sol1(3));

% Se le asignan los parámetros encontrados a la función logística
% y se calculan los valores aproximados en cada trimestre y se
% guardan en la variable Pt.
P_sol = @(t1) P(x_sol, t1);
Pt = zeros(16, 1);

for k = 1:24
	Pt(k) = P_sol(t(k));
end

P_sol2 = @(t1) P(x_sol1, t1);
Pt2 = zeros(16, 1);

for k = 1:24
	Pt2(k) = P_sol2(t(k));
end

% Se grafican ambas soluciones finales junto con los datos que se usaron.
subplot(2, 2, 1)
plot(t, Pt, 'or', t, prodt, 'db', 'linewidth', 3)
xlabel('Año')
ylabel('Producción por hectárea')
legend('Producción estimada GN', 'Producción observada')

subplot(2, 2, 2)
plot(t, Pt2, 'ok', t, prodt, 'db', 'linewidth', 3)
xlabel('Año')
ylabel('Producción por hectárea')
legend('Producción estimada LM', 'Producción observada')

subplot(2, 2, 3)
plot(t, Pt, 'or', t, Pt2, 'ok', t, prodt, 'db', 'linewidth', 3)
xlabel('Año')
ylabel('Producción por hectárea')
legend('Producción estimada GN', 'Producción estimada LM',...
'Producción observada')

hold off;
