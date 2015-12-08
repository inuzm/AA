function [Jx] = jacobiana(fname,x)
% Esta funcion aproxima por diferencia finitas la
% matriz jacobiana de fname en x, donde
% fnama:R^n --> R^n.

% In
% fname.- cadena con el nombre de la funcion.
% x    .- vector de dimension n.
% Out
% Jx.- matriz cuadrada de orden n.


h = 1.e-06;

Fx = feval(fname,x);
n = length(x);
m = length(Fx);

Jx = zeros(m,n);

for k = 1:n
    xh = x; 
    xh(k) = xh(k) + h;
    Fxh = feval(fname,xh);
    Jx(:,k) = (Fxh - Fx)/ h;
end
