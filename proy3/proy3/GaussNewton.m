function [x,iter] = GaussNewton(fname,x0)

maxiter = 1000;
maxjter = 10;
tol = 1.e-08;
x = x0;
n = length(x);
size(x)
rx = feval(fname,x);
size(rx)
Jx = jacobiana(fname,x0);
size(Jx)
gx = Jx'*rx;
iter = 0;
disp( 'iter  ||rgx||    ' )
disp( '---------------------' )

% Iteraciones
while(iter < maxiter && norm(gx) > tol)
    
    rango = rank(Jx);
    if rango < n
        B = Jx'*Jx + norm(gx)*eye(n);
        p = -B\gx;
    else
%       Minimos cuadrados lienales
        [Q,R] = QRHouseholder(Jx);
        bp = Q'*rx;
        p = -R(1:n,1:n)\bp(1:n);
    end
    alfa = 1.0;
    jter = 0;
    xt = x + alfa*p;
    fxt = feval(fname,xt);
    ped = gx'*p;
    
    while ( fxt' * fxt / 2 > rx' * rx / 2+alfa*.1*ped && jter < maxjter)
        alfa = alfa/2;
        xt = x + alfa *p;
        fxt = feval(fname,xt);
        jter = jter + 1;
    end
    
    x = x + alfa*p;
    
%     x = x + p;
    rx = feval(fname,x);
    Jx = jacobiana(fname,x);
    gx = Jx'*rx;
    iter = iter + 1;
    
    disp(sprintf('%2.0f         %2.6f       %2.1e',iter,norm(gx), alfa))
    
end