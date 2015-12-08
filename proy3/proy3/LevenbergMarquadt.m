function [x,iter] = LevenbergMarquadt(fname,x0)

maxiter = 1000;
maxjter = 10;
tol = 1.e-08;
x = x0;
n = length(x);
size(x);
rx = feval(fname,x);
size(rx);
Jx = jacobiana(fname,x0);
size(Jx);
gx = Jx'*rx;
iter = 0;
disp( 'iter  ||rgx||    ' )
disp( '---------------------' )
delta0 = 1;
delta = delta0;

m = @(p) p'*Jx'*Jx*p/2 + p'*gx+rx'*rx/2;

% Iteraciones
while(iter < maxiter && norm(gx) > tol)
   
    B = Jx' * Jx; 
    p = Doblez(B, gx, delta);
    
    ro = (rx'*rx/2 - feval(fname, x + p)' * feval(fname, x + p) / 2)...
        /(m(zeros(n, 1)) - m(p));
    
    if ro < 1/4
        delta = delta/2;
    else
        if ro > 3/4 && norm(p) == delta
            delta = min(2*delta,delta0);
        else
            delta = delta;
        end
    end
    
    if ro > 0.4
        x = x + p;
    else
        x = x;
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
    m = @(p) p'*Jx'*Jx*p/2 + p'*gx+rx'*rx/2;
    iter = iter + 1;
    
    disp(sprintf('%2.0f         %2.6f       %2.1e',iter,norm(gx), alfa))
    
end
