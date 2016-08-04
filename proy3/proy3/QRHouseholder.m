function [ Q,R ] = QRHouseholder( A )

% Programa que calcula la factorizacion QR de la Matriz A
% por medio de transformaciones de Househoder.

[m,n] = size (A);
R = A;
Q = eye(m);

if( m == n )
    lim = n-1;
else
    lim = min(n,m);
end

for k = 1:lim
    a = R(k:m,k);
    e = zeros(m-k+1,1);
    e(1) = 1;
    v = a + sign(a(1)) * norm(a,2) * e;
    v = v / norm(v,2);
    Qk = eye(m-k+1) - (2 * (v * v'));
    Qk = [eye(k-1) zeros(k-1,m-k+1) ;zeros(m-k+1,k-1) Qk];
    R = Qk * R;
    Q = Qk * Q;
end
Q = Q';
end

