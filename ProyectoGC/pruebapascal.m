% pruebapascal.m
S = char('Jacobi', 'Cholesky incompleto');
n = 10:5:20;

fprintf('\nPrueba Pascal\n');
display('	     Matriz 	 n 	 iter 	 ||x_{GCPre}-X_{GC}||')
display('-------------------------------------------------------------')
for k = 1:3
	A = pascal(n(k));
	b = ones(n(k), 1);
	xgc = GC(A, b);
	C = diag(sqrt(diag(A)));
	[xgcp, iter] = GCPre(A, b, zeros(n(k), 1), C);
	fprintf('%19s 	%3i 	%5i 	%1.15e\n', deblank(S(1,:)), n(k), iter, norm(xgc - xgcp, 2));
	m = floor(n(k) / 3);
	L = chol(A(1:m, 1:m), 'lower');
	D = diag(sqrt(diag(A(m+1:n(k), m+1:n(k)))));
	C = [L zeros(m, n(k)-m); zeros(n(k)-m, m) D];
	[xgcp, iter] = GCPre(A, b, zeros(n(k), 1), C);
	fprintf('%19s 	%3i 	%5i 	%1.15e\n', deblank(S(2, :)), n(k), iter, norm(xgc - xgcp, 2));
end