% pruebaminij.m
%Este script file sirve para resolver sistemas de ecuaciones de la forma:
% A * x = b
%Donde A es una matriz de minij
%b es un vector de unos

%Se encuentra la respuesta x por dos métodos diferentes, Gradiente conjugado y 
%Gradiente conjugado Precondicionado.
%El precondicionamiento se da de dos formas diferentes, una es usando una matriz
%de Jacobi y el otro es usando Cholesky incompleto que proporciona matlab.

%La salida es la impresión de el tipo de precondicionamiento que se uso, el tamaño de 
%la matriz de minij, el número de iteraciones que tomo GCPre con ese precondicionamiento
%y la norma del punto encontrado menos el encontrado por Gradiente Conjugado sin precondicionar
S = char('Jacobi', 'Cholesky incompleto');
n = 10:5:20;

fprintf('\nPrueba Minij\n');
display('	     Matriz 	 n 	 iter 	 ||x_{GCPre}-X_{GC}||')
display('-------------------------------------------------------------')
for k = 1:3
	A = gallery('minij', n(k));
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
