%Molécula de Metanal - Distâncias Interatômicas em Angstrom:
d42=d32= 2.006 % entre H e O
d34=d43= 1.884 % entre H e H
d12=1.229 % entre O e C
d13=d41=1.089 % entre H e C
%Coordenadas calculadas:
u1=0
v1=0
u2=1.229
v2=0
u3=((d13**2)-(d32**2)+(u2**2))/(2*u2)
v3=sqrt((d13**2)-(u3**2))
%Vetores de cooredenadas para cada átomo:
x1=[u1 v1]'
x2=[u2 v2]'
x3=[u3 v3]'
%Calculando a Matriz A:
A=2*[u1-u2 v1-v2;
u1-u3 v1-v3]

b=[((norm(x1)**2)-(norm(x2)**2)) - ((d41**2)-(d42**2));
((norm(x1)**2)-(norm(x3)**2)) - ((d41**2)-(d43**2))]

function [L, U, x] = luFactorization(A, b)
    % Verifica se a matriz A é quadrada
    [m, n] = size(A);
    if m ~= n
        error('A matriz A deve ser quadrada');
    end
    
    % Verifica se o vetor b tem o mesmo número de linhas que A
    if length(b) ~= n
        error('O vetor b deve ter o mesmo número de linhas que a matriz A');
    end
    
    % Inicializa matrizes L e U
    L = eye(n); % Inicializa L como a matriz identidade
    U = zeros(n); % Inicializa U como uma matriz nula
    
    % Fatoração LU
    for k = 1:n
        % Atualiza a parte superior de U
        for j = k:n
            U(k,j) = A(k,j) - L(k,1:k-1)*U(1:k-1,j);
        end
        
        % Calcula os elementos da coluna k de L
        for i = k+1:n
            L(i,k) = (A(i,k) - L(i,1:k-1)*U(1:k-1,k)) / U(k,k);
        end
    end
    
    % Resolve Ly = b
    y = zeros(n, 1);
    for i = 1:n
        y(i) = (b(i) - L(i,1:i-1)*y(1:i-1)) / L(i,i);
    end
    
    % Resolve Ux = y
    x = zeros(n, 1);
    for i = n:-1:1
        x(i) = (y(i) - U(i,i+1:n)*x(i+1:n)) / U(i,i);
    end
end
[L, U, x] = luFactorization(A, b);
disp('Matriz L:');
disp(L);
disp('Matriz U:');
disp(U);
disp('Solução do sistema:');
disp(x);
