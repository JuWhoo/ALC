%Molécula de Metano - Distâncias Interatômicas em Angstrom:
d52=d53=d54=d34=d32=d24= 1.778 % entre Hs
d12=d13=d14=d51=1.089 % entre C e Hs
%Coordenadas do átomo 1 de C:
u1=v1=w1=0
%Cordenadas átomo 2 de H:
u2=v2=0
w2=1.089
%Coordenadas átomo 3 de H:
v3=0
w3=((d13**2)-(d32**2)+(w2**2))/(2*w2)
u3=sqrt((d13**2)-(w3**2))
%Coordenadas átomo 4 de H:
w4=((d14**2)-(d24**2)+(w2**2))/(2*w2)
u4=((d14**2)-(d34**2)-(2*w3*w4)+(u3**2)+(w3**2))/(2*u3)
v4= sqrt((d14**2)-(u4**2)-(w4**2))
%Vetores de cooredenadas para cada átomo:
x1=[u1 v1 w1]'
x2=[u2 v2 w2]'
x3=[u3 v3 w3]'
x4=[u4 v4 w4]'
%Calculando a Matriz A:
A=2*[u1-u2 v1-v2 w1-w2;
u1-u3 v1-v3 w1-w3;
u1-u4 v1-v4 w1-w4]

b=[((norm(x1)**2)-(norm(x2)**2)) - ((d51**2)-(d52**2));
((norm(x1)**2)-(norm(x3)**2)) - ((d51**2)-(d53**2));
((norm(x1)**2)-(norm(x4)**2)) - ((d51**2)-(d54**2))]

function [L, U, P, x] = luFactorizationPivot(A, b)
    % Verifica se a matriz A é quadrada
    [m, n] = size(A);
    if m ~= n
        error('A matriz A deve ser quadrada');
    end
    % Verifica se o vetor b tem o mesmo número de linhas que A
    if length(b) ~= n
        error('O vetor b deve ter o mesmo número de linhas que a matriz A');
    end
    % Inicializa matrizes L, U e P
    L = eye(n); % Inicializa L como a matriz identidade
    U = zeros(n); % Inicializa U como uma matriz nula
    P = eye(n); % Inicializa P como a matriz identidade
        % Fatoração LU com pivotamento parcial
    for k = 1:n-1
        % Encontra o pivô (elemento de maior magnitude) na coluna k
        [~, pivot_row] = max(abs(A(k:n, k)));
        pivot_row = pivot_row + k - 1; % Ajusta o índice da linha
        
        % Troca as linhas k e pivot_row em A, L e P
        A([k,pivot_row],:) = A([pivot_row,k],:);
        L([k,pivot_row],1:k-1) = L([pivot_row,k],1:k-1);
        P([k,pivot_row],:) = P([pivot_row,k],:);
        
        % Atualiza a parte superior de U
        for j = k:n
            U(k,j) = A(k,j) - L(k,1:k-1)*U(1:k-1,j);
        end
        
        % Calcula os elementos da coluna k de L
        for i = k+1:n
            L(i,k) = A(i,k) - L(i,1:k-1)*U(1:k-1,k);
            L(i,k) = L(i,k) / U(k,k);
        end
    end
    
    % Resolve Ly = P*b
    y = zeros(n, 1);
    Pb = P*b;
    for i = 1:n
        y(i) = Pb(i) - L(i,1:i-1)*y(1:i-1);
    end
    
    % Resolve Ux = y
    x = zeros(n, 1);
    for i = n:-1:1
        x(i) = (y(i) - U(i,i+1:n)*x(i+1:n)) / U(i,i);
    end
end
disp('Matriz L:');
disp(L);
disp('Matriz U:');
disp(U);
disp('Matriz P:');
disp(P);
disp('Solução do sistema:');
disp(x);
