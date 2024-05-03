%Molécula de Metanal - Distâncias Interatômicas em Angstrom:
d42=d32= 2.006 % entre H e C
d34=d43= 1.884 % entre H e O
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

%Verificando se o determinante é não nulo:
det(A)

%Aplicando Eliminação de Gauss
E = [A b] % matriz ampliada
n=size(E,1)

%para baixo
for i=1:n-1
  if (abs(E(i,i)) <= 0)
    for j=i+1:n
      if (abs(E(j,i)) > 0)
        break
      end
    endfor
    aux=E(i,:);
    E(i,:)=E(j,:);
    E(j,:)=aux;
  end
  E(i+1:n,:) -= E(i+1:n,i)/E(i,i)*E(i,:);
endfor
E

%para cima
for i=n:-1:2
  E(i,:) = E(i,:)/E(i,i);
  E(1:i-1,:) -= E(1:i-1,i)*E(i,:);
endfor
E(1,:) = E(1,:)/E(1,1)