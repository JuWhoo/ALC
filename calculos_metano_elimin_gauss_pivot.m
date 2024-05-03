%Molécula de Metano - Distâncias Interatômicas em Angstrom:
d52=d53=d54=d34=d32=d24= 1.778 % entre Hidrogenios
d12=d13=d14=d51=1.089 % entre C e H
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

%Verificando se o determinante é não nulo:
det(A)

%Aplicando a eliminação de Gauss:
E = [A b] 

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
