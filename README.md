  Nesse repositório estão os códigos utilizados no artigo feito como avaliação final para disciplina de Álgebra Linear Computacional, intitulado "Aplicação de métodos diretos para solução do Problema de Geometria de Distâncias Moleculares (PGDM)".

  Os códigos foram construídos em Octave e estão separados em quatro arquivos, dois métodos diretos aplicados à duas moléculas diferentes, buscando a resolução do cálculo do PGDM para estas.

  No início de todos os arquivos é realizado primeiro o cálculo das matrizes A e b para o sistema da molécula correspondente. 
  
  Os cálculos das matrizes, reduzidos a um sistema linear simples, é baseado no método utilizado por DONG & WU no artigo "A linear-time algorithm for solving the molecular distance geometry problem with exact inter-atomic distances".

  Nos arquivos 'calculos_metanal_elim_gauss_sem_pivot.m' e 'calculos_metano_elimin_gauss_pivot.m' é aplicado o método de Eliminação de Gauss nas duas moléculas estudadas, o Metanal (H2CO) e ao Metano (CH4).
  
  O código utilizado para Eliminação de Gauss foi construído e disponibilizado pelo professor Pedro Konzen, em seu Github, e que podem ser acessado em https://github.com/phkonzen/notas. Nele só é aplicado pivotamento caso seja realmente necessário.
  
  Nos arquivos 'metanal_fatoracao_lu_sem_pivot.m' e 'metano_fatoracao_lu_pivot.m' é aplicado a resolução dos sistemas por Fatoração LU sem pivotamento e com pivotaamento, respectivamente, de acordo com a necessidade do sistema.

  Ambos os métodos aplicados chegaram a soluções iguais e concluiu-se que neste caso é possível solucionar o PGDM com aplicação destes métodos diretos.


# REFERÊNCIA:
DONG, Q., WU, Z. no artigo "A linear-time algorithm for solving the molecular distance geometry problem with exact inter-atomic distances". Journal of Global Optimization}, 22, 365-375, 2002.
