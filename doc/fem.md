# Abordagem com o Método dos Elementos Finitos (FEM)

O Método dos Elementos Finitos (FEM) é uma técnica numérica amplamente utilizada para resolver equações diferenciais parciais em domínios complexos. Neste projeto, utilizamos o FEM para resolver o problema de autovalores associado à propagação de modos eletromagnéticos TE e TM em guias de onda **retangulares metálicos** com proporção $` a = 2b `$. A malha é composta por elementos triangulares.

---

## Formulação Variacional

A equação de Helmholtz governando o campo escalar ($` E_z `$ ou $` H_z `$) é:

$`
\nabla^2 \psi + k_c^2 \psi = 0
`$

onde $` \psi `$ representa a componente longitudinal $` E_z `$ (modo TM) ou $` H_z `$ (modo TE).

A forma fraca desta equação é obtida multiplicando por uma função teste $` v `$ e integrando sobre o domínio:

$`
\int_\Omega \nabla v \cdot \nabla \psi \, d\Omega = k_c^2 \int_\Omega v \psi \, d\Omega
`$

---

## Discretização

- Utilizamos elementos triangulares lineares (P1).
- A malha do domínio retangular é gerada automaticamente com o **GMSH**.
- Cada função de forma é associada a um nó da malha.

---

## Montagem das Matrizes

A forma fraca leva à montagem de duas matrizes:

- Matriz de rigidez $` S `$: representa o termo $` \nabla v \cdot \nabla \psi `$
- Matriz de massa $` T `$: representa o termo $` v \psi `$

O problema de autovalores torna-se:

$`
S x = \lambda T x
`$

com $` \lambda = k_c^2 `$, sendo $` k_c `$ o número de onda transversal.

---

## Condições de Contorno

- **Modo TM**: aplica-se condição de contorno de Dirichlet em $` E_z `$ (campo nulo na parede metálica).
- **Modo TE**: condição natural de Neumann (fluxo normal do gradiente do campo é nulo).

---

## Resolução Numérica

O problema generalizado é resolvido com os métodos do `scipy.sparse.linalg.eigsh` para obter os autovalores $` \lambda `$, de onde se calcula:

- $` k_c = \sqrt{\lambda} `$
- $` f_c = \frac{c}{2\pi} k_c `$
- $` kc \cdot a `$: valor adimensional usado neste projeto
