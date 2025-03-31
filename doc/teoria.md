# Fundamentos Teóricos dos Modos TE e TM em Guias de Onda Retangulares

Guias de onda retangulares são estruturas metálicas ocas, normalmente utilizadas para conduzir ondas eletromagnéticas em frequências de micro-ondas. A solução das equações de Maxwell em um guia retangular com paredes perfeitamente condutoras leva à existência de dois tipos principais de modos propagantes:

- **Modos Transversais Elétricos (TE)**: o campo elétrico não possui componente ao longo do eixo do guia ($` E_z = 0 `$);
- **Modos Transversais Magnéticos (TM)**: o campo magnético não possui componente axial ($` H_z = 0 `$).

A análise é feita assumindo propagação na direção longitudinal ($` z `$), com campos que variam no plano transversal ($` x, y `$). A separação de variáveis na equação de Helmholtz leva a soluções do tipo seno e cosseno, dependentes dos inteiros $` m `$ e $` n `$.

---

## Equação de Helmholtz

Para ambos os casos (TE e TM), a componente longitudinal ($` \psi = E_z `$ ou $` H_z `$) obedece à equação:

$`
\nabla_t^2 \psi + k_c^2 \psi = 0
`$

onde $` \nabla_t^2 `$ é o laplaciano transversal e $` k_c `$ é o número de onda transversal.

---

## Modos TE

Nos modos TE ($` E_z = 0 `$), a solução para $` H_z(x, y) `$ é:

$`
H_z(x, y) = \cos\left( \frac{m\pi x}{a} \right) \cos\left( \frac{n\pi y}{b} \right)
`$

com $` m, n \ge 0 `$ e pelo menos um deles diferente de zero.

---

## Modos TM

Nos modos TM ($` H_z = 0 `$), a solução para $` E_z(x, y) `$ é:

$`
E_z(x, y) = \sin\left( \frac{m\pi x}{a} \right) \sin\left( \frac{n\pi y}{b} \right)
`$

com $` m, n \ge 1 `$. Ambos devem ser diferentes de zero para garantir $` E_z \neq 0 `$ no domínio.

---

## Frequência de Corte

O número de onda transversal e a frequência de corte associada são dados por:

$`
k_c = \pi \sqrt{ \left( \frac{m}{a} \right)^2 + \left( \frac{n}{b} \right)^2 }
`$

$`
f_c = \frac{c}{2} \sqrt{ \left( \frac{m}{a} \right)^2 + \left( \frac{n}{b} \right)^2 }
`$

onde:
- $` c `$ é a velocidade da luz no vácuo,
- $` a, b `$ são as dimensões do guia,
- $` m, n `$ são os índices dos modos.

No caso deste projeto, consideramos guias com proporção $` a = 2b `$, o que influencia diretamente os modos que têm menor frequência de corte (como o modo $` TE_{10} `$).
