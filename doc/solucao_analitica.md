# Solução Analítica dos Modos TE e TM em Guias de Onda Retangulares

A solução analítica dos modos em um guia de onda retangular parte da equação de Helmholtz, resultante da separação das equações de Maxwell no vácuo.

Neste projeto, consideramos guias retangulares metálicos com proporção $a = 2b$ (largura $a$ e altura $b$).

---

## Modos TM (Transversais Magnéticos)

Para os modos TM, a componente longitudinal do campo elétrico $E_z$ é diferente de zero, e o campo magnético não possui componente na direção axial ($H_z = 0$).

A solução geral para $E_z(x, y)$ é:

$$
E_z(x, y) = \sin\left( \frac{m\pi x}{a} \right) \sin\left( \frac{n\pi y}{b} \right)
$$

com $m, n = 1, 2, 3, \dots$ (ambos diferentes de zero para TM).

O número de onda transversal é dado por:

$$
k_c = \pi \sqrt{ \left( \frac{m}{a} \right)^2 + \left( \frac{n}{b} \right)^2 }
$$

As componentes transversais dos campos são:

$$
E_x = -\frac{j\beta}{k_c^2} E_0 \frac{m\pi}{a} \cos\left( \frac{m\pi x}{a} \right) \sin\left( \frac{n\pi y}{b} \right)
$$

$$
E_y = -\frac{j\beta}{k_c^2} E_0 \frac{n\pi}{b} \sin\left( \frac{m\pi x}{a} \right) \cos\left( \frac{n\pi y}{b} \right)
$$

$$
H_x = \frac{j\omega \varepsilon}{k_c^2} E_0 \frac{n\pi}{b} \sin\left( \frac{m\pi x}{a} \right) \cos\left( \frac{n\pi y}{b} \right)
$$

$$
H_y = -\frac{j\omega \varepsilon}{k_c^2} E_0 \frac{m\pi}{a} \cos\left( \frac{m\pi x}{a} \right) \sin\left( \frac{n\pi y}{b} \right)
$$

---

## Modos TE (Transversais Elétricos)

Nos modos TE, o campo elétrico longitudinal é nulo ($E_z = 0$) e o campo magnético $H_z$ possui a forma:

$$
H_z(x, y) = \cos\left( \frac{m\pi x}{a} \right) \cos\left( \frac{n\pi y}{b} \right)
$$

com $m, n = 0, 1, 2, \dots$ e pelo menos um deles diferente de zero.

O número de onda transversal também é:

$$
k_c = \pi \sqrt{ \left( \frac{m}{a} \right)^2 + \left( \frac{n}{b} \right)^2 }
$$

As componentes transversais dos campos são:

$$
H_x = \frac{j\beta}{k_c^2} H_0 \frac{m\pi}{a} \sin\left( \frac{m\pi x}{a} \right) \cos\left( \frac{n\pi y}{b} \right)
$$

$$
H_y = \frac{j\beta}{k_c^2} H_0 \frac{n\pi}{b} \cos\left( \frac{m\pi x}{a} \right) \sin\left( \frac{n\pi y}{b} \right)
$$

$$
E_x = -\frac{j\omega \mu}{k_c^2} H_0 \frac{n\pi}{b} \cos\left( \frac{m\pi x}{a} \right) \sin\left( \frac{n\pi y}{b} \right)
$$

$$
E_y = \frac{j\omega \mu}{k_c^2} H_0 \frac{m\pi}{a} \sin\left( \frac{m\pi x}{a} \right) \cos\left( \frac{n\pi y}{b} \right)
$$

---

## Frequência de Corte

A frequência de corte $f_c$ está relacionada a $k_c$ por:

$$
f_c = \frac{c}{2\pi} k_c = \frac{c}{2} \sqrt{ \left( \frac{m}{a} \right)^2 + \left( \frac{n}{b} \right)^2 }
$$

onde:
- $m, n$ são inteiros que caracterizam o modo,
- $a, b$ são as dimensões do guia (com $a = 2b$ neste caso),
- $c$ é a velocidade da luz no vácuo.

Essas expressões fornecem os valores teóricos que podem ser comparados com os obtidos numericamente pelo método dos elementos finitos (FEM).

