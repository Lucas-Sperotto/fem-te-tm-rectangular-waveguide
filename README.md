# FEM TE/TM Rectangular Waveguide (Python)

Este repositório contém uma implementação em **Python** para simulação numérica dos **modos eletromagnéticos TE e TM** em **guias de onda retangulares metálicos**, com proporção entre lados $` a = 2b `$, utilizando o **método dos elementos finitos (FEM)**.

## 🔧 Recursos

- Geração automática da malha retangular via **GMSH**;
- Cálculo dos autovalores e frequências de corte para os modos TE e TM;
- Visualização dos modos escalares (Ez ou Hz);
- Visualização vetorial do campo transversal com escala e cor ajustadas;
- Exportação dos dados em `.csv` para análise;
- Compatível com Python 3.8+.

## 📦 Dependências

As bibliotecas necessárias estão listadas em `requirements.txt`. Instale com:

```bash
pip install -r requirements.txt
```

## ▶️ Como executar

Certifique-se de ter o **GMSH** instalado e acessível no PATH.  
Depois, execute o script principal:

```bash
python rectangular_waveguide_modes.py
```

## 📁 Estrutura de saída

- `out/img/`: imagens dos modos escalares e vetoriais.
- `out/results/`: tabelas com autovalores, frequências e kc·a.
- `meshes/`: arquivos de malha `.msh` gerados pelo GMSH.

## 📑 Documentação Técnica

A seguir estão os arquivos com explicações teóricas e metodológicas:

- [Fundamentos Teóricos dos Modos TE/TM](/doc/teoria.md)
- [Solução Analítica para Guia Retangular](/doc/solucao_analitica.md)
- [Abordagem pelo Método dos Elementos Finitos](/doc/fem.md)

## 📊 Resultados Numéricos e Imagens

- [Modos TE – Potencial Hz](/doc/resultados_te_potencial.md)
- [Modos TE – Campo Transversal](/doc/resultados_te_transversal.md)
- [Modos TM – Potencial Ez](/doc/resultados_tm_potencial.md)
- [Modos TM – Campo Transversal](/doc/resultados_tm_transversal.md)

## 📜 Licença

Este projeto está licenciado sob a **Licença MIT**.
