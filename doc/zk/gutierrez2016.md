---
title: "Efficient and automated large-scale detection of structural relationships in proteins with a flexible aligner"
name: Gutierrez_2016
doi: 10.1186/s12859-015-0866-8
year: 2016
author: Fernando I. Gutiérrez and Felipe Rodriguez-Valenzuela and Ignacio L. Ibarra and Damien P. Devos and Francisco Melo
tags: refract structural_alignment protein moma
---

# Efficient and automated large-scale detection of structural relationships in proteins with a flexible aligner

MOMA es una herramiento de alineamiento estructural que utiliza pares de estructuras secundarias (SS, en corto) como información para llevar a cabo el alineamiento.

Utilizando el ángulo ($\omega$) y distancia ($d$) para cada par de SS podemos describir la dispocision topologica de la proteína, permitiendo asi su alineamiento.
Estas caracteristicas son codificadas en una matriz de la forma:

$$
\begin{equation}
  \begin{bmatrix}
  SS_1SS_2     & d_{1,2}      & d_{1,3}      & d_{1,4} \\
  \omega_{1,2} & SS_2SS_3     & d_{2,3}      & d_{2,4} \\
  \omega_{1,3} & \omega_{2,3} & SS_2SS_3     & d_{3,4} \\
  \omega_{1,4} & \omega_{2,4} & \omega_{3,4} & SS_3SS_4 \\
  \end{bmatrix}
\end{equation}
$$

donde $\omega_{i,j}$ es en angulo entro un par de SS, $d_{i,j}$ es la distancia entre un par de SS y $SS_{n}$ corresponde a una estructura secundaria en la posicion $n$.
Esta ultima puede ser de 2 tipos: $A$ si la SS corresponde a una $\alpha$-hélice o $B$ si la SS corresponde a una hebra-$\beta$.

A partir de esta matriz descriptiva podemos generar el alineamiento estructural via un comparasion de a pares para cada par de SS entre las 2 proteinas.
El procedimiento para esto es el siguiente:
1. Comparación estructural entre las 2 proteinas a partir de las matrices descriptivas por medio de una comparasion a pares de cada par de SS presentes.
    - Aqui hay una serie de parametros libres que se pueden optimizar y seleccionar dependiendo del tipo de busqueda que se esta haciendo.
2. Extraccion de "bloques" de SS que son significativamente similares entre ambas proteinas.
3. Calculo de matrices de transformacion para alinear de manera flexible cada segmento de las proteinas.

MOMA es equivalente en rendimiento predictivo a otros metodos del estado del arte, pero es entre 2 y 3 ordenes de magnitud más rápido.
