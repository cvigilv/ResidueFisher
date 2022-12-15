---
title: "Foldseek: fast and accurate protein structure search"
name: van_Kempen_2022
doi: 10.1101/2022.02.07.479398
year: 2022
author: Michel van Kempen and Stephanie Kim and Charlotte Tumescheit and Milot Mirdita and Cameron L.M. Gilchrist and Johannes Soeding and Martin Steinegger
tags: refract structural_alignment protein foldseek
---

# Foldseek: fast and accurate protein structure search

Foldseek corresponde a un algoritmo de alineamiento estructural tridimensional que implementa un alineamiento 1D por medio de un descriptor estructural obtenido de un VQ-VAE.

El descriptor estructural (llamado librería o diccionario 3Di) describe las interacciones ternarias entre residuos de una proteína, las cuales se discretizan por medio de este VQ-VAE que esta entrenado a priori.

  - Esto le permite a Foldseek tener mejor rendimiento gracias a una disminución de la densidad de información, problema clásico de los diccionarios basados en carbonos alfa del backbone.
  - El VQ-VAE se entrena con un descriptor de 10 dimensiones (5 secuencias de aminoácidos, 7 descriptores de ángulos y una medida de distancia entre los 2 carbonos alineados) para discretizar a un alfabeto de 20 estados, trade-off velocidad-precisión.

Foldseek funciona de la siguiente manera:

1. Discretización de la estructura consulta a una secuencia de 3Di
2. Prefiltrado de la secuencia usando un alineamiento k-mer de doble diagonal y alineamiento gapless con MMseq2
  - Implementan una matriz de sustitución tipo BLOSUM para este paso, la cual obtienen del output del VQ-VAE y validaciones cruzadas.
3. Alineamiento local usando 3Di o global usando TM-Align
  - Alineamiento local combina 3Di con sustitución de aminoácidos.

Foldseek presenta resultados que pueden resumirse en 2 puntos principales:
- Rendimiento comparable con el estado el arte.
- Tiempo de computo comparable con métodos basados en alineamiento de secuencia.
  * Eliminando la brecha entre este tipo de métodos! ESTE ES EL PUNTO MAS IMPORTANTE DEL PAPER

El método, básicamente, es:

1. Conversión de 3D a 1D vía VQ-VAE con el diccionario de 3Di
2. Alineamiento de secuencia de 3Di vía MMseq2
3. Superposición estructural con TMalign utilizando el alineamiento obtenido en el paso anterior
