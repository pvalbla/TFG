# Modelos de Regresión para Datos de Conteos Multivariantes

Este repositorio contiene la aplicación práctica y el desarrollo computacional de mi **Trabajo de Fin de Grado (TFG) en Estadística**.

## Estructura del Proyecto

El análisis se desglosa en tres scripts de `R` independientes, correspondientes a los capítulos metodológicos del trabajo:

* **`01_modelo_multinomial.R` (Capítulo 1):** Ajuste del modelo Multinomial clásico usando el dataset `Heating` del paquete `mlogit`, que contiene las decisiones de $n=900$ familias de California sobre su sistema de calefacción centralizado según sus ingresos, edad y tamaño de la vivienda.
* **`02_modelo_dirichlet-multinomial.R` (Capítulo 2):** Ajuste del modelo Dirichlet-Multinomial estándar usando el dataset `us` del paquete `dirmult`, que recopila las frecuencias de los alelos observados en una muestra de población estadounidense para analizar la sobredispersión genética.
* **`03_modelo_dirichlet-multinomial_generalizado.R` (Capítulo 3):** Ajuste del modelo Dirichlet-Multinomial Generalizado usando el dataset `rnaseq` del paquete `MGLM`, que contiene conteos de lecturas de expresión genética (ARN-Seq) distribuidos en 6 exones para evaluar el impacto de covariables clínicas.

## Requisitos e Instalación

Para replicar los análisis y ejecutar los scripts locales, asegúrate de tener instalado **R** (versión $\ge 4.0$) junto con las siguientes librerías especializadas:

```R
install.packages("mlogit")
install.packages("dirmult")
install.packages("MGLM")
