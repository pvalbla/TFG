# ==============================================================================
# CAPÍTULO 3: MODELO DIRICHLET-MULTINOMIAL GENERALIZADO (GDM)
# Objetivo: Regresión multivariante, contraste de especificación y diagnóstico Wald
# ==============================================================================

# --- 1. CONFIGURACIÓN INICIAL Y CARGA DE DATOS ---
library(MGLM)

# Carga del dataset biológico (Conteo de lecturas de ARN-Seq)
data("rnaseq")

# Definición de la fórmula del modelo de regresión multivariante:
# Bloque izquierdo: Matriz de variables de respuesta (conteos de los 6 exones)
# Bloque derecho: Covariables clínicas y de control del paciente
formula_conteo <- cbind(X1, X2, X3, X4, X5, X6) ~ treatment + age + gender + log(totalReads)


# --- 2. AJUSTE DE MODELOS ---
# A) Ajuste del modelo completo: Dirichlet-Multinomial Generalizado (GDM)
cat("Ajustando el modelo GDM completo (puede tardar unos segundos)... \n")
modelo_gdm <- MGLMreg(formula = formula_conteo, data = rnaseq, dist = "GDM", LRT = TRUE)

# B) Ajuste del modelo restringido: Dirichlet-Multinomial estándar (DM)
cat("Ajustando el modelo DM restringido para la comparación... \n\n")
modelo_dm <- MGLMreg(formula = formula_conteo, data = rnaseq, dist = "DM", LRT = FALSE)


# --- 3. TEST DE RAZÓN DE VEROSIMILITUDES (GDM VS DM) ---
loglik_gdm <- logLik(modelo_gdm)
loglik_dm  <- logLik(modelo_dm)

# Cálculo del estadístico del contraste (G^2)
G2 <- 2 * (loglik_gdm - loglik_dm)

# CORRECCIÓN DE BUG DEL PAQUETE: Cálculo manual de los Grados de Libertad.
# MGLM suele errar en el conteo automático debido a la estructura indexada.
# Contamos el número exacto de parámetros estimados en las matrices de coeficientes:
k_gdm <- length(coef(modelo_gdm))  # Estructura doble (50 parámetros: Alfas y Betas)
k_dm  <- length(coef(modelo_dm))   # Estructura clásica (30 parámetros)
delta_df <- k_gdm - k_dm           # Diferencia exacta (20 grados de libertad)

# Obtención del p-valor mediante la distribución Chi-cuadrado clásica
p_valor_TRV <- pchisq(G2, df = delta_df, lower.tail = FALSE)


# --- 4. IMPRESIÓN DE MÉTRICAS GLOBALES DE SELECCIÓN ---
cat("==================================================\n")
cat("          RESULTADOS GLOBALES (CAPÍTULO 3)        \n")
cat("==================================================\n")
cat("Log-Verosimilitud Modelo Restringido (DM) : ", loglik_dm, "\n")
cat("Log-Verosimilitud Modelo Completo (GDM)   : ", loglik_gdm, "\n")
cat("Estadístico del Contraste (G^2)           : ", G2, "\n")
cat("Grados de Libertad del Test (DF)          : ", delta_df, "\n")
cat("p-valor del TRV (GDM vs DM)               : ", p_valor_TRV, "\n")
cat("Criterio de Información (AIC) - GDM       : ", AIC(modelo_gdm), "\n")
cat("Criterio de Información (BIC) - GDM       : ", BIC(modelo_gdm), "\n")
cat("==================================================\n\n")


# --- 5. EXTRACCIÓN DE ESTIMACIONES E INFERENCIA INDIVIDUAL ---
# CORRECCIÓN DE BUG DEL PAQUETE: El método print(modelo_gdm) genérico genera un 
# conflicto de dimensiones ('dimnames') debido a la secuencialidad del modelo.
# Se realiza una extracción manual directa y segura:

cat("==================================================\n")
cat("    MATRIZ DE COEFICIENTES ESTIMADOS (GDM)        \n")
cat("==================================================\n")
print(coef(modelo_gdm))
cat("\n")

cat("==================================================\n")
cat("     P-VALORES DEL TEST DE WALD (GDM)             \n")
cat("==================================================\n")
# Extracción robusta compatible con la naturaleza del objeto S4 de MGLM
if (isS4(modelo_gdm)) {
  print(modelo_gdm@wald.p)
} else {
  print(modelo_gdm$wald.p)
}
cat("==================================================\n")