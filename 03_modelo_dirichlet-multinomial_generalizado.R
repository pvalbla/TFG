# ==============================================================================
# CAPÍTULO 3: MODELO DIRICHLET-MULTINOMIAL GENERALIZADO (GDM)
# Objetivo: Regresión multivariante, contraste de especificación y diagnóstico Wald
# ==============================================================================

# --- 1. CONFIGURACIÓN INICIAL Y CARGA DE DATOS ---
library(MGLM)

# Carga del dataset biológico (Conteo de lecturas de ARN-Seq)
data("rnaseq")

# Definición de la fórmula de regresión multivariante
# Respuesta: Matriz de conteos de los 6 exones
# Explicativas: Tratamiento, edad, género y log(lecturas totales)
formula_conteo <- cbind(X1, X2, X3, X4, X5, X6) ~ treatment + age + gender + log(totalReads)

# --- 2. ESTIMACIÓN DE MODELOS (COMPLETO Y RESTRINGIDO) ---
# Paso 1: Ajuste del modelo Dirichlet-Multinomial Generalizado (GDM)
cat("\n>>> AJUSTANDO MODELO GDM COMPLETO (Espere unos segundos...) <<<\n")
modelo_gdm <- MGLMreg(formula = formula_conteo, data = rnaseq, dist = "GDM", LRT = TRUE)

# Paso 2: Ajuste del modelo Dirichlet-Multinomial estándar (DM)
cat("\n>>> AJUSTANDO MODELO DM RESTRINGIDO PARA COMPARATIVA... <<<\n")
modelo_dm <- MGLMreg(formula = formula_conteo, data = rnaseq, dist = "DM", LRT = FALSE)

# --- 3. CONTRASTE DE RAZÓN DE VEROSIMILITUDES (GDM vs DM) ---
loglik_gdm <- logLik(modelo_gdm)
loglik_dm  <- logLik(modelo_dm)

# Estadístico del contraste (G^2)
G2 <- 2 * (loglik_gdm - loglik_dm)

# CORRECCIÓN METODOLÓGICA: Cálculo manual de los Grados de Libertad (delta_df)
# Se calcula la diferencia exacta de parámetros estimados en las matrices de coeficientes
k_gdm <- length(coef(modelo_gdm))  # Estructura GDM: 50 parámetros (Alfas y Betas)
k_dm  <- length(coef(modelo_dm))   # Estructura DM: 30 parámetros
delta_df <- k_gdm - k_dm           # Diferencia exacta: 20 grados de libertad

# P-valor asintótico (Distribución Chi-cuadrado clásica)
p_valor_TRV <- pchisq(G2, df = delta_df, lower.tail = FALSE)

# --- 4. DIAGNÓSTICO Y COMPARATIVA DE MODELOS ---
cat("\n>>> RESULTADOS GLOBALES Y CONTRASTE DE ESTRUCTURA <<<\n")
cat("Log-Verosimilitud Mod. Restringido (DM) :", loglik_dm, "\n")
cat("Log-Verosimilitud Mod. Completo (GDM)   :", loglik_gdm, "\n")
cat("Estadístico del Contraste (G^2)         :", G2, "\n")
cat("Grados de Libertad del Test (DF)        :", delta_df, "\n")
cat("P-valor del TRV (GDM vs DM)             :", p_valor_TRV, "\n")
cat("Criterio de Información (AIC) - GDM     :", AIC(modelo_gdm), "\n")
cat("Criterio de Información (BIC) - GDM     :", BIC(modelo_gdm), "\n")

# --- 5. INFERENCIA INDIVIDUAL (TEST DE WALD CONJUNTO) ---
cat("\n>>> MATRIZ DE COEFICIENTES ESTIMADOS <<<\n")
# Extracción directa para evitar conflictos de dimensiones ('dimnames') originados por el modelo secuencial
print(coef(modelo_gdm))

cat("\n>>> SIGNIFICATIVIDAD INDIVIDUAL: P-VALORES TEST DE WALD <<<\n")
# Extracción robusta compatible con la naturaleza del objeto S4/S3 del paquete MGLM
if (isS4(modelo_gdm)) {
  print(modelo_gdm@wald.p)
} else {
  print(modelo_gdm$wald.p)
}

# --- 6. INTERVALOS DE CONFIANZA ---
# Obtención manual basada en la estimación y los errores estándar asintóticos
cat("\n>>> DESGLOSE DEL IC PARA EL TRATAMIENTO (EXÓN X2) <<<\n")

# 1. Extracción para el parámetro Gamma (vínculo de la probabilidad)
# Fila 2 (treatment), Columna 2 (X2) en la matriz de coeficientes y SE
beta_gamma_trat_X2 <- coef(modelo_gdm)[2, 2]
ee_gamma_trat_X2   <- modelo_gdm@SE[2, 2]

# 2. Extracción para el parámetro Delta (vínculo de la dispersión)
# Fila 2 (treatment), Columna 7 (X2 en la segunda matriz)
# Nota: MGLM junta las dos matrices por columnas en @SE: cols 1:5 son Gamma, cols 6:10 son Delta
beta_delta_trat_X2 <- coef(modelo_gdm)[2, 7]
ee_delta_trat_X2   <- modelo_gdm@SE[2, 7]

z_critico <- qnorm(0.975)  # Cuantil normal estándar para el 95% de confianza

# Cálculo e impresión de los límites del intervalo
cat("IC 95% Gamma Tratamiento (X2): [", 
    round(beta_gamma_trat_X2 - z_critico * ee_gamma_trat_X2, 4), ", ", 
    round(beta_gamma_trat_X2 + z_critico * ee_gamma_trat_X2, 4), "]\n")

cat("IC 95% Delta Tratamiento (X2): [", 
    round(beta_delta_trat_X2 - z_critico * ee_delta_trat_X2, 4), ", ", 
    round(beta_delta_trat_X2 + z_critico * ee_delta_trat_X2, 4), "]\n")