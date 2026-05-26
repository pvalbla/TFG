# ==============================================================================
# CAPÍTULO 2: MODELO DIRICHLET-MULTINOMIAL ESTÁNDAR
# Objetivo: Estimación de la sobredispersión, contraste TRV y criterios de selección
# ==============================================================================

# --- 1. CONFIGURACIÓN INICIAL Y CARGA DE DATOS ---
library(dirmult)

# Carga del dataset empírico y selección de la matriz de conteos
data(us)
datos_locus <- us[[1]]  # Matriz de conteos del primer locus (D3S1358)

# --- 2. ESTIMACIÓN DEL MODELO DIRICHLET-MULTINOMIAL (COMPLETO) ---
# Ajuste por máxima verosimilitud (tolerancia epsilon fijada según la memoria)
modelo_dm <- dirmult(datos_locus, epsilon = 10^(-4), trace = FALSE)
loglik_dm <- modelo_dm$loglik

# Extracción de la tabla resumen con parámetros y errores estándar asintóticos
resumen_dm <- dirmult.summary(datos_locus, modelo_dm)

# --- 3. ESTIMACIÓN DEL MODELO MULTINOMIAL (RESTRINGIDO / NULO) ---
# Bajo H0 no hay sobredispersión (theta = 0). Las probabilidades esperadas
# convergen a las frecuencias relativas globales de la muestra.
pi_nulo <- colSums(datos_locus) / sum(datos_locus)

# Log-verosimilitud del modelo restringido para el contraste
loglik_nulo <- sum(datos_locus * log(pi_nulo), na.rm = TRUE)

# --- 4. CONTRASTE DE RAZÓN DE VEROSIMILITUDES (TRV) Y DIAGNÓSTICO ---
# Estadístico de contraste (G^2)
G2 <- 2 * (loglik_dm - loglik_nulo)

# P-valor corregido por contraste en la frontera del espacio paramétrico (theta >= 0)
# Distribución mixta asintótica de Self y Liang (1987)
p_valor_TRV <- 0.5 * pchisq(G2, df = 0, lower.tail = FALSE) + 
               0.5 * pchisq(G2, df = 1, lower.tail = FALSE)

# Cálculo del Criterio de Información de Akaike (AIC)
# k = (d - 1) parámetros de probabilidad + 1 parámetro de dispersión global (theta)
k_parametros <- ncol(datos_locus)
AIC_dm <- 2 * k_parametros - 2 * loglik_dm

# --- 5. IMPRESIÓN FORMAL DE RESULTADOS ---
cat("\n>>> ESTIMACIÓN PARAMÉTRICA (MODELO DM) <<<\n")
print(resumen_dm)

cat("\n>>> RESULTADOS GLOBALES Y DIAGNÓSTICO (CAPÍTULO 2) <<<\n")
cat("Log-Verosimilitud Mod. Restringido (Nulo) :", loglik_nulo, "\n")
cat("Log-Verosimilitud Mod. Completo (DM)      :", loglik_dm, "\n")
cat("Estadístico del Contraste (G^2)           :", G2, "\n")
cat("P-valor (Mixtura Chi-cuadrado)            :", p_valor_TRV, "\n")
cat("Criterio de Información (AIC)             :", AIC_dm, "\n")