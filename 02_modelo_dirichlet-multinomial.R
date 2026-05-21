# ==============================================================================
# CAPÍTULO 2: MODELO DIRICHLET-MULTINOMIAL ESTÁNDAR
# Objetivo: Estimación de la sobredispersión, contraste TRV y criterios de selección
# ==============================================================================

# --- 1. CONFIGURACIÓN INICIAL Y CARGA DE DATOS ---
library(dirmult)

# Carga del dataset empírico y selección de la matriz de conteos
data(us)
datos_locus <- us[[1]]  # Matriz de conteos del primer locus (D3S1358)


# --- 2. AJUSTE DEL MODELO COMPLETO (DIRICHLET-MULTINOMIAL) ---
# Estimación del modelo por máxima verosimilitud (epsilon según especificación técnica)
modelo_dm <- dirmult(datos_locus, epsilon = 10^(-4), trace = FALSE)
loglik_dm <- modelo_dm$loglik

# Extracción de la tabla resumen con errores estándar asintóticos
resumen_dm <- dirmult.summary(datos_locus, modelo_dm)


# --- 3. ESTIMACIÓN DEL MODELO RESTRINGIDO (MODELO NULO MULTINOMIAL) ---
# Bajo H0, no hay sobredispersión (theta = 0). Las probabilidades esperadas
# coinciden con las frecuencias relativas globales del locus.
pi_nulo <- colSums(datos_locus) / sum(datos_locus)

# Cálculo de la Log-verosimilitud del modelo restringido (escala comparable)
loglik_nulo <- sum(datos_locus * log(pi_nulo), na.rm = TRUE)


# --- 4. CONTRASTE DE HIPÓTESIS Y SELECCIÓN DE MODELO ---
# Cálculo del estadístico del Test de Razón de Verosimilitudes (G^2)
G2 <- 2 * (loglik_dm - loglik_nulo)

# Obtención del p-valor mediante la distribución mixta de Self y Liang (1987)
# Corrige el test al evaluar el parámetro theta en la frontera del espacio (theta >= 0)
p_valor_TRV <- 0.5 * pchisq(G2, df = 0, lower.tail = FALSE) + 
               0.5 * pchisq(G2, df = 1, lower.tail = FALSE)

# Criterio de Información de Akaike (AIC)
# k = (d - 1) probabilidades independientes + 1 parámetro de dispersión (theta)
k_parametros <- ncol(datos_locus)
AIC_dm <- 2 * k_parametros - 2 * loglik_dm


# --- 5. IMPRESIÓN FORMAL DE RESULTADOS ---
# Bloque de salida unificado con los marcadores clave para la redacción del TFG
cat("\n==================================================\n")
cat("          RESULTADOS DIAGNÓSTICOS (CAPÍTULO 2)    \n")
cat("==================================================\n")
cat("Log-Verosimilitud Modelo Restringido (Nulo) : ", loglik_nulo, "\n")
cat("Log-Verosimilitud Modelo Completo (DM)       : ", loglik_dm, "\n")
cat("Estadístico del Contraste (G^2)              : ", G2, "\n")
cat("p-valor del Test (Mixtura Chi-cuadrado)      : ", p_valor_TRV, "\n")
cat("Criterio de Información (AIC)                : ", AIC_dm, "\n")
cat("==================================================\n\n")