# ==============================================================================
# CAPÍTULO 1: MODELO MULTINOMIAL CLÁSICO
# Objetivo: Estimación, diagnóstico y predicción del sistema de calefacción
# ==============================================================================

# --- 1. CONFIGURACIÓN INICIAL Y CARGA DE DATOS ---
library(mlogit)

# Carga del dataset original de sistemas de calefacción
data("Heating", package = "mlogit")


# --- 2. PREPARACIÓN DE LOS DATOS ---
# Transformación del dataframe a formato indexado requerido por mlogit
H <- dfidx(Heating, choice = "depvar", varying = c(3:12))


# --- 3. AJUSTE DEL MODELO MULTINOMIAL CLÁSICO ---
# Se asume 0 atributos de las alternativas y se evalúan solo variables del individuo
modelo_clasico <- mlogit(depvar ~ 0 | income + agehed + rooms, data = H)


# --- 4. EVALUACIÓN Y DIAGNÓSTICO DEL MODELO ---
# Resumen general del modelo estimado (Test de Wald)
summary(modelo_clasico)

# Extracción de los Criterios de Selección de Modelos (AIC y BIC)
cat("\n=========================================\n")
cat("      CRITERIOS DE SELECCIÓN DE MODELOS    \n")
cat("=========================================\n")
cat("AIC:", AIC(modelo_clasico), "\n")
cat("BIC:", BIC(modelo_clasico), "\n\n")

# Cálculo de los Intervalos de Confianza al 95%
# Basados en la matriz de varianzas-covarianzas (Inversa de la Información de Fisher)
cat("=========================================\n")
cat("   INTERVALOS DE CONFIANZA (NIVEL 95%)   \n")
cat("=========================================\n")
intervalos <- confint(modelo_clasico)
print(round(intervalos, 4))
cat("\n")


# --- 5. SIMULACIÓN Y PREDICCIÓN DE PROBABILIDADES ---
# Propósito: Evaluar el impacto aislado de la edad del cabeza de familia

# 1. Selección de dos perfiles reales de la base de datos para la simulación
casas_simuladas <- Heating[1:2, ]

# 2. Control de variables: Se fijan los ingresos y habitaciones en sus valores medios
casas_simuladas$income <- mean(Heating$income)
casas_simuladas$rooms  <- mean(Heating$rooms)

# 3. Escenarios: Forzamos la diferencia de edad para comparar los perfiles
casas_simuladas$agehed[1] <- 30  # Escenario A: Cabeza de familia joven
casas_simuladas$agehed[2] <- 70  # Escenario B: Cabeza de familia mayor

# 4. Indexación del nuevo mini-dataset simulado
H_simulado <- dfidx(casas_simuladas, choice = "depvar", varying = c(3:12))

# 5. Obtención y formateo de las cuotas de probabilidad predichas
probabilidades <- predict(modelo_clasico, newdata = H_simulado)
rownames(probabilidades) <- c("Familia Joven (30 años)", "Familia Mayor (70 años)")

# Impresión final de resultados predictivos en formato porcentual
cat("=========================================\n")
cat("       PROBABILIDADES PREDICHAS (%)      \n")
cat("=========================================\n")
print(round(probabilidades * 100, 2))
cat("=========================================\n")