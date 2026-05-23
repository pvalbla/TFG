# ==============================================================================
# CAPÍTULO 1: MODELO MULTINOMIAL CLÁSICO
# Selección de variables (Backward Elimination) y predicción final
# ==============================================================================

# --- 1. CONFIGURACIÓN INICIAL Y CARGA DE DATOS ---
library(mlogit)
data("Heating", package = "mlogit")

# --- 2. PREPARACIÓN DE LOS DATOS ---
# Conversión al formato indexado long/wide requerido por mlogit
H <- dfidx(Heating, choice = "depvar", varying = c(3:12))

# --- 3. PROCESO DE SELECCIÓN HACIA ATRÁS (BACKWARD ELIMINATION) ---

# Paso 1: Ajuste del modelo completo saturado
modelo_completo <- mlogit(depvar ~ 0 | income + agehed + rooms, data = H)
cat("\n>>> RESUMEN: MODELO COMPLETO <<<\n")
print(summary(modelo_completo))

# Paso 2: Exclusión de la variable menos significativa ('income') y reajuste
modelo_paso2 <- mlogit(depvar ~ 0 | agehed + rooms, data = H)
cat("\n>>> RESUMEN: PASO 2 (SIN INCOME) <<<\n")
print(summary(modelo_paso2))

# Paso 3: Exclusión de la siguiente variable no significativa ('rooms') -> Modelo óptimo
modelo_final <- mlogit(depvar ~ 0 | agehed, data = H)
cat("\n>>> RESUMEN: MODELO FINAL DEPURADO (SOLO AGEHED) <<<\n")
print(summary(modelo_final))

# --- 4. DIAGNÓSTICO Y COMPARATIVA DE MODELOS ---

# Evaluación global mediante Criterio de Información de Akaike (AIC)
cat("\n>>> COMPARATIVA DE DIAGNÓSTICO (AIC) <<<\n")
cat("Modelo Completo :", AIC(modelo_completo), "\n")
cat("Modelo Paso 2   :", AIC(modelo_paso2), "\n")
cat("Modelo Final    :", AIC(modelo_final), "\n")

# Intervalos de confianza asintóticos (95%) del modelo parsimonioso
cat("\n>>> INTERVALOS DE CONFIANZA - MODELO FINAL (95%) <<<\n")
print(round(confint(modelo_final), 4))

# --- 5. VERIFICACIÓN METODOLÓGICA DE LOS INTERVALOS DE CONFIANZA ---
# Obtención manual basada en la estimación y la inversa de la matriz de Fisher
cat("\n>>> DESGLOSE ASINTÓTICO DEL IC (ALTERNATIVA ER) <<<\n")

# Extracción de coeficientes y errores estándar en el paso de su respectiva evaluación
beta_inc_er  <- coef(modelo_completo)["income:er"]
ee_inc_er    <- sqrt(diag(vcov(modelo_completo)))["income:er"]

beta_room_er <- coef(modelo_paso2)["rooms:er"]
ee_room_er   <- sqrt(diag(vcov(modelo_paso2)))["rooms:er"]

beta_age_er  <- coef(modelo_final)["agehed:er"]
ee_age_er    <- sqrt(diag(vcov(modelo_final)))["agehed:er"]

z_critico    <- qnorm(0.975) # Cuantil normal estándar

cat("IC 95% Income (er) [Mod. Completo]: [", 
    round(beta_inc_er - z_critico * ee_inc_er, 4), ", ", 
    round(beta_inc_er + z_critico * ee_inc_er, 4), "]\n")

cat("IC 95% Rooms (er)  [Mod. Paso 2]:   [", 
    round(beta_room_er - z_critico * ee_room_er, 4), ", ", 
    round(beta_room_er + z_critico * ee_room_er, 4), "]\n")

cat("IC 95% Agehed (er) [Mod. Final]:    [", 
    round(beta_age_er - z_critico * ee_age_er, 4), ", ", 
    round(beta_age_er + z_critico * ee_age_er, 4), "]\n")

# --- 6. SIMULACIÓN Y PREDICCIÓN DE PROBABILIDADES ---
# Análisis de escenarios para evaluar el impacto aislado de la edad (30 vs 70 años)
casas_simuladas <- Heating[1:2, ]
casas_simuladas$agehed[1] <- 30  # Escenario A: Decisor joven
casas_simuladas$agehed[2] <- 70  # Escenario B: Decisor mayor

H_simulado <- dfidx(casas_simuladas, choice = "depvar", varying = c(3:12))

# Generación de cuotas de mercado condicionales con el modelo final óptimo
probabilidades <- predict(modelo_final, newdata = H_simulado)
rownames(probabilidades) <- c("Familia Joven (30 años)", "Familia Mayor (70 años)")

cat("\n>>> PROBABILIDADES PREDICHAS DEL SISTEMA (%) <<<\n")
print(round(probabilidades * 100, 2))
