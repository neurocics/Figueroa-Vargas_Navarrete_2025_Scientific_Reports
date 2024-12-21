# 
# Figueroa-Vargas, Navarrete et al., 2025
#    
#
# Cargar las librerías necesarias
#library(tidyverse)  # Para la manipulación de datos

rm(list = ls())
library(stats)      # Para el PCA
library(ggplot2)    # Para la visualización
library(dplyr)



setwd("~/Documents/GitHub/Figueroa-Vargas_Navarrete_2025_Scientific_Reports/MRI") # WD 

# Cargar datos
DATA <- read.csv("../DATA/data_seg_GM_WM_test.csv")

str(DATA)  # Muestra la estructura de 'DATA' incluyendo el tipo de cada columna



# Verificar las dimensiones del objeto original "DATA"
dime=dim(DATA)  # Esto debería mostrar 'n x m'

####################### PCA para Seg subcortical ####################

# Seleccionar solo las columnas que comienzan con "aseg_"
Data_aseg <- DATA %>% select(starts_with("aseg_"))


Data_aseg_filtered <- Data_aseg %>%
  select(-contains("Ventricle"),
         -contains("Vent"),
         -contains("WM"),
         -contains("White"),
         -contains("CC"),
         -contains("optic"),
         -contains("CSF"),
         -contains("vessel"),
         -contains("plexus"),
         -contains("hypointensities"))


# Definir la columna que quieres "regress out"
columna_a_regresar <- DATA$EstimatedTotalIntraCranialVol

# Aplicar la regresión de columna_x sobre cada columna y obtener los residuos
for (col_name in colnames(Data_aseg)) {
  # Ajustar el modelo de regresión lineal
  model <- lm(Data_aseg[[col_name]] ~ columna_a_regresar)
  
  # Reemplazar la columna con los residuos (efecto de columna_x eliminado)
  Data_aseg[[col_name]] <- resid(model)
}

# Mostrar la tabla resultante
print(Data_aseg_filtered)

# Guardar la nueva tabla en un archivo CSV (opcional)
write.csv(Data_aseg_filtered, "Data_aseg.csv", row.names = FALSE)

# Normalizar todas las columnas de Data_aseg_filtered
Data_aseg_scaled <- as.data.frame(scale(Data_aseg_filtered))

# Mostrar la tabla normalizada
print(Data_aseg_scaled)

# Guardar la nueva tabla normalizada en un archivo CSV (opcional)
write.csv(Data_aseg_scaled, "../DATA/Data_aseg_scaled.csv", row.names = FALSE)

# Realizar el análisis de componentes principales (PCA)
pca_result_aseg <- prcomp(Data_aseg_scaled, center = TRUE, scale. = TRUE)

# Resumir los resultados del PCA
summary(pca_result_aseg)

# Obtener los primeros 10 componentes
pca_10_components_grupal_aseg <- pca_result_aseg$x[, 1:10]

# Guardar la matriz en un archivo RData
save(pca_10_components_grupal_aseg, file = "pca_10_components_grupal_aseg.RData")

# Visualizar la varianza explicada por los primeros 10 componentes
explained_variance <- summary(pca_result_aseg)$importance[2, 1:10] * 100
explained_variance_df <- data.frame(
  Component = paste0("PC", 1:10),
  Variance = explained_variance
)

ggplot(explained_variance_df, aes(x = Component, y = Variance)) +
  geom_bar(stat = "identity") +
  xlab("Principal Components") +
  ylab("Explained Variance (%)") +
  ggtitle("Explained Variance by the First 10 Principal Components") +
  theme_minimal()

# Guardar los resultados en un archivo CSV (opcional)
#write.csv(pca_20_components_grupal, "PCA_20_components_grupal.csv", row.names = FALSE)

# Obtener los loadings (rotaciones) para PC1
pca_loadings_aseg <- pca_result_aseg$rotation[, 1]

# Ordenar los loadings por su magnitud
pca_loadings_sorted_aseg <- sort(pca_loadings_aseg, decreasing = FALSE)

# Visualizar las contribuciones en un gráfico de barras
barplot(pca_loadings_sorted_aseg, main = "Contribution of Cortical Regions to PC1", 
        xlab = "Regions", ylab = "Loading Value", las = 2, cex.names = 0.7)

# Imprimir las regiones que más contribuyen al PC1 aseg 
print(pca_loadings_sorted_aseg[1:10])  # Imprime las 10 regiones con mayor contribución







###### crear PCA solo sustancia Blanca ###########



Data_WM <- DATA %>%
  select(starts_with("WMparc_"),
         contains("aseg_CC"),
         contains("Cerebellum_White_Matter"),
         contains("Chiasm"))

# Normalizar los datos para las columnas 10 a 108
nombres = names(Data_WM)

data_scaled_wm <- scale(Data_WM)

# Verificar las dimensiones de la matriz escalada
dim(data_scaled_wm) 

# Realizar el análisis de componentes principales (PCA)
pca_result_wm <- prcomp(data_scaled_wm, center = TRUE, scale. = TRUE)

# Resumir los resultados del PCA
summary(pca_result_wm)

# Obtener los primeros 20 componentes
pca_20_components_grupal_wm <- pca_result_wm$x[, 1:20]

# Guardar la matriz en un archivo RData
save(pca_20_components_grupal_wm, file = "pca_20_components_grupal_wm.RData")

# Visualizar la varianza explicada por los primeros 10 componentes WM
explained_variance <- summary(pca_result_wm)$importance[2, 1:20] * 100
explained_variance_df <- data.frame(
  Component = paste0("PC", 1:20),
  Variance = explained_variance
)

ggplot(explained_variance_df, aes(x = Component, y = Variance)) +
  geom_bar(stat = "identity") +
  xlab("Principal Components") +
  ylab("Explained Variance (%)") +
  ggtitle("Explained Variance by the First 10 Principal Components") +
  theme_minimal()

# Guardar los resultados en un archivo CSV (opcional)
#write.csv(pca_20_components_grupal, "PCA_20_components_grupal.csv", row.names = FALSE)

# Obtener los loadings (rotaciones) para PC1
pca_loadings_wm <- pca_result_wm$rotation[, 1]

# Ordenar los loadings por su magnitud
pca_loadings_sorted_wm <- sort(pca_loadings_wm, decreasing = FALSE)

# Visualizar las contribuciones en un gráfico de barras wm
barplot(pca_loadings_sorted_wm, main = "Contribution of Cortical Regions to PC1", 
        xlab = "Regions", ylab = "Loading Value", las = 2, cex.names = 0.7)

# Imprimir las regiones que más contribuyen al PC1 WM 
print(pca_loadings_sorted_wm[1:20])  # Imprime las 10 regiones con mayor contribución



# Obtener los loadings (rotaciones) para PC15
pca_loadings_wm <- pca_result_wm$rotation[, 15]

# Ordenar los loadings por su magnitud
pca_loadings_sorted_wm <- sort(abs(pca_loadings_wm), decreasing = TRUE)

# Visualizar las contribuciones en un gráfico de barras wm
barplot(pca_loadings_sorted_wm, main = "Contribution of Cortical Regions to PC15", 
        xlab = "Regions", ylab = "Loading Value", las = 2, cex.names = 0.7)

# Imprimir las regiones que más contribuyen al PC15 WM 
print(pca_loadings_sorted_wm[1:20])  # Imprime las 10 regiones con mayor contribución
df <- data.frame(
  #Name = names(pca_loadings_sorted_wm[1:20]), # Obtiene los nombres
  Value = pca_loadings_sorted_wm[1:20]        # Obtiene los valores
)

# Imprimir el data frame
print(df)

# Obtener los loadings (rotaciones) para PC8
pca_loadings_wm8 <- pca_result_wm$rotation[, 8]

# Ordenar los loadings por su magnitud
pca_loadings_sorted_wm8 <- sort(abs(pca_loadings_wm8), decreasing = TRUE)

# Visualizar las contribuciones en un gráfico de barras wm
barplot(pca_loadings_sorted_wm8, main = "Contribution of Cortical Regions to PC15", 
        xlab = "Regions", ylab = "Loading Value", las = 2, cex.names = 0.7)

# Imprimir las regiones que más contribuyen al PC15 WM 
print(pca_loadings_sorted_wm8[1:20])  # Imprime las 10 regiones con mayor contribución
df <- data.frame(
  #Name = names(pca_loadings_sorted_wm[1:20]), # Obtiene los nombres
  Value = pca_loadings_sorted_wm8[1:20]        # Obtiene los valores
)

# Imprimir el data frame
print(df)


# Obtener los loadings (rotaciones) para PC8
pca_loadings_wm13 <- pca_result_wm$rotation[, 13]

# Ordenar los loadings por su magnitud
pca_loadings_sorted_wm13 <- sort(abs(pca_loadings_wm13), decreasing = TRUE)

# Visualizar las contribuciones en un gráfico de barras wm
barplot(pca_loadings_sorted_wm13, main = "Contribution of Cortical Regions to PC15", 
        xlab = "Regions", ylab = "Loading Value", las = 2, cex.names = 0.7)

df <- data.frame(
  #Name = names(pca_loadings_sorted_wm[1:20]), # Obtiene los nombres
  Value = pca_loadings_sorted_wm13[1:20]        # Obtiene los valores
)

# Imprimir el data frame
print(df)# Imprimir las regiones que más contribuyen al PC15 WM 

###### crear PCA solo sustancia gris ###########



Data_GM <- DATA %>% select(starts_with("DKT"))

# Normalizar los datos para las columnas 10 a 108
nombres = names(Data_GM)

data_scaled_gm <- scale(Data_GM)

# Verificar las dimensiones de la matriz escalada
dim(data_scaled_gm) 

# Realizar el análisis de componentes principales (PCA)
pca_result_gm <- prcomp(data_scaled_gm, center = TRUE, scale. = TRUE)

# Resumir los resultados del PCA
summary(pca_result_gm)

# Obtener los primeros 20 componentes
pca_20_components_grupal_gm <- pca_result_gm$x[, 1:20]

# Guardar la matriz en un archivo RData
save(pca_20_components_grupal_gm, file = "pca_20_components_grupal_gm.RData")

# Visualizar la varianza explicada por los primeros 10 componentes GM
explained_variance <- summary(pca_result_gm)$importance[2, 1:20] * 100
explained_variance_df <- data.frame(
  Component = paste0("PC", 1:20),
  Variance = explained_variance
)

ggplot(explained_variance_df, aes(x = Component, y = Variance)) +
  geom_bar(stat = "identity") +
  xlab("Principal Components") +
  ylab("Explained Variance (%)") +
  ggtitle("Explained Variance by the First 10 Principal Components") +
  theme_minimal()

# Guardar los resultados en un archivo CSV (opcional)
#write.csv(pca_20_components_grupal, "PCA_20_components_grupal.csv", row.names = FALSE)

# Obtener los loadings (rotaciones) para PC1
pca_loadings_gm <- pca_result_gm$rotation[, 1]

# Ordenar los loadings por su magnitud
pca_loadings_sorted_gm <- sort(pca_loadings_gm, decreasing = FALSE)

# Visualizar las contribuciones en un gráfico de barras gm
barplot(pca_loadings_sorted_gm, main = "Contribution of Cortical Regions to PC1", 
        xlab = "Regions", ylab = "Loading Value", las = 2, cex.names = 0.7)

# Imprimir las regiones que más contribuyen al PC1 GM 
print(pca_loadings_sorted_gm[1:20])  # Imprime las 10 regiones con mayor contribución

save.image(file="../DATA/PCA_sorted.RData")

