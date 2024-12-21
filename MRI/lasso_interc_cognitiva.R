# LASSO implementation fatiga cognitiva
rm(list = ls())

library(rjags)           # rjags library allows R to interface with JAGS 
library(coda)            # coda package provides tools for summarizing and visualizing MCMC output 
library(ggmcmc)          # ggmcmc is used for diagnostics of MCMC chains and plots 
library(MASS)            # MASS for generating multivariate normal data 
library(Matrix)          # Matrix package is for matrix computations, particularly for sparse matrices 

setwd("~/Documents/GitHub/Figueroa-Vargas_Navarrete_2025_Scientific_Reports/MRI") # WD 
source("HDIofMCMC.r") 
# Cargar datos
DATA <- read.csv("../DATA/data_seg_GM_WM_test.csv")

# Verificar las dimensiones del objeto original "DATA"
dime=dim(DATA)  # Esto debería mostrar 'n x m'

load("../DATA/pca_20_components_grupal_wm.RData") # wm 20 componentes

# Modelo en JAGS
model_string_Lasso <- " 
model { 
        # Likelihood 
        for (i in 1:N) { 
            y[i] ~ dnorm(mu[i], pre_sig2)  # Distribución normal para TOTAL
            mu[i] <- beta0 + sum(beta[] * X[i,])  # Modelo para la media
        } 
 
        # Prior for beta0 (intercept) 
        beta0 ~ dnorm(0, 1.0E-2)  # Asumiendo un prior débilmente informativo 
 
        # Prior for sigma2 
        pre_sig2 ~ dgamma(0.01, 0.01)  # Prior gamma para la precisión 
        sigma2 <- 1 / pre_sig2  # Conversión de precisión a varianza
 
        # Priors for beta coefficients 
        for (j in 1:p) { 
            beta[j] ~ dnorm(0, tau_beta[j]) 
            tau_beta[j] <- 1 / (sigma2 * tau_sq[j]) 
            tau_sq[j] ~ dexp(lambda^2 / 2)  # Lasso penalty
        } 
} 
"

# Cargar la matriz de PCA (asumimos que se llama 'X' tras la carga)
X <- pca_20_components_grupal_wm # La matriz cargada desde el archivo
X = cbind(DATA$Grupo, X, DATA$EstimatedTotalIntraCranialVol)

# Crear variable de interacción (Grupo * PCA)
interaction_terms <- X[, -1] * DATA$Grupo  # Excluyendo la primera columna que es 'Grupo'
X <- cbind(X, interaction_terms,DATA$Fatiga_Fis)  # Unir los términos de interacción a la matriz de diseño

# Inicialización de valores para los parámetros
N = dim(X)[1] #(the number of observations) 
p = dim(X)[2]  #(the number of predictors)

# Variable dependiente
y <- as.numeric(DATA$Fatiga_Cog)

# Preparación de los datos para JAGS
lambda <- p * sqrt(var(lm(y ~ X - 1)$residuals)) / sum(abs(lm(y ~ X - 1)$coefficients))
data_lasso <- list(N = N, 
                   p = p, 
                   X = X, 
                   y = as.numeric(y), 
                   lambda = lambda,  
                   a = 2,  # Hyperparameter for the Beta distribution prior on theta 
                   b = p,  # Hyperparameter for the Beta distribution prior on theta 
                   lambda1 = 0.1,  # Hyperparameter for the lambda calculation 
                   lambda0 = 4) # Hyperparameter for the lambda calculation 

# Ejecutar el modelo en JAGS
fit <- jags.model(
  textConnection(model_string_Lasso),  # Especificación del modelo
  data = data_lasso,  # Los datos
  n.chains = 3,  # Número de cadenas de Markov
  n.adapt = 1000  # Número de iteraciones de adaptación
)

# Parámetros a monitorear
params <- c("beta0", "beta", "sigma2")
update(fit, 5000)  # Burn-in

# Muestra de las cadenas posteriori
samples <- coda.samples(fit, variable.names = params, n.iter = 20000)

# Guardar los resultados en 'model.res_cog'
model.res_cog <- as.mcmc(do.call(rbind, samples))

# Pruebas de Geweke
geweke.diag(model.res_cog)

# Convertir los resultados a formato 'ggmcmc'
ggmcmc_object <- ggs(model.res_cog)

# Gráfico caterpillar para visualizar las distribuciones posteriores de los parámetros
ggs_caterpillar(D = ggmcmc_object, family = "beta")

# Resumen de las muestras
summary(samples)

# Calcular la interacción entre el grupo de pacientes y controles
chains = cbind(samples[[1]], samples[[2]], samples[[3]])

# Ejemplo de cálculo de efectos de interacción
interaction_effect_patient = mean(chains[,"beta[37]"])  # Por ejemplo, beta para la interacción de grupo de pacientes
interaction_effect_control = mean(chains[,"beta[16]"])  # Por ejemplo, beta para la interacción de grupo de controles

# Imprimir resultados
cat("Efecto de interacción para pacientes:", interaction_effect_patient, "\n")
cat("Efecto de interacción para controles:", interaction_effect_control, "\n")

# Comprobar si los efectos son significativos
significant_interaction_patient = mean(chains[,"beta[37]"] < 0) * 2 ## component 15
significant_interaction_control = mean(chains[,"beta[16]"] > 0) * 2 ##

mean( (chains[,"beta[37]"] + chains[,"beta[16]"])< 0) * 2
HDIofMCMC((chains[,"beta[37]"] + chains[,"beta[16]"]))
mean((chains[,"beta[37]"] + chains[,"beta[16]"]))

mean( chains[,"beta[37]"] < 0) * 2
HDIofMCMC((chains[,"beta[37]"]))
mean((chains[,"beta[37]"]))

cat("Significancia de la interacción en pacientes (p-valor):", significant_interaction_patient, "\n")


##plot



library(patchwork)
library(ggcharts)
source("HDIofMCMC.r") 
load("../DATA/PCA_sorted.RData")
data_summary <- function(x) {
  m <- median(x)
  hdi = HDIofMCMC(x , credMass=0.95)
  ymin <- hdi[1]
  ymax <-  hdi[2]
  return(c(y=m,ymin=ymin,ymax=ymax))
}

pv2str =  function(pval){
  if (pval<0.001){
    AST="***"} else {if (pval<0.01){
      AST="**"} else {if (pval<0.05){
        AST="*"}else{AST="n.s."}}
    }  
}

# Cargar las librerías necesarias
library(ggplot2)
library(reshape2)
#install.packages("remotes")
#remotes::install_github("LCBC-UiO/ggseg")

library(ggseg)
library(ggseg3d)
library(dplyr)
library(tidyr)


regions = names(pca_loadings_sorted_wm)
regions = gsub("WMparc_wm_", "", regions)
regions = gsub("WMparc_wm_", "", regions)
regions = gsub("WMparc_wm_", "", regions)

nn=20
paraplot = data.frame(label=regions[1:nn], 
                      loading = pca_loadings_sorted_wm[1:nn]
)
regions = names(pca_loadings_sorted_wm8)
regions = gsub("WMparc_wm_", "", regions)
regions = gsub("WMparc_wm_", "", regions)
regions = gsub("WMparc_wm_", "", regions)

nn=20
paraplot8 = data.frame(label=regions[1:nn], 
                      loading = pca_loadings_sorted_wm8[1:nn]
)
regions = names(pca_loadings_sorted_wm13)
regions = gsub("WMparc_wm_", "", regions)
regions = gsub("WMparc_wm_", "", regions)
regions = gsub("WMparc_wm_", "", regions)

nn=20
paraplot13 = data.frame(label=regions[1:nn], 
                       loading = pca_loadings_sorted_wm13[1:nn]
)





# Suponiendo que ya tienes tus datos
interaction_patient15 <- chains[,"beta[37]"] 
beta_15 <- chains[,"beta[16]"] 
interaction_patient13 <- chains[,"beta[35]"] 
beta_13 <- chains[,"beta[14]"] 
interaction_patient8 <- chains[,"beta[30]"] 
beta_8 <- chains[,"beta[9]"] 
# Crear un data frame con los datos
data <- data.frame(
  value = c(interaction_patient15, beta_15),
  variable = factor(rep(c("MS i", "HC"), each = length(interaction_patient15)))
  
)
# Crear el gráfico de violín
AA= ggplot(data, aes(x = variable, y = value,color=variable, fill=variable)) +
  geom_violin(trim = FALSE) +
  coord_flip() + # Para hacerlo horizontal
  labs(title = "PC15 WM",
       x = "Lasso Regressors",
       y = "Posterior Distribution") +
  #theme(legend.position = "none") +
  geom_hline(yintercept=0,col="red") +
  stat_summary(fun.data=data_summary, geom="pointrange", color="black")+
  annotate(geom="text", y=mean(interaction_patient15)*1.0, x=2.4, 
           label=pv2str(mean(interaction_patient15<0)*2),
           color="black") +
  annotate(geom="text", y=mean(beta_15)*1.0, x=1.5, 
           label=pv2str(mean(beta_15>0)*2),
           color="black")



# Suponiendo que ya tienes tus datos

# Crear un data frame con los datos
data <- data.frame(
  value = c(interaction_patient13, beta_13),
  variable = factor(rep(c("MS i", "HC"), each = length(interaction_patient13)))
  
)
# Crear el gráfico de violín
BB = ggplot(data, aes(x = variable, y = value,color=variable, fill=variable)) +
  geom_violin(trim = FALSE) +
  coord_flip() + # Para hacerlo horizontal
  labs(title = "PC13 WM",
       x = "Lasso Regressors",
       y = "Posterior Distribution") +
  scale_fill_discrete(name = "Groups") + 
  scale_color_discrete(name = "Groups") + 
  # Cambiar el título de la leyenda
   theme(legend.position = "none") +
  geom_hline(yintercept=0,col="red") +
  stat_summary(fun.data=data_summary, geom="pointrange", color="black")+
  annotate(geom="text", y=mean(interaction_patient13)*1, x=2.4, 
           label=pv2str(mean(interaction_patient13<0)*2),
           color="black") +
  annotate(geom="text", y=mean(beta_13)*1, x=1.5, 
           label=pv2str(mean(beta_13>0)*2),
           color="black")

# Crear un data frame con los datos
data <- data.frame(
  value = c(interaction_patient8, beta_8),
  variable = factor(rep(c("MS i", "HC"), each = length(interaction_patient8)))
  
)
# Crear el gráfico de violín
CC = ggplot(data, aes(x = variable, y = value,color=variable, fill=variable)) +
  geom_violin(trim = FALSE) +
  coord_flip() + # Para hacerlo horizontal
  labs(title = "PC8 WM",
       x = "Lasso Regressors",
       y = "Posterior Distribution") +
  scale_fill_discrete(name = "Groups") + 
  scale_color_discrete(name = "Groups") + 
  # Cambiar el título de la leyenda
   theme(legend.position = "none") +
  geom_hline(yintercept=0,col="red") +
  stat_summary(fun.data=data_summary, geom="pointrange", color="black")+
  annotate(geom="text", y=mean(interaction_patient8)*1, x=2.4, 
           label=pv2str(mean(interaction_patient8<0)*2),
           color="black") +
  annotate(geom="text", y=mean(beta_8)*1, x=1.5, 
           label=pv2str(mean(beta_8>0)*2),
           color="black")


#colours = rev(heat.colors(10))
# para plot esta en el 
DD = ggseg(paraplot8,  colour = "gray", mapping = aes(fill = loading)) + #, position = "stacked"
  scale_fill_gradient(low="white", high="red", na.value = "white",limits = c(0.15, 0.32)) +  # Paleta de colores "hot"
  theme_minimal() + 
  labs(fill = "Loading values", title = "White Matter Loading for PC8")
EE = ggseg(paraplot13,  colour = "gray", mapping = aes(fill = loading)) + #, position = "stacked"
  scale_fill_gradient(low="white", high="red", na.value = "white",limits = c(0.15, 0.32)) +  # Paleta de colores "hot"
  theme_minimal() + 
  labs(fill = "Loading values", title = "White Matter Loading for PC13")
layout <- "
AABBCC
DDDDDD
DDDDDD
EEEEEE
EEEEEE"

CC +  BB + AA +DD+ EE + plot_layout(design = layout) +plot_annotation(tag_levels = list(c('A','B','C','D','E')))


####  Control con GAD& y PH9


library(readr)
Datos_paper <- read_delim("Datos_paper.csv", 
                          delim = ";", escape_double = FALSE, trim_ws = TRUE)

md = glm(PHQ9 ~  Sexo + Edad + Escolaridad + TOTAL, data = Datos_paper )

summary(md)

# Filtrar las filas sin valores faltantes en PHQ9
Datos_completos <- Datos_paper[!is.na(Datos_paper$PHQ9), ]

# Ajustar el modelo con los datos completos
md <- glm(PHQ9 ~ Sexo + Edad + Escolaridad + TOTAL, data = Datos_completos)

predicciones <- predict(md, newdata = Datos_paper)

# Redondear las predicciones al entero más cercano
predicciones_redondeadas <- round(predicciones)

# Reemplazar los valores NA en PHQ9 con las predicciones redondeadas
Datos_paper$PHQ9[is.na(Datos_paper$PHQ9)] <- predicciones_redondeadas[is.na(Datos_paper$PHQ9)]

# Ajustar el modelo con los datos completos
md <- glm(GAD7 ~  Sexo + Edad + Escolaridad + TOTAL, data = Datos_completos)

predicciones <- predict(md, newdata = Datos_paper)

# Redondear las predicciones al entero más cercano
predicciones_redondeadas <- round(predicciones)

# Reemplazar los valores NA en PHQ9 con las predicciones redondeadas
Datos_paper$GAD7[is.na(Datos_paper$GAD7)] <- predicciones_redondeadas[is.na(Datos_paper$GAD7)]


Datos_paper$S = Datos_paper$Sujeto
tabla_unida <- merge(DATA, Datos_paper, by = "S")

##

X2 <- cbind(X, Datos_paper$PHQ9, Datos_paper$GAD7) 
# Inicialización de valores para los parámetros
N = dim(X2)[1] #(the number of observations) 
p = dim(X2)[2]  #(the number of predictors)

data_lasso <- list(N = N, 
                   p = p, 
                   X = X2, 
                   y = as.numeric(y), 
                   lambda = lambda,  
                   a = 2,  # Hyperparameter for the Beta distribution prior on theta 
                   b = p,  # Hyperparameter for the Beta distribution prior on theta 
                   lambda1 = 0.1,  # Hyperparameter for the lambda calculation 
                   lambda0 = 4) # Hyperparameter for the lambda calculation 

# Ejecutar el modelo en JAGS
fit <- jags.model(
  textConnection(model_string_Lasso),  # Especificación del modelo
  data = data_lasso,  # Los datos
  n.chains = 3,  # Número de cadenas de Markov
  n.adapt = 1000  # Número de iteraciones de adaptación
)

# Ejecutar el modelo en JAGS
fit <- jags.model(
  textConnection(model_string_Lasso),  # Especificación del modelo
  data = data_lasso,  # Los datos
  n.chains = 3,  # Número de cadenas de Markov
  n.adapt = 1000  # Número de iteraciones de adaptación
)

# Parámetros a monitorear
params <- c("beta0", "beta", "sigma2")
update(fit, 5000)  # Burn-in

# Muestra de las cadenas posteriori
samples <- coda.samples(fit, variable.names = params, n.iter = 20000)

# Guardar los resultados en 'model.res'
model.res <- as.mcmc(do.call(rbind, samples))

# Pruebas de Geweke
geweke.diag(model.res)

# Convertir los resultados a formato 'ggmcmc'
ggmcmc_object <- ggs(model.res)

# Gráfico caterpillar para visualizar las distribuciones posteriores de los parámetros
ggs_caterpillar(D = ggmcmc_object, family = "beta")

# Resumen de las muestras
summary(samples)

# Calcular la interacción entre el grupo de pacientes y controles
chains = cbind(samples[[1]], samples[[2]], samples[[3]])

# Ejemplo de cálculo de efectos de interacción
interaction_effect_patient = mean(chains[,"beta[37]"])  # Por ejemplo, beta para la interacción de grupo de pacientes
interaction_effect_control = mean(chains[,"beta[16]"])  # Por ejemplo, beta para la interacción de grupo de controles


# Imprimir resultados
cat("Efecto de interacción para pacientes:", interaction_effect_patient, "\n")
cat("Efecto de interacción para controles:", interaction_effect_control, "\n")

# Comprobar si los efectos son significativos
significant_interaction_patient = mean(chains[,"beta[37]"] < 0) * 2 ## component 15
significant_interaction_control = mean(chains[,"beta[16]"] > 0) * 2 ##

significant_patientes = mean((chains[,"beta[16]"] + chains[,"beta[37]"]   ) < 0) * 2 ##

cat("Significancia de la interacción en pacientes (p-valor):", significant_interaction_patient, "\n")
cat("Significancia de la interacción en controles (p-valor):", significant_interaction_control, "\n")



