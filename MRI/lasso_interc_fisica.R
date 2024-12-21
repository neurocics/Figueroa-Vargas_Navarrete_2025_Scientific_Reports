# LASSO implementation fatiga fisica
rm(list = ls())

library(rjags)           # rjags library allows R to interface with JAGS 
library(coda)            # coda package provides tools for summarizing and visualizing MCMC output 
library(ggmcmc)          # ggmcmc is used for diagnostics of MCMC chains and plots 
library(MASS)            # MASS for generating multivariate normal data 
library(Matrix)          # Matrix package is for matrix computations, particularly for sparse matrices 

setwd("~/Documents/GitHub/Figueroa-Vargas_Navarrete_2025_Scientific_Reports/MRI") # WD 

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
X <- cbind(X, interaction_terms)  # Unir los términos de interacción a la matriz de diseño

# Inicialización de valores para los parámetros
N = dim(X)[1] #(the number of observations) 
p = dim(X)[2]  #(the number of predictors)

# Variable dependiente
y <- as.numeric(DATA$Fatiga_Fis)

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

# Guardar los resultados en 'model.res_fis'
model.res_fis <- as.mcmc(do.call(rbind, samples))

# Pruebas de Geweke
geweke.diag(model.res_fis)

# Convertir los resultados a formato 'ggmcmc'
ggmcmc_object <- ggs(model.res_fis)

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

cat("Significancia de la interacción en pacientes (p-valor):", significant_interaction_patient, "\n")
cat("Significancia de la interacción en controles (p-valor):", significant_interaction_control, "\n")


##plot

library(patchwork)
library(ggcharts)
load("../DATA/PCA_sorted.RData")
source("~/Documents/GitHub/LU/Behave/HDIofMCMC.r") 
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

# Suponiendo que ya tienes tus datos
interaction_patient <- chains[,"beta[37]"] 
beta_16 <- chains[,"beta[16]"] 
interaction_patient1 <- chains[,"beta[23]"] 
beta_1 <- chains[,"beta[2]"] 
# Crear un data frame con los datos
data <- data.frame(
  value = c(interaction_patient, beta_16),
  variable = factor(rep(c("MS i", "HC"), each = length(interaction_patient)))
  
)
# Crear el gráfico de violín
AA= ggplot(data, aes(x = variable, y = value,color=variable, fill=variable)) +
  geom_violin(trim = FALSE) +
  coord_flip() + # Para hacerlo horizontal
  labs(title = "PC15 White matter",
       x = "Lasso Regressors",
       y = "Posterior Distribution") +
  theme(legend.position = "none") +
  geom_hline(yintercept=0,col="red") +
  stat_summary(fun.data=data_summary, geom="pointrange", color="black")+
  annotate(geom="text", y=mean(interaction_patient)*1.0, x=2.5, 
           label=pv2str(mean(interaction_patient<0)*2),
           color="black") +
  annotate(geom="text", y=mean(beta_16)*1.0, x=1.5, 
           label=pv2str(mean(beta_16>0)*2),
           color="black")



# Suponiendo que ya tienes tus datos
interaction_patient <- chains[,"beta[37]"] 
beta_16 <- chains[,"beta[16]"] 
interaction_patient1 <- chains[,"beta[23]"] 
beta_1 <- chains[,"beta[2]"] 
# Crear un data frame con los datos
data <- data.frame(
  value = c(interaction_patient1, beta_1),
  variable = factor(rep(c("MS i", "HC"), each = length(interaction_patient)))
  
)
# Crear el gráfico de violín
BB = ggplot(data, aes(x = variable, y = value,color=variable, fill=variable)) +
  geom_violin(trim = FALSE) +
  coord_flip() + # Para hacerlo horizontal
  labs(title = "PC1 White Matter",
       x = "Lasso Regressors",
       y = "Posterior Distribution") +
  scale_fill_discrete(name = "Groups") + 
  scale_color_discrete(name = "Groups") + 
  # Cambiar el título de la leyenda
  # theme(legend.position = "none") +
  geom_hline(yintercept=0,col="red") +
  stat_summary(fun.data=data_summary, geom="pointrange", color="black")+
  annotate(geom="text", y=mean(interaction_patient1)*1, x=2.5, 
           label=pv2str(mean(interaction_patient1<0)*2),
           color="black") +
  annotate(geom="text", y=mean(beta_1)*1, x=1.5, 
           label=pv2str(mean(beta_1>0)*2),
           color="black")

#colours = rev(heat.colors(10))
# para plot esta en el 
CC = ggseg(paraplot,  colour = "gray", mapping = aes(fill = loading)) + #, position = "stacked"
  scale_fill_gradient(low="white", high="red", na.value = "white") +  # Paleta de colores "hot"
  theme_minimal() + 
  labs(fill = "Loading values", title = "White Matter Loading for PC15")

layout <- "
AABB
AABB
CCCC
CCCC
CCCC"

AA + BB + CC + plot_layout(design = layout) +plot_annotation(tag_levels = list(c('A','B','C')))


