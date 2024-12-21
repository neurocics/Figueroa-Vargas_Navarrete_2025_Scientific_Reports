# modelo regresion lineal de FA y tractos

#install.packages("readxl")   # en caso de no estar instalada
library(readxl)

setwd("~/Documents/GitHub/Figueroa-Vargas_Navarrete_2025_Scientific_Reports/MRI")

# anotar nombre del archivo y la hoja de los datos
datos <- read_excel("../DATA/Pacientes_ROI_data.xlsx", sheet = "Mean FA")
names(datos) <- c("...1"   ,                    "lh.caudalanteriorcingulate",
     "lh.parstriangularis"   ,     "Right.Amygdala"  ,          
      "lh.transversetemporal"   ,   "rh.temporalpole"  ,         
       "lh.entorhinal"      ,        "rh.inferiortemporal"   ,    
        "rh.precentral"    ,          "lh.paracentral"     ,       
          "lh.posteriorcingulate"     , "lh.caudalmiddlefrontal"  ,  
         "Right.Accumbens.area",       "lh.parsopercularis"   ,     
          "rh.entorhinal"   ,           "rh.lateraloccipital"     ,  
         "rh.frontalpole"   ,          "lh.isthmuscingulate"       ,
              "rh.bankssts"    ,            "Rh.medialorbitofrontal"    ,
           "fatigue"        ,            "SDMT"   ,                   
             "Mean"                      )
# Realiza la regresión lineal utilizando Mean para predecir fatigue
modelo <- lm(fatigue ~ lh.caudalanteriorcingulate +
               lh.parstriangularis +lh.posteriorcingulate+rh.bankssts+
               lh.transversetemporal, data = datos[datos$lh.caudalanteriorcingulate != 0 ,])


datos$lh.caudalanteriorcingulate == 0
datos$lh.parstriangularis == 0
datos$lh.entorhinal == 0
datos$lh.posteriorcingulate== 0



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
model_string_SSL <- " 
model { 
   #----------------- 
   #Likelihood 
   #----------------- 
  for (i in 1:N) { 
    y[i] ~ dnorm(mu[i], pre_sig2) 
    mu[i] <- beta0 + sum(beta[] * X[i,]) 
  } 
 
  # Priors for beta and reparameterization for Laplace prior 
  for (i in 1:p) { 
    # Prior for gamma (Bernoulli distribution) 
    gamma[i] ~ dbern(theta)   # Binary variable following Bernoulli distribution. 
 
    # Definition of lambda 
    lambda[i] <- gamma[i] * lambda1 + (1 - gamma[i]) * lambda0   #Mixing two lambda values based on gamma 
 
    # Prior for beta (Normal distribution with mean 0 and precision tau) 
    tau_sq[i] ~ dexp(lambda[i]^2 / 2) 
    beta[i] ~ dnorm(0, tau_beta[i]) 
        tau_beta[i] <- 1 / (sigma2 * tau_sq[i]) 
  } 
 
  theta ~ dbeta(a, b)  # Beta prior for theta (Bernoulli parameter) 
  beta0 ~ dnorm(0, 0.01)  # Normal prior for the intercept 
  pre_sig2 ~ dgamma(0.01,0.01)  # Precision follows a gamma distribution 
  sigma2 <- 1/pre_sig2 # Definition of sigma squared (variance) as the reciprocal of precision 
  } 
"


# Data preparation 



X = cbind(datos$lh.caudalanteriorcingulate,
          datos$rh.bankssts,
            datos$lh.parstriangularis,datos$lh.transversetemporal)

N = dim(X)[1] #(the number of observations) 
p = dim(X)[2]  #(the number of predictors)
X = X  #(the design matrix)

#y =  as.vector(pca_20_components[,21])
y <- as.numeric((datos$fatigue))


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
                   lambda0 = 4) # Hyperparameter for the lambda calculation )

# Ejecutar el modelo en JAGS
fit <- jags.model(
  textConnection(model_string_Lasso),  # Especificación del modelo
  data = data_lasso,  # Los datos
  #inits = inits_lasso,  # Valores iniciales para los parámetros
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

# Trace plot
# muy ento por los 20
#trace_plot <- traceplot(samples)

# Pruebas de Geweke
geweke.diag(model.res)

# Convertir los resultados a formato 'ggmcmc'
ggmcmc_object <- ggs(model.res)

# Gráfico caterpillar para visualizar las distribuciones posteriores de los parámetros
ggs_caterpillar(D = ggmcmc_object, family = "beta")

# Resumen de las muestras
summary(samples)

chains = cbind(samples[[1]],samples[[2]],samples[[3]])
mean(chains[,"beta[3]"]<0)*2
mean(chains[,"beta[3]"])
HDIofMCMC(chains[,"beta[3]"])

##
##

library(patchwork)
library(ggcharts)
source("HDIofMCMC.r") 
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


# Crear un data frame con los datos
data <- data.frame(
  value = c(chains[,"beta[3]"]),
  variable = factor(rep(c("IFG pars triangularis"), each = length(chains[,"beta[3]"])))
  
)
# Crear el gráfico de violín
AA= ggplot(data, aes(x = variable, y = value,color=variable, fill=variable)) +
  geom_violin(trim = FALSE) +
  
  #coord_flip() + # Para hacerlo horizontal
  labs(title = "FA ~ Fatigue",
       x = "SSL Regressors",
       y = "Posterior Distribution") +
  theme(legend.position = "none") +
  geom_hline(yintercept=0,col="red") +
  stat_summary(fun.data=data_summary, geom="pointrange", color="black")+
  annotate(geom="text", y=mean(chains[,"beta[3]"])*1.5, x=1.15, 
           label=pv2str(mean(chains[,"beta[3]"]<0)*2),
           color="black") 

AA


##

# Preparación de los datos para JAGS
#lambda <- p * sqrt(var(lm(y ~ X - 1)$residuals)) / sum(abs(lm(y ~ X - 1)$coefficients))
data_lasso <- list(N = N, 
                   p = p, 
                   X = X, 
                   y = as.numeric(y), 
                  # lambda = lambda,  
                   a = 2,  # Hyperparameter for the Beta distribution prior on theta 
                   b = p,  # Hyperparameter for the Beta distribution prior on theta 
                   lambda1 = 0.1,  # Hyperparameter for the lambda calculation 
                   lambda0 = 4) # Hyperparameter for the lambda calculation )

# Ejecutar el modelo en JAGS
fit <- jags.model(
  textConnection(model_string_SSL),  # Especificación del modelo
  data = data_lasso,  # Los datos
  #inits = inits_lasso,  # Valores iniciales para los parámetros
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

# Trace plot
# muy ento por los 20
#trace_plot <- traceplot(samples)

# Pruebas de Geweke
geweke.diag(model.res)

# Convertir los resultados a formato 'ggmcmc'
ggmcmc_object <- ggs(model.res)

# Gráfico caterpillar para visualizar las distribuciones posteriores de los parámetros
ggs_caterpillar(D = ggmcmc_object, family = "beta")

# Resumen de las muestras
summary(samples)

chains = cbind(samples[[1]],samples[[2]],samples[[3]])
mean(chains[,"beta[3]"]<0)*2
mean(chains[,"beta[3]"])
HDIofMCMC(chains[,"beta[3]"])

