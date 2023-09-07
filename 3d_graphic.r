cost_function <- function(mu, reference_data, N0, N1, d0, fmutant, p) {
  commande <- sprintf("./atreyu_forward_simulator %e,%e %e %e %e %e 384 | cut -d' ' -f1 > TauxMutation2.txt", N0, N1, d0, mu, fmutant, p)
  system(commande)
  generated_data <- scan("TauxMutation2.txt", quiet = TRUE)
  ks_statistic <- suppressWarnings(ks.test(generated_data, reference_data)$statistic)
  return(as.double(ks_statistic))
}

inference_retrec_interval_2d <- function(N0, N1, d0, p, sample_file, cost_function) {
  reference_data = scan(sample_file, quiet = TRUE)
  
  mu_init = median(reference_data)/N1
  fmutant_init = 1
  
  
  pas_mu = 1
  pas_fmutant = 0.35
  
  cost <- 1
  
  x <- 0
  y <- 0
  z <- 0
  
  for (j in 1:1) {
    x <- seq(log10(mu_init) - pas_mu, log10(mu_init) + pas_mu, by = 2*pas_mu/10)
    y1 <- 10^x
    x2 <- seq(log10(fmutant_init) - pas_fmutant, log10(fmutant_init) + pas_fmutant, by = 2*pas_fmutant/10)
    y2 <- 10^x2
    
    matrice_parametres = expand.grid(y1, y2)
    
    cost_function_values  = c(1:length(matrice_parametres[, 1]))
    
    x <- c(1:length(matrice_parametres[, 1]))
    y <- c(1:length(matrice_parametres[, 1]))
    z <- c(1:length(matrice_parametres[, 1]))
    
    for (i in 1:length(matrice_parametres[, 1])) {
      cost_function_values[i] = cost_function(matrice_parametres[i, 1], reference_data, N0, N1, d0, matrice_parametres[i, 2], p)
      x[i] = matrice_parametres[i, 1]
      y[i] = matrice_parametres[i, 2]
      z[i] = cost_function_values[i]
    }
    
    indice_min = which.min(cost_function_values)
    
    mu_init = matrice_parametres[indice_min, 1]
    fmutant_init = matrice_parametres[indice_min, 2]
    
    cost <- min(cost_function_values)
    
    print(cost)
    print(c(mu_init, fmutant_init))
    pas_mu = max(pas_mu/2, 0.1)
    if ((pas_fmutant > 0.1) && cost < 0.4 ) {
      pas_fmutant = pas_fmutant/2
    }
  }
  
  return(list(mu_init, fmutant_init, x, y, z))
}

vector <- inference_retrec_interval_2d(10, 1e7, 0, 1, "datachallengeA", cost_function)

reference_data <- scan("challengeC1e7_0_1.5e-6_1.2_1_96.txt", quiet = TRUE)

mu_init <- 1.5e-6
fmutant_init <- 1

pas_mu = 0.5
pas_fmutant = 0.35

library(rgl)
library(akima)


v <- seq(log10(mu_init) - pas_mu, log10(mu_init) + pas_mu, by = 2*pas_mu/10)
x <- 10^v
x2 <- seq(log10(fmutant_init) - pas_fmutant, log10(fmutant_init) + pas_fmutant, by = 2*pas_fmutant/10)
y <- 10^x2

# grid <- expand.grid(x = x, y = y)
# grid$z <- with(grid, cost_function(x, reference_data, 10, 1e7, 0, y, 1))  

z = matrix(nrow = 11, ncol = 11)

# z <- outer(x1, y1, FUN = function(x, y) cost_function(x, reference_data, 10, 1e7, 0, y, 1))

for (i in 1:length(x)) {
  for (j in 1:length(y)) {
    z[j, i] = cost_function(x[i], reference_data, 10, 1e7, 0, y[j], 1)
  }
}



# # Convert data to a regular grid
# xyz <- interp2xyz(x, y, z, maxpixels = 50000)
# 
# # Interpolate data on the regular grid
# grid <- with(xyz, interp(x, y, z, xo = x, yo = y))

# plot3d(x, y, z, type = "n", xlab = "mu", ylab = "fitness_mutant", zlab = "ditance_to_curve")
# 
# surface3d(x, y, z, alpha = 0.6, color = "blue")

library(ggplot2)
library(reshape2)

# Reshape the matrix to a data frame
df <- reshape2::melt(z, varnames = c("x", "y"), value.name = "z")

# Use ggplot2 to create a heatmap
ggplot(df, aes(x, y, fill = z)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red", limits = range(df$z),
                      breaks = seq(0, 1, by = 0.2), name = "Distance to Curve") +
  labs(x = "mu", y = "fitness_mutant") +
  theme_bw()

