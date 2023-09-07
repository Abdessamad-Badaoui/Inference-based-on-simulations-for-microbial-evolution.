cost_function <- function(mu, reference_data, N0, N1, d0, fmutant, p) {
  commande <- sprintf("./atreyu_forward_simulator %e,%e %e %e %e %e 384 | cut -d' ' -f1 > TauxMutation2.txt", N0, N1, d0, mu, fmutant, p)
  system(commande)
  generated_data <- suppressMessages(scan("TauxMutation2.txt", quiet = TRUE))
  ks_statistic <- suppressWarnings(ks.test(generated_data, reference_data)$statistic)
  return(as.double(ks_statistic))
}


inference_retrec_interval <- function(N0, N1, d0, fmutant, p, sample_file, cost_function) {
  reference_data = suppressMessages(scan(sample_file, quiet = TRUE))
  
  mu_init = median(reference_data)/N1
  
  pas = 1
  
  for (i in 1:6) {
    x <- seq(log10(mu_init) - pas, log10(mu_init) + pas, by = 2*pas/10)
    y <- 10^x
    
    cost_function_values = c(1:length(y))
    
    for (i in 1:length(y)) {
      cost_function_values[i] = cost_function(y[i], reference_data, N0, N1, d0, fmutant, p)
    }
    
    mu_init = y[which.min(cost_function_values)]
    
    pas = pas/2
  }
  
  return(mu_init)
}

setwd(sprintf("%s%s", path.expand("~"),"/2ammis/Semestre2/ProjetSpe/pspeb2023"))

dataChallenge2 <- readLines("./challenges/challenge2.data")
fitnessChallenge2 <- readLines("./challenges/challenge2.f")
deathChallenge2 <- readLines("./challenges/challenge2.d")
samplingChallenge2 <- readLines("./challenges/challenge2.p")
populationChallenge2 <- readLines("./challenges/challenge2.N")

output2 <- "./challenges/challenge2_retrec_intervalle.mrate"

mrate_vector2 <- c(1:length(dataChallenge2))

if (TRUE) {    
  for (i in 1:length(dataChallenge2)) {
    data <- suppressMessages(scan(text = dataChallenge2[i], sep = ",", quiet = TRUE))
    sample_file <- writeLines(as.character(data), "temporary.txt")
    fmutant <- as.double(fitnessChallenge2[i])
    d0 <- as.double(deathChallenge2[i])
    p <- as.double(samplingChallenge2[i])
    N <- as.double(populationChallenge2[i])
    
    mu_infere <- inference_retrec_interval(10, N, d0, fmutant, p, "temporary.txt", cost_function)
    print(mu_infere)
    mrate_vector2[i] <- mu_infere 
  }
  
  writeLines(as.character(mrate_vector2), output2)
}
