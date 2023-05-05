suppressPackageStartupMessages(
  {
    library(fitdistrplus)
    library(actuar)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(glue)
  }
)

# PDF = probability density function


# Here we load the static data file giving the list of distributions and
# associated parameters that have been fitted using LoS observations
# This data is required to inform parameter value choices for later simulation
# in this script
fitted_distributions <- 
  readr::read_csv(
    "/WholeSystMod/Syntax/WSSR/static file/data/statistics_distributions_Adult_2023-03-09.csv",
    show_col_types = F)

fitted_distributions.summary <- 
  fitted_distributions %>% 
  group_by(description) %>% 
  summarise(n=n(),
            p1.mean = mean(`parameter 1`),
            p2.mean = mean(`parameter 2`),
            p3.mean = mean(`parameter 3`),
            p1.sd = sd(`parameter 1`),
            p2.sd = sd(`parameter 2`),
            p3.sd = sd(`parameter 3`)) %>% arrange(-n)


fitted_distributions.summary[is.na(fitted_distributions.summary$description),"description"] <- "NOT FITTED"




# Function to fit a two parameter PDF to a univariate data set
# using maximum goodness-of-fit with Kolmogor-Smirnov stat as distance to be 
# minimised
fit.ks <- function(data,dist){
  r <- mgedist(data = data, distr = dist,gof = "KS")
  
  names(r$estimate) <- NULL
  rl <- list(ks = r$value,
             p1 = r$estimate[1],
             p2 = r$estimate[2])
  
  return(rl)
}
# Example
fit.ks(rgamma(1e3,5,2),"gamma")

# Run multiple repetitions (# = reps) of drawing n samples from a given PDF 
# (dist) and fitting via fit.ks()
simulate_fit <- function(n,reps,dist,p1,p2){
  df <- tibble(n=rep(n,reps),dist=dist,p1.true=p1,p2.true=p2)
  rngf <- paste0("r",dist) # name of random number generating function e.g.
  # rgamma
  df <- df %>% mutate(data=Map(rngf,n,p1,p2))
  df <- df %>% mutate(fit=Map(fit.ks,data,dist)) %>% unnest_wider(fit)
  return(df)
}

simulate_fit(50,3,"gamma",5,2) 

# Vectorised version of simulate_fit whereby n can be a vector
simulate_fit.multiple_n <- function(n,reps,dist,p1,p2){
  f <- function(x) simulate_fit(x,reps,dist,p1,p2)
  df <- bind_rows(lapply(n,f))
  return(df)
   
}

simulate_fit.multiple_n(50,3,"gamma",5,2)
simulate_fit.multiple_n(c(50,100),3,"gamma",5,2)

# get mean of fitted parameters
gamma_p1 <- fitted_distributions.summary %>% filter(description == "gamma") %>% pull(p1.mean)
gamma_p2 <- fitted_distributions.summary %>% filter(description == "gamma") %>% pull(p2.mean)
weibull_p1 <- fitted_distributions.summary %>% filter(description == "weibull") %>% pull(p1.mean)
weibull_p2 <- fitted_distributions.summary %>% filter(description == "weibull") %>% pull(p2.mean)
lnorm_p1 <- fitted_distributions.summary %>% filter(description == "lnorm") %>% pull(p1.mean)
lnorm_p2 <- fitted_distributions.summary %>% filter(description == "lnorm") %>% pull(p2.mean)

df_fitted_params <- 
  tibble(description = c("gamma","lnorm","weibull"),
         p1 = c(gamma_p1,lnorm_p1,weibull_p1),
         p2 = c(gamma_p2,lnorm_p2,weibull_p2))

# A plot showing the distribution of fitted parameters for a given
# PDF
fitted_distributions %>% 
  select(description,`parameter 1`,`parameter 2`) %>% 
  filter(description %in% c("gamma","lnorm","weibull")) %>% 
  ggplot(aes(`parameter 1`,`parameter 2`)) + 
  geom_point() + 
  geom_point(data = df_fitted_params, aes(p1,p2), colour = "red", size = 3) +
  facet_wrap(~description)

if (refresh_data <- F){
  N <- c(50,100,200,400,800,1600)
  Reps <- 1000
  df_sim_fit <- 
    bind_rows(
      simulate_fit.multiple_n("gamma", n=N ,reps = Reps,
                              p1 = gamma_p1, p2=gamma_p2),
      simulate_fit.multiple_n("weibull", n=N,reps = Reps,
                              p1 = weibull_p1, p2 = weibull_p2),
      simulate_fit.multiple_n("lnorm", n=N, reps = Reps, 
                              p1 = lnorm_p1, p2=lnorm_p2)
    )
  saveRDS(df_sim_fit,"data/simulated_fit.rds")  
} else {
  df_sim_fit <- readRDS("data/simulated_fit.rds")  
}


format(object.size(df_sim_fit),"MB")

if (F) {
  # playing around with z-scores
  x <- rnorm(1e3,50,5)
  x_scaled <- scale(x)[,1]
  hist(x)
  hist(x_scaled)
  y <- rnorm(1e3,50,25)
  y_scaled <- scale(y)[,1]
  hist(y)
  hist(y_scaled)
  sum( (x-mean(x))/sd(x) - scale(x)[,1] )  
}



