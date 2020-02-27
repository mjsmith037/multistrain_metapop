library(magrittr)
library(tidyverse)
library(tidygraph)
library(JuliaCall)
library(ggraph)
library(facetscales)

julia_setup()
julia_source("../../code/multipop_MANTIS.jl")

## plotting parameters
my_cols <- c("#f74700", "#016394", "#b6003b", "#005342")
scales_y <- list(
  `currently infectious (y)`    = scale_y_continuous(limits = c(0, NA)),
  `specific immunity (z)`       = scale_y_continuous(limits = c(0, 1)),
  `cross-reactive immunity (w)` = scale_y_continuous(limits = c(0, 1))
)

## structural parameters
set.seed(0)
struct <- c(2, 2)
strains <- expand.grid(lapply(struct, seq, from=1, by=1) %>%
                         setNames(str_c("locus", 1:length(.))))
n_populations <- 2
initial_conditions <- runif(prod(struct) * n_populations) %>%
  matrix(n_populations, prod(struct)) %>%
  apply(1, . %>% {. / (5 * sum(.))}) %>%
  t()
## dynamical parameters
gamma <- 0.66                   # partial cross-protective immunity (cpi)
sigma <- 8                      # recovery rate
mu <- 0.1                       # death rate
R0 <- c(2, 5)                   # (approximate) reproductive number
movement_rate <- 0.05

## integration parameters
maxtime <- 1000
tsteps <- round(0.85 * maxtime):maxtime

## run the simulation
timeseries <- bind_rows(
  julia_call("runMANTIS", strainstructure=struct, tstep=tsteps, tmax=maxtime,
             beta=R0 * (sigma + mu + movement_rate), gamma=gamma, sigma=sigma, mu=mu + movement_rate,
             chi=diag(0, 2), initvals=initial_conditions)$timeseries %>%
    set_colnames(c(expand.grid(str_c("Population_", 1:n_populations),
                               str_c("Strain_", apply(strains, 1, str_c, collapse="")),
                               c("y", "z", "w")) %>% apply(1, str_c, collapse="."),
                   "time")) %>%
    gather("details", "prevalence", -time) %>%
    separate(details, c("population", "strain", "equation"), "\\.") %>%
    mutate(variable = factor(equation, levels=c("z", "w", "y"),
                             labels=c("specific immunity (z)",
                                      "cross-reactive immunity (w)",
                                      "currently infectious (y)")),
           network = "A      B"),
  julia_call("runMANTIS", strainstructure=struct, tstep=tsteps, tmax=maxtime,
             beta=R0 * (sigma + mu + c(movement_rate, 0)), gamma=gamma, sigma=sigma, mu=mu,
             chi=matrix(c(-movement_rate, 0, movement_rate, 0), 2, 2),
             initvals=initial_conditions)$timeseries %>%
    set_colnames(c(expand.grid(str_c("Population_", 1:n_populations),
                               str_c("Strain_", apply(strains, 1, str_c, collapse="")),
                               c("y", "z", "w")) %>% apply(1, str_c, collapse="."),
                   "time")) %>%
    gather("details", "prevalence", -time) %>%
    separate(details, c("population", "strain", "equation"), "\\.") %>%
    mutate(variable = factor(equation, levels=c("z", "w", "y"),
                             labels=c("specific immunity (z)",
                                      "cross-reactive immunity (w)",
                                      "currently infectious (y)")),
           network = "A -> B"),
  julia_call("runMANTIS", strainstructure=struct, tstep=tsteps, tmax=maxtime,
             beta=R0 * (sigma + mu + c(0, movement_rate)), gamma=gamma, sigma=sigma, mu=mu,
             chi=matrix(c(0, movement_rate, 0, -movement_rate), 2, 2),
             initvals=initial_conditions)$timeseries %>%
    set_colnames(c(expand.grid(str_c("Population_", 1:n_populations),
                               str_c("Strain_", apply(strains, 1, str_c, collapse="")),
                               c("y", "z", "w")) %>% apply(1, str_c, collapse="."),
                   "time")) %>%
    gather("details", "prevalence", -time) %>%
    separate(details, c("population", "strain", "equation"), "\\.") %>%
    mutate(variable = factor(equation, levels=c("z", "w", "y"),
                             labels=c("specific immunity (z)",
                                      "cross-reactive immunity (w)",
                                      "currently infectious (y)")),
           network = "A <- B")
) %>%
  mutate(network = factor(network, levels=c("A      B", "A -> B", "A <- B")),
         variable = factor(variable, levels=c("currently infectious (y)",
                                              "specific immunity (z)",
                                              "cross-reactive immunity (w)")),
         population = factor(population, levels=c("Population_1", "Population_2"),
                             labels=c("Population A", "Population B"))) %>%
  as_tibble()

ggplot(timeseries %>% filter(strain == "Strain_11")) +
  aes(colour=population, y=prevalence, x=time) +
  geom_line(size=0.75) +
  facet_grid_sc(rows=vars(variable), cols=vars(network), scales=list(y=scales_y)) +
  scale_x_continuous(expand=c(0,0), breaks=c(860, 900, 940, 980)) +
  scale_colour_manual(values=my_cols) +
  ylab("Proportion of population") +
  theme_bw() +
  theme(legend.position="bottom",
        legend.background=element_blank(),
        legend.title=element_blank())
ggsave("two-pop-comparison.png", width=8, height=6)
