library(magrittr)
library(tidyverse)
library(tidygraph)
library(JuliaCall)
library(ggraph)
library(facetscales)

julia_setup()
julia_source("../../Code/multipop_MANTIS.jl")

my_cols <- c("#f74700", "#016394", "#b6003b", "#005342")
scales_y <- list(
  `currently infectious (y)`    = scale_y_continuous(),
  `specific immunity (z)`       = scale_y_continuous(limits = c(0, 1)),
  `cross-reactive immunity (w)` = scale_y_continuous(limits = c(0, 1))
)

## structural parameters
set.seed(0)
struct <- c(2, 2)
strains <- expand.grid(lapply(struct, seq, from=1, by=1) %>%
                         setNames(str_c("locus", 1:length(.))))
n_populations <- 3
# initial_conditions <- rep(0.00001, prod(struct) * n_populations)
initial_conditions <- runif(prod(struct) * n_populations) %>%
  matrix(n_populations, prod(struct)) %>%
  apply(1, . %>% {. / (5 * sum(.))}) %>%
  t()
## dynamical parameters
beta <- 40                   # Infection rate
gamma <- c(0.25, 0.75, 0.25) # partial cross-protective immunity (cpi)
sigma <- 10                  # recovery rate
mu <- 0.05                   # disease induced mortality
# delta <- NA                  # increase in cpi per allele (not yet implemented)
# epsilon <- 0                 # seasonality

movement_rate <- 0.05
chi <- matrix(c(-movement_rate, 0, 0,
                0, -movement_rate, 0,
                movement_rate, movement_rate, -movement_rate), 3, 3)

timeseries <- julia_call("runMANTIS", strainstructure=struct, tstep=1, tmax=1000,
                         beta=beta, gamma=gamma, sigma=sigma, mu=mu, chi=chi,
                         initvals=initial_conditions)$timeseries %>%
  set_colnames(c(expand.grid(str_c("Population_", 1:n_populations),
                             str_c("Strain_", apply(strains, 1, str_c, collapse="")),
                             c("y", "z", "w")) %>% apply(1, str_c, collapse="."),
                 "time")) %>%
  gather("details", "prevalence", -time) %>%
  separate(details, c("population", "strain", "equation"), "\\.") %>%
  mutate(variable = factor(equation, levels=c("y", "z", "w"),
                           labels=c("currently infectious (y)",
                                    "specific immunity (z)",
                                    "cross-reactive immunity (w)")),
         population = factor(population, levels=c("Population_1", "Population_2", "Population_3"),
                             labels=c("Population A", "Population B", "Population C"))) %>%
  as_tibble()

ggplot(timeseries %>% filter(time > 0.85*max(time), strain == "Strain_11")) +
  aes(colour=population, y=prevalence, x=time) +
  geom_line() +
  facet_grid_sc(rows=vars(variable), cols=vars(population), scales=list(y=scales_y)) +
  scale_colour_manual(values=my_cols) +
  ylab("Proportion of population") +
  theme_bw() +
  theme(legend.position="bottom",
        legend.background=element_blank(),
        legend.title=element_blank())
# ggsave("dynamics-hierarchy.pdf", width=8, height=6)
ggsave("dynamics-hierarchy.png", width=8, height=6)
