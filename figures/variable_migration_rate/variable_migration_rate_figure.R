library(magrittr)
library(tidyverse)
library(JuliaCall)

julia_setup()
julia_source("state_space_exploration.jl")

# chi       <- c(5^-seq(0.5, 6, length.out=10), 10^-seq(0.5, 6, length.out=10))
# structure <- c(2, 2)
# maxtime   <- 10000
# gamma     <- c(0.25, 0.5, 0.67, 0.75)
# mu        <- 0.05
# sigma     <- c(1, 5, 10, 50)
# rnaught   <- c(2, 4, 6, 8)

chi       <- 10^-seq(0.5, 6, length.out=10)
structure <- c(2, 2)
maxtime   <- 10000
gamma     <- c(0.25, 0.75)
mu        <- 0.05
sigma     <- c(1, 5, 10, 50)
rnaught   <- c(2, 4, 6, 8)

mantis_output <- julia_call("variable_movement_rate", strainstruct=structure,
                            chi=chi, maxtime=maxtime, gamma=gamma, sigma=sigma, mu=mu, rnaught=rnaught)

save.image("~/Desktop/variable_migration_rate_round12.RData")
# load("~/Desktop/variable_migration_rate_round6.RData")

mantis_output %>%
  as.data.frame() %>%
  mutate_all(unlist) %>%
  as_tibble() %>%
  mutate(dynamics = case_when(V3 == 0 ~ "extinct",
                              V3 == 1 ~ "steady",
                              V3 > 1  ~ "cycles")) %>%
  rename(strain = V1, movement = V4, population = V2, gamma = V5, n_extrema = V3, sigma = V6, R0 = V7) %>%
  group_by(strain, movement, sigma, R0) %>%
  do(bind_cols(filter(., population == 1) %>%
                 transmute(gamma_1=gamma, n_extrema_1=n_extrema, dynamics_1=dynamics),
               filter(., population == 2) %>%
                 transmute(gamma_2=gamma, n_extrema_2=n_extrema, dynamics_2=dynamics))) %>%
  filter(gamma_1 == 0.25 & gamma_2 == 0.75 | gamma_1 == 0.75 & gamma_2 == 0.25) %>%
  group_by(movement, sigma, R0, gamma_1, gamma_2) %>%
  summarise(prop_transferred = sum(dynamics_1 == dynamics_2) / n()) %>%
  {ggplot(.) +
      aes(x=movement, y=prop_transferred, colour=str_c(gamma_1, " -> ", gamma_2)) +
      geom_point(size=2) + geom_line() +
      xlab("Instantaneous movement rate (proportion of source population)") +
      ylab("Proportion of cases where sink population has\ndynamics similar to those in source population") +
      scale_x_log10(breaks=10^-(6:1)) +
      # scale_shape_manual(values=c(4, 3)) +
      facet_grid(R0~sigma) +
      theme_bw() +
      theme(legend.title=element_blank(),
            legend.position=c(0.1,0.1),
            legend.background=element_rect(colour="black"))}
# ggsave("../Figures/variable_migration_rate.pdf", width=10, height=5.5)

# mantis_output %>%
#   as.data.frame() %>%
#   mutate_all(unlist) %>%
#   as_tibble() %>%
#   mutate(dynamics = case_when(V3 == 0 ~ "extinct",
#                               V3 == 1 ~ "steady",
#                               V3 > 1  ~ "cycles")) %>%
#   rename(strain = V1, movement = V4, population = V2, gamma = V5, n_extrema = V3, sigma = V6, R0 = V7) %>%
#   group_by(strain, movement, sigma, R0) %>%
#   do(bind_cols(filter(., population == 1) %>%
#                  transmute(gamma_1=gamma, n_extrema_1=n_extrema, dynamics_1=dynamics),
#                filter(., population == 2) %>%
#                  transmute(gamma_2=gamma, n_extrema_2=n_extrema, dynamics_2=dynamics))) %>%
#                  {ggplot(.) + aes(x=n_extrema_1, y=n_extrema_2) + geom_hex() + scale_fill_viridis_c(option="D")}
