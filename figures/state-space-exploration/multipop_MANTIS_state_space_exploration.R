library(magrittr)
library(tidyverse)
library(JuliaCall)

julia_setup()
julia_source("state_space_exploration.jl")

mantis_output <- julia_call("explore_parameter_space_one_pop")
save.image("~/Desktop/tmp22_wR0_1pop.RData")

clean_output <- . %>%
  as.data.frame() %>%
  mutate_all(unlist) %>%
  as_tibble() %>%
  rename(strain=V1, n_extrema=V2, gamma=V3, sigma=V4, R0=V5, mu=V6) %>%
  mutate(dynamics = case_when(n_extrema == 0 ~ "extinct",
                              n_extrema == 1 ~ "steady",
                              n_extrema > 1 ~ "cycles")) %>%
  group_by(gamma, sigma, R0, mu) %>%
  summarise(dynamics = dynamics %>% sort() %>% table() %>% {str_c(., names(.), collapse="")})

plot_data <- clean_output(mantis_output)

ggplot(plot_data) +
  aes(x=gamma, y=sigma, fill=dynamics) +
  geom_raster() +
  scale_y_log10() +
  scale_fill_viridis_d(option="C") +
  coord_cartesian(expand=FALSE) +
  facet_grid(mu~R0)
  

redone_params <- plot_data %>%
  filter(dynamics %>% is_in(c("4steady", "4cycles", "2extinct2steady")) %>% not()) %>%
  select(-dynamics) %>%
  t() %>%
  as.data.frame() %>% 
  as.list() %>%
  lapply(function(params) julia_call("explore_parameter_space_one_pop", tstep=0.5,
                                     gamma=params[1], sigma=params[2], rnaught=params[3], mu=params[4]))

new_plot_data <- lapply(redone_params, clean_output) %>% bind_rows()
plot_data <- bind_rows(plot_data, new_plot_data)

redone_params <- new_plot_data %>%
  filter(dynamics %>% is_in(c("4steady", "4cycles", "2extinct2steady")) %>% not()) %>%
  select(-dynamics) %>%
  distinct() %>%
  filter(R0 == 4, mu == 0.05) %>%
  t() %>%
  as.data.frame() %>% 
  as.list() %>%
  lapply(function(params) julia_call("explore_parameter_space_one_pop", tstep=0.05,
                                     gamma=params[1], sigma=params[2], rnaught=params[3], mu=params[4]))

new_plot_data <- lapply(redone_params, clean_output) %>% bind_rows()
plot_data <- bind_rows(plot_data, new_plot_data)


ggplot(plot_data %>%
         filter(dynamics %>% is_in(c("4steady", "4cycles", "2extinct2steady")))) +
  aes(x=gamma, y=sigma, fill=dynamics) +
  geom_raster() +
  scale_y_log10() +
  scale_fill_manual(values=c("#f74700", "#016394", "#b6003b", "#005342")) +
  coord_cartesian(expand=FALSE) +
  facet_grid(mu~R0)
