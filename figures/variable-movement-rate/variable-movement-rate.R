library(magrittr)
library(tidyverse)
library(JuliaCall)

# have to change the working directory for julia to find multipop_mantis.jl
setwd("../../code/")
julia_setup()
julia_source("state_space_exploration.jl")
setwd("../figures/variable-movement-rate/")

my_cols <- c("#f74700", "#016394", "#b6003b", "#005342")

movement_rates <- 10^-seq(1, 9, by=0.5)
# list of parameter combinations (gamma, sigma, R0_1, R0_2)
parameter_sets <- list(c(0.55, 32, 2, 5),
                       c(0.55, 32, 5, 2),
                       c(0.66,  8, 2, 5),
                       c(0.66,  8, 5, 2),
                       c(0.77,  4, 3, 5),
                       c(0.77,  4, 5, 3))
# note: other parameters not set here: strainstruct, maxtime, mu (see state_space_exploration.jl)

# WARNING: this can take a long time to run
data_file <- "../../data/variable-movement-rate.RData"
if (file.exists(data_file)) {
  load(data_file)
  old_sim_out <- simulation_output
  # check whether the loaded file already contains all of the parameters specified
  if (movement_rates %>% is_in(params$movement_rates) %>% all() %>% not() |
      parameter_sets %>% is_in(params$parameter_sets) %>% all() %>% not()) {
    # if it is only one of the two, we can easily run for the difference
    if (parameter_sets %>% is_in(params$parameter_sets) %>% all()) {
      movement_rates <- movement_rates[movement_rates %>% is_in(params$movement_rates) %>% not()]
    }
    if (movement_rates %>% is_in(params$movement_rates) %>% all()) {
      parameter_sets <- parameter_sets[parameter_sets %>% is_in(params$parameter_sets) %>% not()]
    }
    # else we just re-run everything
    simulation_output <- julia_call("variable_movement_rate",
                                    param_combinations=parameter_sets,
                                    chi=movement_rates)
    params <- list("movement_rates"=unique(c(movement_rates, params$movement_rates)),
                   "parameter_sets"=unique(c(parameter_sets, params$parameter_sets)))
    save(rbind(old_sim_out, simulation_output), params, file=data_file)
  }
} else {
  # if there was no data file, just proceed with the simulation
  simulation_output <- julia_call("variable_movement_rate",
                                  param_combinations=parameter_sets,
                                  chi=movement_rates)
  params <- list("movement_rates"=movement_rates, "parameter_sets"=parameter_sets)
  save(simulation_output, params, file=data_file)
}

cleaned_output <- simulation_output %>%
  set_colnames(c("strain", "population", "n_extrema", "movement", "gamma", "sigma", "R0")) %>%
  as_tibble() %>%
  mutate(dynamics = case_when(n_extrema == -1 ~ "extinct",
                              n_extrema == 0 ~ "steady",
                              n_extrema == 1 ~ "unconverged",
                              n_extrema > 1  ~ "cycles")) %>%
  group_by(strain, population, movement, gamma, sigma) %>% mutate(rep = 1:n()) %>% ungroup() %>%
  pivot_wider(names_from="population", values_from=c("n_extrema", "dynamics", "R0")) %>%
  select(-rep, -starts_with("n_extrema")) %>%
  group_by(movement, gamma, sigma, R0_1, R0_2) %>%
  summarise_at(vars(starts_with("dynamics")), . %>% sort() %>% table() %>% {str_c(., names(.), collapse="")}) %>%
  mutate_at(vars(starts_with("dynamics")), ~case_when(. == "4steady" ~ "steady",
                                                      . == "4cycles" ~ "cycles",
                                                      . == "2extinct2steady" ~ "extinct",
                                                      TRUE ~ "unconverged"))

## what dynamics are expected?
expected_dynamics <- cleaned_output %>%
  select(-R0_2, -dynamics_2) %>%
  group_by(gamma, sigma, R0_1) %>%
  filter(movement == min(movement)) %>%
  select(-movement) %>%
  distinct() %>%
  rename(expected_1 = dynamics_1)

plot_data <- cleaned_output %>%
  left_join(expected_dynamics, by=c("gamma", "sigma", "R0_1")) %>%
  filter_at(vars(starts_with("expected"), starts_with("dynamics")),
            all_vars(is_in(., c("cycles", "steady")))) %>%
  left_join(expected_dynamics %>% rename(R0_2 = R0_1, expected_2 = expected_1),
            by=c("gamma", "sigma", "R0_2")) %>%
  # filter(expected_1 != expected_2) %>%
  mutate(trt = str_c(expected_1, " %->% ", expected_2))

gammas <- levels(factor(plot_data$gamma))
sigmas <- levels(fct_inorder(factor(plot_data %>% arrange(gamma) %>% use_series(sigma))))
rnaughts <- plot_data %>% group_by(gamma, sigma) %>% distinct(R0_1, R0_2) %>%
  filter(R0_1 < R0_2) %>% unite("x", R0_1, R0_2, sep=", ") %>% use_series(x)

labs <- list(bquote(gamma==.(gammas[1])*.(";")~sigma==.(sigmas[1])*.(";")~tilde(R[0])==.("[")*.(rnaughts[1])*.("]")),
             bquote(gamma==.(gammas[2])*.(";")~sigma==.(sigmas[2])*.(";")~tilde(R[0])==.("[")*.(rnaughts[2])*.("]")),
             bquote(gamma==.(gammas[3])*.(";")~sigma==.(sigmas[3])*.(";")~tilde(R[0])==.("[")*.(rnaughts[3])*.("]")))

ggplot(plot_data) +
  aes(x=movement, y=as.integer(factor(dynamics_2)) - 1,
      colour=str_c(gamma, sigma), shape=str_c(gamma, sigma)) +
  geom_line(stat="smooth", method="glm", method.args=list(family="binomial"),
            se=FALSE, size=1.25, show.legend=FALSE, n=500, na.rm=TRUE) +
  geom_point(aes(size=str_c(gamma, sigma)), stroke=1.25, alpha=0.75) +
  xlab("Movement rate (proportion of source population per unit time)") +
  ylab("Destination disease dynamics") +
  scale_colour_manual(values=my_cols, labels=labs) +
  scale_shape_manual(values=c(1,3,4), labels=labs) +
  scale_size_manual(values=c(3,2.5,2.5), labels=labs) +
  scale_x_log10(breaks=c(1e-9, 1e-7, 1e-5, 0.001, 0.1),
                labels=c("1e-9", "1e-7", "1e-5", "0.001", "0.1")) +
  scale_y_continuous(breaks=0:1, labels=c("cycles", "steady")) +
  # facet_grid(rounding_regime~trt, labeller=label_parsed) +
  facet_wrap(.~trt, labeller=label_parsed, ncol=1) +
  theme_bw() +
  theme(legend.title=element_blank(),
        legend.position="bottom")
ggsave(filename="variable-movement-rate.png", width=7, height=7)
