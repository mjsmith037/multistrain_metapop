library(magrittr)
library(tidyverse)
library(JuliaCall)

julia_setup()
julia_source("state_space_exploration.jl")

movement_rates <- 10^-seq(1, 7, by=0.5)
structure <- c(2, 2)
maxtime   <- 10000
mu        <- 0.15

# list of parameter combinations (gamma, sigma, R0_1, R0_2)
parameter_sets <- [(0.55, 32, 2, 5),
                   (0.55, 32, 5, 2),
                   (0.66, 8, 2, 5),
                   (0.66, 8, 5, 2),
                   (0.77, 4, 2, 5),
                   (0.77, 4, 5, 3)]

for (movement_rate in movement_rates) {
    for (parameter_set in parameter_sets) {

    }
}

## what dynamics are expected?
expected_dynamics <- all_results %>%
  select(-R0_2, -dynamics_2, -rep) %>%
  group_by(gamma, sigma, R0_1, rounding_regime) %>%
  filter(movement == min(movement)) %>%
  select(-movement) %>%
  distinct() %>%
  rename(expected_1 = dynamics_1)

plot_data <- all_results %>%
  left_join(expected_dynamics,
            by=c("gamma", "sigma", "R0_1", "rounding_regime")) %>%
  filter_at(vars(starts_with("expected"), starts_with("dynamics")),
            all_vars(is_in(., c("cycles", "steady")))) %>%
  left_join(expected_dynamics %>% rename(R0_2 = R0_1, expected_2 = expected_1),
            by=c("gamma", "sigma", "R0_2", "rounding_regime")) %>%
  filter(expected_1 != expected_2) %>%
  mutate(trt = str_glue("{expected_1} %->% {expected_2}"))

gammas <- levels(factor(plot_data$gamma))
sigmas <- levels(fct_inorder(factor(plot_data %>% arrange(gamma) %>% use_series(sigma))))
rnaughts <- plot_data %>% group_by(gamma, sigma) %>% distinct(R0_1, R0_2) %>%
  filter(R0_1 < R0_2) %>% unite("x", R0_1, R0_2, sep=", ") %>% use_series(x)

labs <- list(bquote(gamma==.(gammas[1])*.(";")~sigma==.(sigmas[1])*.(";")~tilde(R[0])==.("[")*.(rnaughts[1])*.("]")),
             bquote(gamma==.(gammas[2])*.(";")~sigma==.(sigmas[2])*.(";")~tilde(R[0])==.("[")*.(rnaughts[2])*.("]")),
             bquote(gamma==.(gammas[3])*.(";")~sigma==.(sigmas[3])*.(";")~tilde(R[0])==.("[")*.(rnaughts[3])*.("]")))

ggplot(plot_data %>% filter(rounding_regime == "round8")) +
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
