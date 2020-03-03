library(magrittr)
library(tidyverse)
library(JuliaCall)

my_cols <- c("#f74700", "#016394", "#b6003b", "#005342")

julia_setup()
julia_source("../../code/state_space_exploration.jl")

mantis_output <- julia_call("state_space_exploration",
                            rnaught=seq(1, 5, length.out=5),
                            m=seq(0.01, 0.25, length.out=5))
save.image("~/Desktop/tmp22_wR0_1pop.RData")

all_results <- list.files("../../Results/SSE1p_R0,mu/", full.names=TRUE) %>%
  lapply(read_csv,
         col_names=c("strain", "n_extrema_noround", "n_extrema_round8", "gamma", "sigma", "R0", "mu"),
         col_types="ddddddd") %>%
  bind_rows() %>%
  mutate_at(vars(starts_with("n_extrema")), get_dynamics) %>%
  group_by(gamma, sigma, mu, R0) %>%
  summarise_at(vars(starts_with("n_extrema")),
               . %>% sort() %>% table() %>% {str_c(., names(.), collapse="")}) %>%
  mutate_at(vars(starts_with("n_extrema")), ~case_when(. == "4steady" ~ "steady",
                                                       . == "4cycles" ~ "cycles",
                                                       str_detect(., "unconverged") ~ "unconverged",
                                                       . == "2extinct2steady" ~ "extinct",
                                                       TRUE ~ "unconverged"))

plot_data <- all_results %>% ungroup() %>%
  select(gamma, sigma, mu, R0, dynamics=n_extrema_round8) %>%
  complete(sigma, gamma, mu, R0, fill=list("dynamics"="unconverged"))

ggplot(plot_data) +
  aes(x=R0, y=mu, fill=dynamics) +
  xlab(expression(tilde(R)[0])) + ylab(expression(mu)) +
  geom_raster(alpha=0.7) +
  coord_cartesian(expand=FALSE) +
  scale_fill_manual(values=my_cols) +
  facet_grid(str_glue("{expression(sigma)}=={format(sigma, digits=3)}") ~
               str_glue("{expression(gamma)}=={round(gamma, 3)}"),
             labeller=label_parsed) +
  theme_bw() +
  theme(panel.grid=element_blank())
ggsave("state_space_[2,2]_R0xmu.png", width=10, height=10)
