library(magrittr)
library(tidyverse)
library(JuliaCall)

my_cols <- c("#f74700", "#016394", "#b6003b", "#005342")

# WARNING: this can take a long time to run
rnaught <- seq(1, 5, length.out=40)
mu <- seq(0.01, 0.25, length.out=40)
# note: other parameters not set here: strainstruct, gamma, sigma, mu, maxtime (see state_space_exploration.jl)

data_file <- str_glue("../../data/state-space-exploration_",
                      "rnaught={min(rnaught)}:{max(rnaught)}:{length(rnaught)}_",
                      "mu={min(mu)}:{max(mu)}:{length(mu)}.RData")

if (file.exists(data_file)) {
  load(data_file)
} else {
  julia_setup()
  julia_source("../../code/state_space_exploration.jl")
  simulation_output <- julia_call("state_space_exploration", rnaught=rnaught, mu=mu) %>%
    set_colnames(c("strain", "n_extrema", "gamma", "sigma", "R0", "mu")) %>%
    as_tibble()
  save(simulation_output, file=data_file)
}

plot_data <- simulation_output %>%
  mutate(dynamics = case_when(n_extrema == -1 ~ "extinct",
                              n_extrema == 0 ~ "steady",
                              n_extrema == 1 ~ "unconverged",
                              n_extrema > 1  ~ "cycles")) %>%
  group_by(gamma, sigma, mu, R0) %>%
  summarise(dynamics = dynamics %>% sort() %>% table() %>%
              {str_c(., names(.), collapse="")} %>%
  {case_when(. == "4steady" ~ "steady",
             . == "4cycles" ~ "cycles",
             str_detect(., "unconverged") ~ "unconverged",
             . == "2extinct2steady" ~ "extinct",
             TRUE ~ "unconverged")}) %>% ungroup() %>%
  select(gamma, sigma, mu, R0, dynamics) %>%
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
ggsave("state-space-exploration.png", width=10, height=10)
