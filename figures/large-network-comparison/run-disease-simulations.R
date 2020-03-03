library(tidyverse)
library(magrittr)
library(JuliaCall)
library(patchwork)
library(ggraph)

################################################################################
#### PARAMETER/MODEL SET-UP: ADJUST THESE TO ALTER UNDERLYING DISEASE MODEL ####
################################################################################
## structural parameters
struct <- c(2, 2)
strains <- expand.grid(lapply(struct, seq, from=1, by=1) %>%
                         setNames(str_c("locus", 1:length(.))))
## dynamical parameters
gamma <- 0.66            # degree of cross-protective immunity
sigma <- 8               # recovery rate
mu <- 0.05               # mortality rate
movement_rate <- 0.01
################################################################################

my_cols <- c("er"="#016394", "sb"="#f74700", "tr"="#b6003b", "ba"="#005342", "ws"="#1b264f",
             "Erdős-Rényi"="#016394", "Stochastic Block"="#f74700", "Tree"="#b6003b", "Barabasi-Albert"="#005342", "Watts-Strogatz"="#1b264f")

count_unique_local_minima <- function(vector, digits=NA) {
  # -1 = extinct, 0 = stable, 1 = unconverged, >1 = chaos/cycles
  if (!is.na(digits)) vector <- round(vector, digits=digits)
  # if the last element is 0, the strain is extinct, regardless of other dynamics
  if (tail(vector, 1) == 0) return(-1)
  # pick out the local minima (use isless bc handles missing values introduced by lead/lag
  # silently: missing > all numbers)
  minima <- vector[vector < lag(vector) & vector < lead(vector)]
  number_unique_minima <- length(unique(minima))
  # if there are no minima, the strain is constant or has a period proportional to the stepsize (unlikely)
  # if only one minima, the vector is monotonic (implies failure to converge)
  if (number_unique_minima <= 1) return(number_unique_minima)
  # now check for monotonicity (read: failure to converge) in cyclical dynamics (i.e. length(minima) > 1)
  if (minima %>% subtract(lag(.)) %>% is_weakly_less_than(0) %>% all(na.rm=TRUE)) return(-1)
  if (minima %>% subtract(lag(.)) %>% is_weakly_greater_than(0) %>% all(na.rm=TRUE)) {
    maxima <- vector[vector > lag(vector) & vector > lead(vector)]
    if (maxima %>% subtract(lag(.)) %>% is_weakly_greater_than(0) %>% all(na.rm=TRUE)) {
      return(0)
    } else {
      return(1)
    }
  }
  # if not caught by any of the earlier conditions, then the dynamics are chaotic/cyclical
  return(number_unique_minima)
}

julia_setup()
julia_source("../../code/multipop_MANTIS.jl")

## load the networks
load("../../data/large-network-degree-distributions.RData")

run_all_sims <- function() {
  lapply(all_networks, function(one_type_of_network) {
    lapply(one_type_of_network, function(one_network) {

      set.seed(one_network$i)
      n_populations <- nrow(one_network$degree_distribution)

      initial_conditions <- runif(prod(struct) * n_populations) %>%
        matrix(n_populations, prod(struct)) %>%
        apply(1, . %>% {. / (5 * sum(.))}) %>%
        t()

      chi <- one_network$movement_network %>%
        activate(edges) %>%
        mutate(weight = movement_rate) %>%
        full_join(activate(., nodes) %>%
                    mutate(indeg = centrality_degree(weights=weight, mode="in"),
                           outdeg = centrality_degree(weights=weight, mode="out"),
                           n=1:n()) %>%
                    transmute(from=n, to=n, weight=min(indeg - outdeg, -indeg)) %>%
                    filter(weight != 0) %>%
                    as_tibble(), by=c("from", "to", "weight")) %>%
        as_tibble() %>%
        complete(from=1:n_populations, to=1:n_populations, fill=list("weight"=0)) %>%
        mutate(from = str_c("Population_", from),
               to   = str_c("Population_", to)) %>%
        spread(to, weight) %>%
        column_to_rownames("from") %>%
        as.matrix()

      R0 <- runif(n_populations, 1, 6)      # basic reproduction number
      beta <- R0 * (sigma + mu - diag(chi)) # infection rate

      ## run the simulation
      mantis_output <- julia_call("runMANTIS", strainstructure=struct, tstep=1000:1500, tmax=1500,
                                  beta=beta, gamma=gamma, sigma=sigma, chi=chi, mu=mu, initvals=initial_conditions)

      timeseries <- mantis_output$timeseries %>%
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
               prevalence = round(prevalence, 8)) %>%
        group_by(population, strain) %>%
        do(summarise(., ave_inf = sum(.$prevalence[.$equation == "y"]) / (max(time) - min(time)),
                     var_inf = var(.$prevalence[.$equation == "y"]),
                     ave_res = sum(.$prevalence[.$equation == "w"]) / (max(time) - min(time)),
                     var_res = var(.$prevalence[.$equation == "w"]),
                     ave_cycle = filter(., equation == "y", lead(prevalence) > prevalence, lag(prevalence) > prevalence) %>%
                       use_series(time) %>% {ifelse(length(.) > 1, mean(diff(.), na.rm=TRUE), NA)},
                     var_cycle = distinct(., prevalence, .keep_all=TRUE) %>%
                       filter(equation == "y", lead(prevalence) > prevalence, lag(prevalence) > prevalence) %>%
                       use_series(time) %>% {ifelse(length(.) > 1, var(diff(.), na.rm=TRUE), NA)},
                     num_nonzero = filter(., equation == "z") %>%
                       summarise(num_nonzero_extrema = count_unique_local_minima(prevalence)) %>%
                       use_series(num_nonzero_extrema)))

      return(list("network"=one_network[1:2], "timeseries"=timeseries,
                  "initial_conditions"=initial_conditions,
                  parameters=list("gamma"=gamma, "sigma"=sigma, "beta"=beta, "mu"=mu, "delta"=movement_rate)))
    })
  })
}

## Running simulations takes ~90 minutes to run on 8 cores
start_time <- Sys.time()
data_file <- str_glue("../../data/large-network-simulations_g{gamma},s{sigma},m{mu},d{movement_rate}.RData")
if (file.exists(data_file)) {
  load(data_file)
} else {
  all_sims <- run_all_sims()
  save(all_sims, file=data_file)
}
print(diff(c(start_time, Sys.time())))

layout_matrix <- rbind(example_er %>% as.igraph() %>% layout_nicely(),
                       example_sb %>% as.igraph() %>% layout_nicely(),
                       example_ba %>% as.igraph() %>% layout_nicely(),
                       example_ws %>% as.igraph() %>% layout_nicely(),
                       example_tr %>% as.igraph() %>% layout_as_tree()) %>%
  set_colnames(letters[24:25]) %>%
  as_tibble()

nets <- bind_graphs(example_er, example_sb, example_ba, example_ws, example_tr) %>%
  mutate(net_type = factor(net_type, levels=c("er", "sb", "tr", "ba", "ws"),
                           labels=c("Erdős-Rényi", "Stochastic Block", "Tree", "Barabasi-Albert", "Watts-Strogatz"))) %>%
  activate(nodes) %>%
  mutate(net_type = factor(net_type, levels=c("er", "sb", "tr", "ba", "ws"),
                           labels=c("Erdős-Rényi", "Stochastic Block", "Tree", "Barabasi-Albert", "Watts-Strogatz"))) %>%
                           {ggraph(., layout=layout_matrix) +
                               geom_edge_link(aes(colour=net_type),
                                              end_cap=circle(5, "pt"),
                                              edge_width=0.33,
                                              edge_alpha=0.5,
                                              arrow=arrow(angle=30, length=unit(3, "pt"), type="closed")) +
                               geom_node_point(aes(colour=net_type), size=1) +
                               facet_nodes(~net_type, scales="free", nrow=1) +
                               scale_colour_manual(values=my_cols) +
                               scale_edge_colour_manual(values=my_cols) +
                               theme_bw() +
                               theme(axis.text=element_blank(),
                                     axis.title=element_blank(),
                                     axis.ticks=element_blank(),
                                     panel.grid=element_blank(),
                                     legend.position="none")}

## plot some disease impact metrics
my_summary <- function(vec) {
  if (is.character(vec)) {
    return(vec %>% table() %>% {str_c(., names(.), collapse=";")})
  } else {
    return(vec %>% mean(na.rm=TRUE))
  }
}

plot_data <- lapply(all_sims, function(x) {
  lapply(x, function(y) {
    y$timeseries  %>%
      mutate(dynamics = case_when(num_nonzero == -1 ~ "extinct",
                                  num_nonzero ==  0 ~ "stable",
                                  num_nonzero ==  1 ~ "unconverged",
                                  TRUE ~ "cycles")) %>%
      ungroup() %>%
      select(-population, -strain) %>%
      summarise_all(my_summary) %>%
      mutate(replicate = y$network$i,
             network = y$network$net_type)
  }) %>% bind_rows()
}) %>% bind_rows()%>%
  gather("metric", "value", starts_with("tot"), starts_with("max"), starts_with("ave"), starts_with("var"), starts_with("mean")) %>%
  mutate(network = factor(network, levels=rev(c("er", "sb", "tr", "ba", "ws")),
                          labels=rev(c("Erdős-Rényi", "Stochastic Block", "Tree", "Barabasi-Albert", "Watts-Strogatz"))),
         metric = factor(metric, levels=c("ave_inf", "var_inf", "ave_res", "var_res", "ave_cycle", "var_cycle")))

met1 <- ggplot(plot_data %>% filter(metric == "ave_inf") %>%
                 mutate(metric = "Mean proportion of population infectious (y)")) +
  aes(x=network, y=value, colour=network) +
  geom_boxplot(outlier.alpha=0, size=0.75) +
  geom_jitter(width=0.25, alpha=0.15) +
  facet_wrap(~metric, scales="free_y") +
  scale_colour_manual(values=my_cols) +
  theme_bw() +
  theme(legend.position="none",
        axis.title=element_blank(),
        axis.text.x=element_text(angle=30, hjust=1))

met2 <- ggplot(plot_data %>% filter(metric == "ave_res") %>%
                 mutate(metric = "Mean cross-protective immunity (w)")) +
  aes(x=network, y=value, colour=network) +
  geom_boxplot(outlier.alpha=0, size=0.75) +
  geom_jitter(width=0.25, alpha=0.15) +
  facet_wrap(~metric, scales="free_y") +
  scale_colour_manual(values=my_cols) +
  theme_bw() +
  theme(legend.position="none",
        axis.title=element_blank(),
        axis.text.x=element_text(angle=30, hjust=1),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank())

met3 <- ggplot(plot_data %>% filter(metric == "ave_cycle") %>%
                 mutate(metric = "Mean time between epidemic peaks")) +
  aes(x=network, y=value, colour=network) +
  geom_boxplot(outlier.alpha=0, size=0.75) +
  geom_jitter(width=0.25, alpha=0.15) +
  facet_wrap(~metric, scales="free_y") +
  scale_colour_manual(values=my_cols) +
  scale_y_log10() +
  theme_bw() +
  theme(legend.position="none",
        axis.title=element_blank(),
        axis.text.x=element_text(angle=30, hjust=1),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank())

nets + ((met1 + met2 + met3) * coord_flip() + plot_layout(nrow=1)) + plot_layout(ncol=1, heights=1:2)
ggsave(str_replace(basename(data_file), ".RData$", ".png"), width=9, height=6)
