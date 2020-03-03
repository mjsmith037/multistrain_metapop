library(igraph)
library(tidygraph)
library(tidyverse)
library(magrittr)
library(JuliaCall)
library(patchwork)
library(e1071)
library(ggraph)

source("random_tree_networks.R")
source("directed_watts_strogatz_networks.R")

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

## structural parameters
n_populations <- 25
n_sims <- 100
struct <- c(2, 2)
strains <- expand.grid(lapply(struct, seq, from=1, by=1) %>%
                         setNames(str_c("locus", 1:length(.))))

## dynamical parameters
gamma <- 0.66            # degree of cross-protective immunity
sigma <- 8               # recovery rate
mu <- 0.05               # mortality rate
movement_rate <- 0.01

run_all_sims <- run_all_sims <- function() {
  lapply(c("er", "sb", "ba", "tr", "ws"), function(net_type) {
    lapply(1:n_sims, function(i) {
      set.seed(i)
      if (net_type == "er") {
        movement_network <- play_erdos_renyi(n_populations, m=90)
      } else if (net_type == "sb") {
        movement_network <- play_blocks(n_populations, c(floor(n_populations / 2),
                                                         ceiling(n_populations / 2)),
                                        matrix(c(0.28, 0.03, 0.03, 0.28), 2, 2))
      } else if (net_type == "ba") {
        movement_network <- play_barabasi_albert(n_populations, 0.25, 4) %>%
          activate(edges) %>% distinct(to, from) %>% filter(to != from)
      } else if (net_type == "tr") {
        movement_network <- simple_tree(n_populations, 1, 6, 12)
      } else if (net_type == "ws") {
        movement_network <- directed_watts_strogatz(n_populations, 4, 0.5, 90)
      }

      degree_distribution <- movement_network %>%
        activate(nodes) %>%
        transmute(indeg = centrality_degree(mode="in"),
                  outdeg = centrality_degree(mode="out"),
                  population=str_c("Population_", 1:n_populations)) %>%
        as_tibble()

      initial_conditions <- runif(prod(struct) * n_populations) %>%
        matrix(n_populations, prod(struct)) %>%
        apply(1, . %>% {. / (5 * sum(.))}) %>%
        t()

      chi <- movement_network %>%
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

      degree_distribution <- movement_network %>%
        activate(nodes) %>%
        transmute(indeg = centrality_degree(mode="in"),
                  outdeg = centrality_degree(mode="out"),
                  population=str_c("Population_", 1:n_populations)) %>%
        as_tibble()

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

      return(list("net_type"=net_type, "rep"=i,
                  "movement_network"=movement_network, "degree_distribution"=degree_distribution,
                  "initial_conditions"=initial_conditions, "timeseries"=timeseries,
                  parameters=list("gamma"=gamma, "sigma"=sigma, "beta"=beta, "mu"=mu, "delta"=movement_rate)))
    })
  })
}

## Running simulations takes ~90 minutes to run on 8 cores
start_time <- Sys.time()
data_file <- str_glue("large-network-simulations_g{gamma},s{sigma},m{mu},d{movement_rate}.RData")
if (file.exists(data_file)) {
  load(data_file)
} else {
  all_sims <- run_all_sims()
  save(all_sims, file=data_file)
}
print(diff(c(start_time, Sys.time())))


## plot the degree distribution summary
degree_stats <- lapply(all_networks, function(x) {
  lapply(x, function(y) {
    variances <- var(y$degree_distribution %>% select(-population))
    tibble(connectance = mean(y$degree_distribution$indeg) / nrow(y$degree_distribution),
           in_deg_var = variances[1], out_deg_var = variances[4], in_out_cov = variances[2],
           in_skew = skewness(y$degree_distribution$indeg), out_skew = skewness(y$degree_distribution$outdeg),
           replicate = y$i,
           network = y$net_type)
  }) %>% bind_rows()
}) %>% bind_rows() %>%
  gather("metric", "value", connectance, in_deg_var, out_deg_var, in_out_cov, in_skew, out_skew) %>%
  mutate(network = factor(network, levels=rev(c("er", "sb", "tr", "ba", "ws")),
                          labels=rev(c("Erdős-Rényi", "Stochastic Block", "Tree", "Barabasi-Albert", "Watts-Strogatz"))),
         metric = factor(metric, levels=c("connectance", "in_deg_var", "out_deg_var", "in_out_cov", "in_skew", "out_skew"),
                         labels=c("Connectance", "Indegree variance", "Outdegree variance",
                                  "Degree co-variance", "Indegree skewness", "Outdegree skewness")))

degs <- ggplot(degree_stats) +
  aes(x=network, y=value, colour=network) +
  geom_boxplot(outlier.alpha=0, size=0.75) +
  geom_jitter(width=0.25, alpha=0.15) +
  facet_wrap(~metric, scales="free") +
  scale_colour_manual(values=my_cols) +
  theme_bw() +
  theme(legend.position="none",
        axis.title=element_blank(),
        axis.text.x=element_text(angle=30, hjust=1))
ggsave(degs, filename="large-network-degree-distributions.png", width=9, height=7)


## plot some example networks
set.seed(0)
example_er <- play_erdos_renyi(n_populations, m=90) %>% mutate(net_type = "er") %>%
  activate(edges) %>% mutate(net_type = "er")
example_sb <- play_blocks(n_populations, c(floor(n_populations / 2), ceiling(n_populations / 2)),
                          matrix(c(0.28, 0.03, 0.03, 0.28), 2, 2)) %>% mutate(net_type = "sb") %>%
  activate(edges) %>% mutate(net_type = "sb")
example_ba <- play_barabasi_albert(n_populations, 0.5, 4) %>% mutate(net_type = "ba") %>%
  activate(edges) %>% distinct(to, from) %>% filter(to != from) %>%
  activate(edges) %>% mutate(net_type = "ba")
example_tr <- simple_tree(n_populations, 1, 6, 12) %>% mutate(net_type = "tr") %>%
  activate(edges) %>% mutate(net_type = "tr")
example_ws <- directed_watts_strogatz(n_populations, 4, 0.5, 90) %>% mutate(net_type = "ws") %>%
  activate(edges) %>% mutate(net_type = "ws")

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
ggsave(str_replace(data_file, "RData$", ".png"), width=9, height=6)
