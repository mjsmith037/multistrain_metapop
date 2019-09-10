library(igraph)
library(magrittr)
library(tidyverse)
library(tidygraph)
library(parallel)
library(JuliaCall)
library(kableExtra)
library(assertthat)
library(patchwork)
library(e1071)
library(ggraph)
library(janitor)

source("random_tree_networks.R")
source("directed_watts_strogatz_networks.R")

my_cols <- c("er"="#016394", "sb"="#f74700", "tr"="#b6003b", "ba"="#005342", "ws"="#1b264f",
             "Erdős-Rényi"="#016394", "Stochastic Block"="#f74700", "Tree"="#b6003b", "Barabasi-Albert"="#005342", "Watts-Strogatz"="#1b264f")

julia_setup()
julia_source("multipop_MANTIS.jl")

## structural parameters
set.seed(1)
struct <- c(2, 2)
strains <- expand.grid(lapply(struct, seq, from=1, by=1) %>%
                         setNames(str_c("locus", 1:length(.))))
n_populations <- 25

## dynamical parameters
R0 <- 4
sigma <- 40               # recovery rate
beta <- R0 * sigma                # Infection rate
mu <- 0.05                # mortality rate

movement_rate <- 0.05

n_sims <- 100

times <- Sys.time()

all_sims <- lapply(c("er", "sb", "ba", "tr", "ws"), function(net_type) {
  lapply(1:n_sims, function(i) {

    gamma <- runif(n_populations, 0.05, 0.95)    # partial cross-protective immunity (cpi)
    initial_conditions <- runif(prod(struct) * n_populations) %>%
      matrix(n_populations, prod(struct)) %>%
      apply(1, . %>% {. / (5 * sum(.))}) %>%
      t()

    if (net_type == "er") {
      movement_network <- play_erdos_renyi(n_populations, m=90)
    } else if (net_type == "sb") {
      movement_network <- play_blocks(n_populations, c(floor(n_populations / 2), ceiling(n_populations / 2)), matrix(c(0.28, 0.03, 0.03, 0.28), 2, 2))
    } else if (net_type == "ba") {
      movement_network <- play_barabasi_albert(n_populations, 0.25, 4) %>%
        activate(edges) %>% distinct(to, from) %>% filter(to != from)
    } else if (net_type == "tr") {
      movement_network <- simple_tree(n_populations, 1, 6, 12)
    } else if (net_type == "ws") {
      movement_network <- directed_watts_strogatz(n_populations, 4, 0.5, 90)
    }

    if (!with_graph(movement_network, graph_is_simple())) {
      if (not(exists("complex_graphs"))) complex_graphs <<- tibble()
      complex_graphs <<- bind_rows(complex_graphs, tibble(net_type, i))
      return(NULL)
    }


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

    ## run the simulation
    mantis_output <- julia_call("runMANTIS", strainstructure=struct, tstep=1, tmax=1000,
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
                                        "currently infectious (y)")))

    # assert_that(not(any(timeseries$prevalence > 1)))
    # assert_that(not(any(timeseries$prevalence < 0)))

    if (any(timeseries$prevalence > 1 | timeseries$prevalence < 0)) {
      if (not(exists("out_of_bounds"))) out_of_bounds <<- tibble()
      out_of_bounds <<- bind_rows(out_of_bounds,
                                 tibble(net_type, i,
                                        under=any(timeseries$prevalence < 0),
                                        over=any(timeseries$prevalence > 1)))
      return(NULL)
    }

    timeseries <- timeseries %>%
      filter(time > 0.75 * max(time)) %>%
      group_by(population, strain) %>%
      do(summarise(., ave_inf = sum(.$prevalence[.$equation == "y"]) / (max(time) - min(time)),
                   var_inf = var(.$prevalence[.$equation == "y"]),
                   ave_res = sum(.$prevalence[.$equation == "w"]) / (max(time) - min(time)),
                   var_res = var(.$prevalence[.$equation == "w"]),
                   ave_cycle = mutate(., prevalence = round(prevalence, 12)) %>% distinct(prevalence, .keep_all=TRUE) %>%
                     filter(equation == "y", lead(prevalence) > prevalence, lag(prevalence) > prevalence) %>%
                     use_series(time) %>% {ifelse(length(.) > 1, mean(diff(.)), NA)},
                   var_cycle = mutate(., prevalence = round(prevalence, 12)) %>% distinct(prevalence, .keep_all=TRUE) %>%
                     filter(equation == "y", lead(prevalence) > prevalence, lag(prevalence) > prevalence) %>%
                     use_series(time) %>% {ifelse(length(.) > 1, var(diff(.)), NA)},
                   num_nonzero = mutate(., prevalence = round(prevalence, 12)) %>%
                     filter(equation == "z", (lead(prevalence, default=2) > prevalence & lag(prevalence, default=2) > prevalence) |
                              (lead(prevalence, default=-1) < prevalence & lag(prevalence, default=-1) < prevalence)) %>%
                     summarise(num_nonzero_extrema = sum(prevalence > 0)) %>%
                     use_series(num_nonzero_extrema)))

    degree_distribution <- movement_network %>%
      activate(nodes) %>%
      transmute(indeg = centrality_degree(mode="in"),
                outdeg = centrality_degree(mode="out"),
                population=str_c("Population_", 1:n_populations)) %>%
      as_tibble()

    return(list(net_type, i, timeseries, degree_distribution, initial_conditions, gamma))
  })
})

if (exists("complex_graphs")) print(complex_graphs %>% tabyl(net_type))
if (exists("out_of_bounds")) print(out_of_bounds %>%
                                     mutate(both = under & over,
                                            under = ifelse(both, FALSE, under),
                                            over = ifelse(both, FALSE, over)) %>%
                                     select(-i) %>% group_by(net_type) %>% summarise_all(sum) %>%
                                     mutate(total = under + over + both))


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
                           {ggraph(., layout="manual", node.positions=layout_matrix) +
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


disease_stats <- lapply(all_sims, function(x) {
  lapply(x, function(y) {

    if (is.null(y)) return(NULL)

    net_type <- y[[1]]
    i <- y[[2]]

    output <- y[[3]]
    degree_distribution <- y[[4]]
    initial_conditions <- y[[5]]
    gamma <- y[[6]]

    variances <- var(degree_distribution)

    output  %>%
      ungroup() %>%
      select(-population, -strain) %>%
      summarise_all(mean, na.rm=TRUE) %>%
      mutate(connectance = mean(degree_distribution$indeg) / nrow(degree_distribution),
             in_deg_var = variances[1], out_deg_var = variances[4], in_out_cov = variances[2],
             in_skew = skewness(degree_distribution$indeg), out_skew = skewness(degree_distribution$outdeg),
             replicate = i,
             network = net_type)
  }) %>% bind_rows()
}) %>% bind_rows()

dynamics_stats <- lapply(all_sims, function(x) {
  lapply(x, function(y) {

    if (is.null(y)) return(NULL)

    net_type <- y[[1]]
    i <- y[[2]]

    output <- y[[3]]
    degree_distribution <- y[[4]]
    initial_conditions <- y[[5]]
    gamma <- y[[6]]

    variances <- var(degree_distribution)

    output %>%
      mutate(dynamics = case_when(num_nonzero == 0 ~ "extinct",
                                  num_nonzero == 1 ~ "steady",
                                  num_nonzero > 1  ~ "cycles")) %>%
      full_join(initial_conditions %>%
                  set_colnames(c("Strain_11", "Strain_12", "Strain_21", "Strain_22")) %>%
                  as_tibble() %>%
                  mutate(population = str_c("Population_", 1:25)) %>%
                  gather("strain", "init_y", -population), by = c("population", "strain")) %>%
      full_join(gamma %>% enframe(name=NULL, value="gamma") %>% mutate(population = str_c("Population_", 1:25)), by="population") %>%
      full_join(degree_distribution, by="population") %>%
      ungroup() %>%
      mutate(replicate = i,
             network = net_type)
  }) %>% bind_rows()
}) %>% bind_rows()

disease_stats %>%
  select(-replicate) %>%
  group_by(network) %>%
  summarise_all(. %>% {str_c(round(mean(.), 3), " (", round(sd(.), 3), ")")}) %>%
  rename("Total infected (y)" = tot_inf, "Mean specific immunity (z)" = mean_res) %>%
  kable(format="latex", booktabs=TRUE)

deg_data <- disease_stats %>%
  gather("metric", "value", connectance, in_deg_var, out_deg_var, in_out_cov, in_skew, out_skew) %>%
  mutate(network = factor(network, levels=rev(c("er", "sb", "tr", "ba", "ws")),
                          labels=rev(c("Erdős-Rényi", "Stochastic Block", "Tree", "Barabasi-Albert", "Watts-Strogatz"))),
         metric = factor(metric, levels=c("connectance", "in_deg_var", "out_deg_var", "in_out_cov", "in_skew", "out_skew"),
                         labels=c("Connectance", "Indegree variance", "Outdegree variance", "Degree co-variance", "Indegree skewness", "Outdegree skewness")))

plot_data <- disease_stats %>%
  gather("metric", "value", starts_with("tot"), starts_with("max"), starts_with("ave"), starts_with("var"), starts_with("mean")) %>%
  mutate(network = factor(network, levels=rev(c("er", "sb", "tr", "ba", "ws")),
                          labels=rev(c("Erdős-Rényi", "Stochastic Block", "Tree", "Barabasi-Albert", "Watts-Strogatz"))),
         metric = factor(metric, levels=c("ave_inf", "var_inf", "ave_res", "var_res", "ave_cycle", "var_cycle")))#,
                         # labels=c("Mean infected (y)", "Mean specific immunity (z)", "Mean time between epidemic peaks")))

degs <- ggplot(deg_data) +
  aes(x=network, y=value, colour=network) +
  geom_boxplot(outlier.alpha=0, size=0.75) +
  geom_jitter(width=0.25, alpha=0.15) +
  facet_wrap(~metric, scales="free") +
  scale_colour_manual(values=my_cols) +
  theme_bw() +
  theme(legend.position="none",
        axis.title=element_blank(),
        axis.text.x=element_text(angle=30, hjust=1))

met1 <- ggplot(plot_data %>% filter(metric == "Total infected")) +
  aes(x=network, y=value, colour=network) +
  geom_boxplot(outlier.alpha=0, size=0.75) +
  geom_jitter(width=0.25, alpha=0.15) +
  facet_wrap(~metric, scales="free_y") +
  scale_colour_manual(values=my_cols) +
  theme_bw() +
  theme(legend.position="none",
        axis.title=element_blank(),
        axis.text.x=element_text(angle=30, hjust=1))

met2 <- ggplot(plot_data %>% filter(metric == "Mean specific immunity (z)")) +
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

met3 <- ggplot(plot_data %>% filter(metric == "Mean time between epidemic peaks")) +
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
ggsave(str_c("large-network-comparison-sigma=", sigma, ".png"), width=7.5, height=6)

(met1 + met2 + met3) * coord_flip() + plot_layout(nrow=1)

# ggplot(timeseries %>% filter(time > 0.75*max(time), strain == "Strain_11")) +
#   aes(colour=population, y=prevalence, x=time) +
#   geom_line() +
#   # facet_grid(variable~strain, scales="free") +
#   facet_grid(variable~population, scales="free") +
#   scale_color_viridis_d() +
#   # scale_y_log10() +
#   theme_bw() +
#   theme(legend.position="bottom",
#         legend.background=element_blank(),
#         legend.title=element_blank())

times <- c(times, Sys.time())
diff(times)

disease_stats %>% select(ave_cycle, network) %>% tabyl(network)
disease_stats %>% select(ave_cycle, network) %>% na.omit() %>% tabyl(network)

