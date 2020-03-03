library(igraph)
library(tidygraph)
library(tidyverse)
library(magrittr)
library(e1071)

source("random-tree-networks.R")
source("directed-watts-strogatz-networks.R")

my_cols <- c("er"="#016394", "sb"="#f74700", "tr"="#b6003b", "ba"="#005342", "ws"="#1b264f",
             "Erdős-Rényi"="#016394", "Stochastic Block"="#f74700", "Tree"="#b6003b", "Barabasi-Albert"="#005342", "Watts-Strogatz"="#1b264f")

## structural parameters
n_populations <- 25
n_sims <- 100

all_networks <- lapply(c("er", "sb", "ba", "tr", "ws"), function(net_type) {
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

    return(list("net_type"=net_type, "i"=i, "movement_network"=movement_network, "degree_distribution"=degree_distribution))
  })
})


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

save(all_networks, example_er, example_sb, example_ba, example_tr, example_ws,
     degree_stats, file="../../data/large-network-degree-distributions.RData")

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
