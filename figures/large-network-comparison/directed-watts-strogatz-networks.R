## custom sample function to avoid bug where x is length 1 and is silently converted to 1:x
my_sample <- function (x, ...) x[sample.int(length(x), ...)]

# modified from: https://www.r-bloggers.com/generating-directed-watts-strogatz-network/
directed_watts_strogatz <- function(size, neighborhood_size, rewire_prob, desired_number_links) {
  stopifnot(neighborhood_size < size)
  edge.list <- vector("list", size)
  for (v in 0:(size-1)) {
    edge.end <- union((v + 1:neighborhood_size) %% size, (v + (-1:-neighborhood_size)) %% size)
    rewire <- (runif(length(edge.end)) < rewire_prob)
    edge.end <- edge.end[!rewire]
    rewired <- my_sample(setdiff(0:(size-1), c(edge.end, v)), sum(rewire))
    edges <- rep(v, 4 * neighborhood_size)
    edges[c(FALSE, TRUE)] <- c(edge.end, rewired)
    edge.list[[v + 1]] <- edges
  }
  graph_from_adj_list(edge.list %>% lapply(add, e2=1)) %>%
    # remove self-loops
    as_tbl_graph() %>% activate(edges) %>% filter(from != to) %>%
    # now subsample (randomly) to desired connectance
    sample_n(desired_number_links, .env=NULL) %>% activate(nodes)
}
