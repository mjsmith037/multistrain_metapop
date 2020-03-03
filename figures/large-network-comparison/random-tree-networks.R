# random_tree <- function(size, max_roots=1, children_min=4, children_max=6) {
#   nodes_left <- 1:size
#   edges <- tibble(from=NA, to=NA)
#   # select the initial parent(s)
#   parents <- sample(size, sample(max_roots, 1))
#   while (length(unique(unlist(edges))) != size) {
#     used_nodes <- c()
#     for (parent in parents) {
#       # choose between children_min and children_max children for each parent
#       nodes_to_choose_from <- nodes_left[(nodes_left != parent) &
#                                            (nodes_left %>%
#                                               is_in(edges %>% filter(to == parent) %>%
#                                                       use_series(from)) %>% not())]
#       children <- sample(nodes_to_choose_from,
#                          min(sample(children_min:children_max, 1),
#                              length(nodes_to_choose_from)))
#       # add the links:
#       edges <- bind_rows(edges, tibble(from = rep(parent, length(children)),
#                                        to   = children))
#       # keep track of used nodes (though allow multiple roots to connect to a
#       # other children in the same level
#       used_nodes <- c(used_nodes, children) %>% unique()
#     }
#     # remove the previous parents from the pool
#     nodes_left <- nodes_left[nodes_left %>% is_in(parents) %>% not()]
#     # set the children from this step to be the parents in the next
#     parents <- used_nodes %>% unique()
#   }
#   # tbl_graph(edges=edges %>% na.omit()) %>% plot()
#   return(tbl_graph(edges = edges %>% na.omit()))
# }

## custom sample function to avoid bug where x is length 1 and is silently converted to 1:x
my_sample <- function (x, ...) x[sample.int(length(x), ...)]

simple_tree <- function(size, max_roots=1, children_min=5, children_max=10) {
  nodes_left <- 1:size
  edges <- tibble()
  # select the initial parent(s)
  parents <- sample(size, sample(max_roots, 1))
  nodes_left <- nodes_left[nodes_left %>% is_in(parents) %>% not()]
  while (nodes_left %>% is_empty() %>% not()) {
    used_nodes <- c()
    for (parent in parents) {
      # choose between children_min and children_max children for each parent
      children <- my_sample(nodes_left, min(sample(children_min:children_max, 1), length(nodes_left)))
      # add the links:
      edges <- bind_rows(edges, tibble(from = rep(parent, length(children)),
                                       to   = children))
      # keep track of used nodes (though allow multiple parents to connect to a
      # other children in the same level
      used_nodes <- c(used_nodes, children) %>% unique()
    }
    # set the children from this step to be the parents in the next
    parents <- used_nodes %>% unique()
    # remove the previous parents from the pool
    nodes_left <- nodes_left[nodes_left %>% is_in(parents) %>% not()]
  }
  # tbl_graph(edges=edges %>% na.omit()) %>% plot()
  return(tbl_graph(edges = edges))
}
