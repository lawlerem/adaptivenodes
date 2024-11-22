#' Find a set of nodes to cover a dataset
#' 
#' @param data A matrix of data points. Each row should be a data point.
#' @param n_nodes The number of nodes to cover the dataset.
#' @param max.it The number of iterations to use when finding the best set of nodes.
#' 
#' @return A list giving the nodes and the objective function trajectory.
#' 
#' @export
adaptive_nodes<- function(
        data,
        n_nodes = 10,
        max.it = 100
    ) {
    if( !("matrix" %in% class(data)) ) data<- as.matrix(data)
    cov<- cov(data)
    precision<- solve(cov)
    data_precision<- precision * nrow(data) ^ (1 / (ncol(data) + 4))
    node_precision<- precision * n_nodes ^ (1 / (ncol(data) + 4))
    node_cov<- solve(node_precision)

    penalty_history<- numeric(max.it + 1)
    nodes<- data[sample(nrow(data), size = n_nodes), , drop = FALSE]
    penalty<- penalty_function(
        data,
        nodes,
        data_precision,
        node_precision
    )
    penalty_history[1]<- penalty
    for( i in seq(max.it) ) {
        temperature<- 0.5 * ((max.it - i) / max.it)^4
        type<- sample(
            c("Jitter", "Replace"),
            size = sample(2, size = 1, prob = c(0.8, 0.2)),
            prob = c(
                0.8,
                0.2
            )
        )
        if( "Jitter" %in% type ) {
            new_nodes<- nodes + MASS::mvrnorm(
                n = nrow(nodes),
                mu = cbind(c(0, 0)),
                Sigma = 0.001 * node_cov
            )
        }
        if( "Replace" %in% type ) {
            new_nodes<- nodes
            new_nodes[
                sample(nrow(new_nodes), size = 1),

            ]<- data[
                sample(nrow(data), size = 1),

            ]
        }
        if( "Add" %in% type ) {
            new_nodes<- rbind(
                nodes,
                data[sample(nrow(data), size = 1), , drop = FALSE]
            )
        }
        new_penalty<- penalty_function(
            data,
            new_nodes,
            data_precision,
            node_precision
        )
        acceptance<- exp(-(new_penalty - penalty) / temperature)
        if( (new_penalty < penalty) | (stats::runif(1) < acceptance) ) {
            nodes<- new_nodes
            penalty<- new_penalty
        }
        penalty_history[i + 1]<- penalty
    }
    return(list(nodes = nodes, penalty_history = penalty_history))
}





