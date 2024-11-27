#' Find a set of nodes to cover a dataset
#' 
#' @param data A matrix of data points. Each row should be a data point.
#' @param n_nodes The number of nodes to cover the dataset.
#' @param max.it The number of iterations to use when finding the best set of nodes.
#' @param distance_type The type of distance to compute.
#'   - mahalanobis-scott takes into account the global covariate structure of the data with a scale adjustment for the dimension and number of data points.
#'   - mahalanobis is the same as mahalanobis-scott without the adjustment for dimension and number of data points.
#'   - independent-scott take into account the global variance of each dimension but not the correlation between coordinates with a scale adjustment for the dimension and number of data points.
#'   - independent is the same as independent-scott without the adjustment for dimension and number of data points.
#'   - euclidean is the standard euclidean distance.
#' @param jitter_magnitude The relative step size for randomly sampled jitters
#' 
#' @return A list giving the nodes and the objective function trajectory.
#' 
#' @export
adaptive_nodes<- function(
        data,
        n_nodes = 10,
        max.it = 100,
        distance_type = c(
            "mahalanobis-scott",
            "mahalanobis",
            "independent-scott",
            "independent",
            "euclidean"
        ),
        jitter_magnitude = 0.0001
    ) {
    if( !("matrix" %in% class(data)) ) data<- as.matrix(data)
    data<- data[complete.cases(data), , drop = FALSE]
    cov<- cov(data)
    precision<- solve(cov)
    if( distance_type[[1]] == "mahalanobis-scott" ) {
        data_precision<- precision * nrow(data) ^ (1 / (ncol(data) + 4))
        node_precision<- precision * n_nodes ^ (1 / (ncol(data) + 4))   
    } else if( distance_type[[1]] == "mahalanobis" ) {
        data_precision<- precision
        node_precision<- precision
    } else if( distance_type[[1]] == "independent-scott" ) {
        data_precision<- diag(diag(cov)^-1) * nrow(data) ^ (1 / (ncol(data) + 4))
        node_precision<- diag(diag(cov)^-1) * n_nodes ^ (1 / (ncol(data) + 4))
    } else if( distance_type[[1]] == "independent" ) {
        data_precision<- diag(diag(cov)^-1)
        node_precision<- data_precision
    } else if( distance_type[[1]] == "euclidean" ) {
        data_precision<- diag(ncol(data))
        node_precision<- diag(ncol(data))
    } else {
        stop("Distance type node recognized.")
    }
    node_cov<- solve(node_precision)

    penalty_history<- data.frame(
        penalty = numeric(max.it + 2),
        operation = ""
    )
    nodes<- data[sample(nrow(data), size = n_nodes), , drop = FALSE]
    penalty<- penalty_function(
        data,
        nodes,
        data_precision,
        node_precision
    )
    best_nodes<- nodes
    best_penalty<- penalty
    penalty_history$penalty[1]<- penalty
    penalty_history$operation[1]<- "Start"
    for( i in seq(max.it) ) {
        T<- temperature(i, max.it, 1, "decay")
        operation<- sample(
            c("Jitter", "Replace"),
            size = sample(2, size = 1, prob = c(0.8, 0.2)),
            prob = c(
                0.8,
                0.2
            )
        )
        new_nodes<- nodes
        if( "Jitter" %in% operation ) {
            new_nodes<- jitter(
                nodes,
                jitter_magnitude,
                node_cov
            )
        }
        if( "Replace" %in% operation ) {
            new_nodes<- replace(
                nodes,
                data,
                data_precision
            )
        }
        new_penalty<- penalty_function(
            data,
            new_nodes,
            data_precision,
            node_precision
        )
        if( accept(new_penalty, penalty, T, 0.1) ) {
            nodes<- new_nodes
            penalty<- new_penalty
        } else {
            operation<- "Keep"
        }
        if( penalty < best_penalty ) {
            best_nodes<- nodes
            best_penalty<- penalty
        } else if( restart(new_penalty, best_penalty, T, 0.3) ) {
            nodes<- best_nodes
            penalty<- best_penalty
            operation<- "Restart"
        }
        penalty_history$penalty[i + 1]<- penalty
        penalty_history$operation[i + 1]<- paste(operation, collapse = "")
    }
    nodes<- best_nodes
    penalty_history$penalty[i + 2]<- best_penalty
    penalty_history$operation[i + 2]<- "Best"
    return(list(nodes = best_nodes, penalty_history = penalty_history))
}

temperature<- function(
        i,
        max.it,
        T0 = 1,
        type = c(
            "default",
            "decay",
            "exp",
            "fast",
            "boltzmann"
        )
    ) {
    T<- switch(
        type[[1]],
        "default" = T0 * 0.5 * (1 - (i - 1) / max.it)^4,
        "decay" = T0 * (1 - (i - 1) / max.it),
        "exp" = T0 * 0.95^i,
        "fast" = T0 * (1 / i),
        "boltzmann" = T0 * (1 / log(i))
    )
    return(T)
}

accept<- function(
        new_penalty,
        old_penalty,
        T,
        max_probability
    ) {
    if( new_penalty < old_penalty ) return(TRUE)
    p<- exp( -(new_penalty - old_penalty) * (1 / T))
    p<- max_probability * p
    return(sample(c(TRUE, FALSE), 1, prob = c(p, 1 - p)))
}

restart<- function(
        new_penalty,
        best_penalty,
        T,
        max_probability
    ) {
    p<- 1 - exp( -(new_penalty - best_penalty) * T)
    p<- max_probability * p
    return(sample(c(TRUE, FALSE), 1, prob = c(p, 1 - p)))
}

jitter<- function(
        nodes,
        magnitude,
        covariance
    ) {
    nodes<- nodes + MASS::mvrnorm(
        n = nrow(nodes),
        mu = cbind(numeric(ncol(nodes))),
        Sigma = magnitude * covariance
    )
    return(nodes)
}

replace<- function(
        nodes,
        data,
        precision
    ) {
    d<- cross_distance(
        data,
        nodes,
        precision
    )
    node_min<- apply(
        d,
        MARGIN = 2,
        FUN = min
    )
    data_min<- apply(
        d,
        MARGIN = 1,
        FUN = min
    )
    delete_idx<- sample(
        nrow(nodes),
        size = 1,
        prob = node_min / sum(node_min)
    )
    add_idx<- sample(
        nrow(data),
        size = 1,
        prob = data_min / sum(data_min)
    )
    nodes[delete_idx, ]<- data[add_idx, ]
    return(nodes)
}