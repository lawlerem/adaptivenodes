#include <RcppEigen.h>

double distance_function(
        const Eigen::VectorXd &x,
        const Eigen::VectorXd &y,
        const Eigen::MatrixXd &precision
    ) {
    Eigen::VectorXd diff = x - y;
    Eigen::MatrixXd d2 = diff.transpose() * precision * diff;
    double d = sqrt(d2(0, 0));
    return d + 1e-9;
}

Eigen::MatrixXd self_distance(
        const Eigen::MatrixXd &x,
        const Eigen::MatrixXd &precision
    ) {
    Eigen::MatrixXd d(x.rows(), x.rows());
    for(int i = 0; i < x.rows(); i++) {
        d(i, i) = 0.0;
        for(int j = 0; j < i; j++) {
            d(i, j) = distance_function(
                x.row(i),
                x.row(j),
                precision
            );
            d(j, i) = d(i, j);
        }
    }
    return d;
}

Eigen::MatrixXd cross_distance(
        const Eigen::MatrixXd &x,
        const Eigen::MatrixXd &y,
        const Eigen::MatrixXd &precision
    ) {
    Eigen::MatrixXd d(x.rows(), y.rows());
    for(int i = 0; i < x.rows(); i++) {
        for(int j = 0; j < y.rows(); j++) {
            d(i, j) = distance_function(
                x.row(i),
                y.row(j),
                precision
            );
        }
    }
    return d;
}

Eigen::VectorXd vector_projection(
        const Eigen::VectorXd a,
        const Eigen::VectorXd b
    ) {
    return (a.dot(b) / b.dot(b)) * b;
}

double data_penalty(
        const Eigen::MatrixXd &data,
        const Eigen::MatrixXd &nodes,
        const Eigen::MatrixXd &precision
    ) {
    Eigen::MatrixXd d = cross_distance(
        data,
        nodes,
        precision
    );
    Eigen::VectorXd data_min = d.rowwise().minCoeff();
    Eigen::VectorXd node_min = d.colwise().minCoeff();
    Eigen::VectorXd weights(4); 
    weights << 0.4, 0.1, 0.4, 0.1;
    Eigen::VectorXd penalties(4);
    penalties << data_min.mean(), 
        data_min.maxCoeff(),
        node_min.mean(),
        node_min.maxCoeff();
    return vector_projection(penalties, weights).sum();
}

double node_penalty(
        const Eigen::MatrixXd &nodes,
        const Eigen::MatrixXd &precision
    ) {
    Eigen::MatrixXd d = self_distance(
        nodes,
        precision
    );
    Eigen::VectorXd weights(2);
    weights << 0.8, 0.2;
    Eigen::VectorXd inv_min_d(nodes.rows() - 1);
    for(int i = 0; i < inv_min_d.size(); i++) {
        inv_min_d(i) = 1.0 / d.col(i + 1).head(i + 1).minCoeff();
    }
    Eigen::VectorXd penalties(2);
    penalties << inv_min_d.mean(),
        inv_min_d.maxCoeff();
    return vector_projection(penalties, weights).sum();
}

// [[Rcpp::export("penalty_function")]]
double penalty_function(
        const Eigen::MatrixXd &data,
        const Eigen::MatrixXd &nodes,
        const Eigen::MatrixXd &data_precision,
        const Eigen::MatrixXd &node_precision
    ) {
    Eigen::VectorXd weights(2);
    weights << 0.5, 0.5;
    Eigen::VectorXd penalties(2);
    penalties << data_penalty(
            data,
            nodes,
            data_precision
        ), 
        node_penalty(
            nodes,
            node_precision
        );
    return vector_projection(penalties, weights).sum();
}



