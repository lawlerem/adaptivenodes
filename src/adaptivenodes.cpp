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
    double penalty = data_min.sum() + node_min.sum() + 10 * node_min.maxCoeff();
    return penalty;
}

double node_penalty(
        const Eigen::MatrixXd &nodes,
        const Eigen::MatrixXd &precision
    ) {
    Eigen::MatrixXd d = self_distance(
        nodes,
        precision
    );
    Eigen::VectorXd min_d(nodes.rows() - 1);
    for(int i = 0; i < min_d.size(); i++) {
        min_d(i) = d.col(i + 1).head(i + 1).minCoeff();
    }
    double penalty = (1 / min_d.array()).sum();
    return penalty;
}

// [[Rcpp::export("penalty_function")]]
double penalty_function(
        const Eigen::MatrixXd &data,
        const Eigen::MatrixXd &nodes,
        const Eigen::MatrixXd &data_precision,
        const Eigen::MatrixXd &node_precision
    ) {
    double ans = 0.0;
    ans += data_penalty(
        data,
        nodes,
        data_precision
    );
    ans += node_penalty(
        nodes,
        node_precision
    );
    return ans;
}



