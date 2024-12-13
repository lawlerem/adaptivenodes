#include <RcppEigen.h>

double point_kde(
        const Eigen::VectorXd &x,
        const Eigen::MatrixXd &data,
        const Eigen::MatrixXd &precision
    ) {
    double pi = 3.141593;
    int k = data.cols();
    // data.rowwise() -= x.transpose();
    Eigen::VectorXd mh2(data.rows());
    for( int i = 0; i < data.rows(); i++ ) {
        mh2(i) = (data.row(i) - x) * precision *  (data.row(i).transpose() - x);
    }
    return exp(-0.5 * k * log(2 * pi) + 0.5 * log(precision.determinant()) + -0.5 * mh2.sum() / mh2.size());
}

// [[Rcpp::export("kde")]]
Eigen::VectorXd kde(
        const Eigen::MatrixXd &x,
        const Eigen::MatrixXd &data,
        const Eigen::MatrixXd &precision
    ) {
    Eigen::VectorXd ans(x.rows());
    for( int i = 0; i < x.rows(); i++ ) {
        ans(i) = point_kde(x.row(i), data, precision);
    }
    return ans;
}




