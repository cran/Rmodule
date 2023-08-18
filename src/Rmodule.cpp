
#include <RcppArmadillo.h>

// [[Rcpp::export]]
Rcpp::NumericVector detR(Rcpp::NumericVector r, int param, const arma::mat& R, arma::mat R_)
{
    arma::uword n = R.n_rows;
    arma::uword m = r.length();
    Rcpp::NumericVector dets(m);
    arma::uword i, j;
    for (arma::uword k = 0; k < m; k++)
    {
        for (i = 0; i < n; i++)
            for (j = i + 1; j < n; j++)
                if (R(i, j) == param)
                {
                    R_(i, j) = r[k];
                    R_(j, i) = r[k];
                }
        dets[k] = det(R_);
    }
    return dets;
}

// [[Rcpp::export]]
double apesOfWrath(double rho, const arma::mat& intervals, double sigma)
{
    arma::uword loc = 0, i;
    arma::vec widths(intervals.n_rows);
    for (i = 0; i < intervals.n_rows; i++)
    {
        if (rho > intervals(i, 0) && rho < intervals(i, 1))
            loc = i;
        widths(i) = intervals(i, 1) - intervals(i, 0);
    }
    arma::rowvec aowInt = intervals.row(loc);
    i = 0;
    while (i < loc)
    {
        aowInt(0) = aowInt(0) - widths(i);
        i++;
    }
    i = loc + 1;
    while (i < intervals.n_rows)
    {
        aowInt(1) = aowInt(1) + widths(i);
        i++;
    }
    double a = aowInt(0);
    double b = aowInt(1);
    double etaOld = tan(M_PI / (b - a) * (rho - (a + b) / 2));
    Rcpp::NumericVector Z = Rcpp::rnorm(1, 0, sigma);
    double etaNew = etaOld + Z(0);
    double r = (b - a) / M_PI * atan(etaNew) + (a + b) / 2;
    if (r < intervals(loc, 0))
    {
        i = loc - 1;
        while (r < intervals(loc, 0) - widths(i))
            i--;
        r = intervals(i, 1) - (intervals(loc, 0) - r);
        if (intervals(i, 1) - r < 0.01)
            r = intervals(i, 1) - 0.01;
    }
    else if (r > intervals(loc, 1))
    {
        i = loc + 1;
        while (r > intervals(loc, 1) + widths(i))
            i++;
        r = intervals(i, 0) + (r - intervals(loc, 1));
        if (r - intervals(i, 0) < 0.01)
            r = intervals(i, 0) + 0.01;
    }
    else
    {
        if (r - intervals(loc, 0) < 0.01)
            r = intervals(loc, 0) + 0.01;
        else if (intervals(loc, 1) - r < 0.01)
            r = intervals(loc, 1) - 0.01;
    }
    return r;
}
