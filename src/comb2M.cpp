#include <Rcpp.h>
using namespace Rcpp;

//' C++ version of combn function
//'
//' @param x Numeric vector
//' @return Numeric matrix with three rows: all pairwise combinations of values
//'         in vector x (rows 1 and 2), and their arithmetic difference (row 3)
//' @export
// [[Rcpp::export]]
NumericMatrix comb2M(NumericVector x) {
    int xlen = x.size();
    int counter = 0;
    double outsize = Rf_choose(xlen,2);
    NumericMatrix out(3,outsize);
    for (int a = 0; a < xlen; ++a) {
        for (int b = a+1; b < xlen; ++b) {
            out(0,counter) = x[a];
            out(1,counter) = x[b];
            if (x[b]>x[a]) {
                out(2,counter) = x[b] - x[a];
            } else {
                out(2,counter) = x[a] - x[b];
            }
            ++counter;
        }
    }
    return out;
}