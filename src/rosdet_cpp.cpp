#include <RcppArmadillo.h>
#include <cmath>
#include <algorithm>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;

// Fast cache-optimal local bicor matrix calculator
// See bmht_cpp.cpp's fast_bicor() for the Pearson-fallback rationale.
inline mat fast_local_bicor(const mat& X) {
    mat Xt = X.t(); 
    int n_genes = Xt.n_cols;
    int n_samples = Xt.n_rows;
    mat X_tilde = zeros(n_samples, n_genes);
    std::vector<int> zero_mad_genes;

    for(int i = 0; i < n_genes; ++i) {
        vec x = Xt.col(i); 
        double med_x = median(x);
        vec abs_dev = abs(x - med_x);
        double mad_x = median(abs_dev); // raw MAD - do NOT scale by 1.4826; the "9" tuning constant already assumes raw MAD (matches WGCNA::bicor and Langfelder&Horvath 2012; verified to machine precision, see tests/testthat/test-bicor-parity.R)

        if(mad_x != 0.0) {
            vec u = (x - med_x) / (9.0 * mad_x);
            vec w = square(1.0 - square(u));
            w.elem(find(abs(u) >= 1.0)).zeros();
            X_tilde.col(i) = (x - med_x) % w;
        } else {
            zero_mad_genes.push_back(i);
        }
    }

    mat cov = X_tilde.t() * X_tilde; 
    colvec sum_sq = sum(square(X_tilde), 0).t(); 
    
    mat denom = sqrt(sum_sq * sum_sq.t()) + 1e-12; 
    mat res = cov / denom;
    res.clean(1e-10);

    if(!zero_mad_genes.empty()) {
        mat pearson = cor(Xt);
        pearson.replace(datum::nan, 0.0);
        for(size_t k = 0; k < zero_mad_genes.size(); ++k) {
            int i = zero_mad_genes[k];
            res.row(i) = pearson.row(i);
            res.col(i) = pearson.col(i);
        }
    }

    res.diag().zeros(); 

    return res; 
}

// Workspace-based variant of fast_local_bicor(), used only inside the hot
// bootstrap loop below. See bmht_cpp.cpp's BicorWorkspace/fast_bicor_ws for
// the rationale and A/B validation - identical pattern applied here.
struct BicorWorkspace {
    mat Xt, X_tilde, cov, denom, res;
    colvec sum_sq;
    std::vector<int> zero_mad_genes;

    BicorWorkspace(int n_samples, int n_genes) :
        Xt(n_samples, n_genes), X_tilde(n_samples, n_genes),
        cov(n_genes, n_genes), denom(n_genes, n_genes), res(n_genes, n_genes),
        sum_sq(n_genes) {
        zero_mad_genes.reserve(n_genes);
    }
};

inline void fast_local_bicor_ws(const mat& X, BicorWorkspace& ws) {
    ws.Xt = X.t();
    int n_genes = ws.Xt.n_cols;
    ws.X_tilde.zeros();
    ws.zero_mad_genes.clear();

    for(int i = 0; i < n_genes; ++i) {
        vec x = ws.Xt.col(i);
        double med_x = median(x);
        vec abs_dev = abs(x - med_x);
        double mad_x = median(abs_dev);

        if(mad_x != 0.0) {
            vec u = (x - med_x) / (9.0 * mad_x);
            vec w = square(1.0 - square(u));
            w.elem(find(abs(u) >= 1.0)).zeros();
            ws.X_tilde.col(i) = (x - med_x) % w;
        } else {
            ws.zero_mad_genes.push_back(i);
        }
    }

    ws.cov = ws.X_tilde.t() * ws.X_tilde;
    ws.sum_sq = sum(square(ws.X_tilde), 0).t();
    ws.denom = sqrt(ws.sum_sq * ws.sum_sq.t()) + 1e-12;
    ws.res = ws.cov / ws.denom;
    ws.res.clean(1e-10);

    if(!ws.zero_mad_genes.empty()) {
        mat pearson = cor(ws.Xt);
        pearson.replace(datum::nan, 0.0);
        for(size_t k = 0; k < ws.zero_mad_genes.size(); ++k) {
            int i = ws.zero_mad_genes[k];
            ws.res.row(i) = pearson.row(i);
            ws.res.col(i) = pearson.col(i);
        }
    }

    ws.res.diag().zeros();
}

// Exported wrapper to match baseline profiles with permutation pipelines
// [[Rcpp::export]]
arma::mat cpp_rosdet_observed_bicor(NumericMatrix expr) {
    mat X(expr.begin(), expr.nrow(), expr.ncol(), false);
    return fast_local_bicor(X);
}

// Thread-Safe, Zero-Contention Parallel Bootstrap Batch Processor (Gate 2 Engine)
// [[Rcpp::export]]
IntegerVector cpp_rosdet_bootstrap_batch(NumericMatrix expr_sub, IntegerVector g1_mapped, 
                                         IntegerVector g2_mapped, NumericVector obs_distance, 
                                         NumericVector c_weights, IntegerMatrix shuffled_labels, int threads) {
    
    mat X(expr_sub.begin(), expr_sub.nrow(), expr_sub.ncol(), false);
    uvec g1 = as<uvec>(g1_mapped);
    uvec g2 = as<uvec>(g2_mapped);
    vec obs_dist = as<vec>(obs_distance);
    vec weights = as<vec>(c_weights);
    // See bmht_cpp.cpp for why this must be a converting cast, not a pointer alias.
    imat S = Rcpp::as<arma::imat>(shuffled_labels);

    int n_pairs = g1.n_elem;
    int n_b_chunk = S.n_cols;

    int actual_threads = 1;
#ifdef _OPENMP
    actual_threads = threads;
    omp_set_num_threads(threads);
#endif

    // PERFORMANCE OPTIMIZATION: Thread-local reduction matrix avoids high-overhead atomic memory locks
    mat thread_counts = zeros<mat>(n_pairs, actual_threads);

    // Group sizes are identical for every permutation (see bmht_cpp.cpp for
    // the same reasoning), so workspace sizes below are valid for the
    // whole loop.
    ivec labels_first = S.col(0);
    int n1 = as_scalar(sum(labels_first == 1));
    int n2 = as_scalar(sum(labels_first == 2));

#pragma omp parallel
{
    BicorWorkspace ws1(n1, (int)X.n_rows);
    BicorWorkspace ws2(n2, (int)X.n_rows);
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif

    #pragma omp for schedule(dynamic)
    for(int b = 0; b < n_b_chunk; ++b) {
        ivec labels = S.col(b);

        uvec idx1 = find(labels == 1);
        uvec idx2 = find(labels == 2);

        fast_local_bicor_ws(X.cols(idx1), ws1);
        fast_local_bicor_ws(X.cols(idx2), ws2);
        mat& r1_null = ws1.res;
        mat& r2_null = ws2.res;

        // Threads write exclusively to their isolated column allocation spaces with zero contention
        for(int p = 0; p < n_pairs; ++p) {
            double null_r1 = r1_null(g1[p], g2[p]);
            double null_r2 = r2_null(g1[p], g2[p]);

            double null_dist = weights[p] * std::abs(null_r1 - null_r2);

            if(null_dist >= obs_dist[p]) {
                thread_counts(p, tid)++;
            }
        }
    }
}

    // Single-threaded aggregation sweep combines columns safely
    IntegerVector out_counts(n_pairs);
    for(int p = 0; p < n_pairs; ++p) {
        double sum_val = 0;
        for(int t = 0; t < actual_threads; ++t) {
            sum_val += thread_counts(p, t);
        }
        out_counts[p] = static_cast<int>(sum_val);
    }
    return out_counts;
}
