#include <RcppArmadillo.h>
#include <cmath> // Explicitly preserves cross-compiler safety for std::sqrt across GCC, Clang, and MinGW
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;

// Hyper-fast cache-optimal Biweight Midcorrelation (bicor) engine
//
// Pearson fallback note: a gene with mad_x == 0 (>=50% of its values tied at
// the median) gets an all-zero X_tilde column below, which forces every
// correlation involving that gene to exactly 0 - not the same as R's
// bicor()/WGCNA::bicor(), both of which fall back to a standard Pearson
// correlation for that gene's row/column (pearsonFallback = "individual").
// The zero-MAD prefilter in the R engines only removes genes constant
// across ALL samples; a gene that is constant within just one condition
// (exactly the case a differential test should care most about) reaches
// this function with mad_x == 0 and was previously silently zeroed out.
// Fixed by computing a Pearson matrix once and splicing it in for exactly
// the affected rows/columns.
inline mat fast_bicor(const mat& X) {
    // Transpose once so each gene's samples are contiguous in memory.
    mat Xt = X.t(); 
    int n_genes = Xt.n_cols;
    int n_samples = Xt.n_rows;
    mat X_tilde = zeros(n_samples, n_genes);
    std::vector<int> zero_mad_genes;

    for(int i = 0; i < n_genes; ++i) {
        vec x = Xt.col(i); // Memory stream sequentializer allows direct CPU L1/L2 cache loads
        double med_x = median(x);
        vec abs_dev = abs(x - med_x);
        double mad_x = median(abs_dev); // raw MAD - do NOT scale by 1.4826; the "9" tuning constant already assumes raw MAD (matches WGCNA::bicor and Langfelder&Horvath 2012; verified to machine precision, see tests/testthat/test-bicor-parity.R)

        if(mad_x != 0.0) {
            vec u = (x - med_x) / (9.0 * mad_x);
            vec w = square(1.0 - square(u));
            w.elem(find(abs(u) >= 1.0)).zeros();
            X_tilde.col(i) = (x - med_x) % w; // SIMD auto-vectorization target
        } else {
            zero_mad_genes.push_back(i);
        }
    }

    mat cov = X_tilde.t() * X_tilde; 
    colvec sum_sq = sum(square(X_tilde), 0).t(); 
    
    mat denom = sqrt(sum_sq * sum_sq.t()) + 1e-12; 
    mat res = cov / denom;
    res.clean(1e-10); // Prunes computational rounding debris near zero bounds

    if(!zero_mad_genes.empty()) {
        // arma::cor() gives NaN for a constant column (0/0); R's bicor()
        // returns exactly 0 in that case (both-constant or one-constant),
        // so NaNs are zeroed after the fact rather than left to propagate.
        mat pearson = cor(Xt);
        pearson.replace(datum::nan, 0.0);
        for(size_t k = 0; k < zero_mad_genes.size(); ++k) {
            int i = zero_mad_genes[k];
            res.row(i) = pearson.row(i);
            res.col(i) = pearson.col(i);
        }
    }

    return res; 
}


// Workspace-based variant of fast_bicor(), used only inside the hot
// permutation loop below. A/B tested directly (see dev/ discussion):
// allocating fresh Xt/X_tilde/cov/denom/res every call costs essentially
// nothing at small-to-moderate gene counts (50-300 genes - noise-level
// difference) but grows to a real ~1.5-1.6x overhead at 600+ genes, which
// matters more now that adequate statistical power can require very large
// permutation counts. Reuses pre-sized buffers across permutations within
// a thread instead. Mathematically identical to fast_bicor() - verified
// via tests/testthat/test-bicor-parity.R style comparison before this was
// wired in (see test-workspace-parity.R).
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

inline void fast_bicor_ws(const mat& X, BicorWorkspace& ws) {
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
}

// Exported wrapper used to match baseline profiles with permutation pipelines
// [[Rcpp::export]]
arma::mat cpp_bmht_observed_bicor(NumericMatrix expr) {
    mat X(expr.begin(), expr.nrow(), expr.ncol(), false);
    mat res = fast_bicor(X);
    res.diag().zeros(); 
    return res;
}

// Thread-Safe, Invariant-Optimized Loop Interface
//
// CALIBRATION NOTE: earlier versions accepted a fixed `mask` computed once
// from the observed (true-label) correlation matrices and reused it for
// every permutation. That leaks the observed group structure into the null:
// the null distribution is only built from edges that were already selected
// for being strong in the real data, which shrinks its spread and makes the
// observed statistic look extreme almost regardless of its true value.
// Empirically this produced ~80% "significant" genes (FDR < 0.05) on pure
// noise data with no true signal. The fix is to recompute the
// informative-edge mask fresh from each permutation's own r1/r2 - the same
// selection rule the observed statistic uses, reapplied under the null.
// [[Rcpp::export]]
NumericMatrix cpp_bmht_permutations(NumericMatrix expr, IntegerMatrix shuffled_conds, double half_threshold, double mad_1, double mad_2, int threads) {
    mat X(expr.begin(), expr.nrow(), expr.ncol(), false);
    // NOTE: IntegerMatrix::begin()/LogicalMatrix::begin() return int*, but with
    // ARMA_64BIT_WORD (RcppArmadillo's default) arma::sword is a 64-bit type.
    // Aliasing an int* buffer as imat via the pointer constructor is a type
    // mismatch, not just a style choice - it either fails to compile or
    // reinterprets 32-bit ints as 64-bit and produces garbage. Use the
    // converting (copying) form instead.
    imat C = Rcpp::as<arma::imat>(shuffled_conds);

    int n_genes = X.n_rows;
    int n_perms = C.n_cols;

    NumericMatrix out_res(n_genes, n_perms);
    mat res(out_res.begin(), n_genes, n_perms, false);

    #ifdef _OPENMP
    omp_set_num_threads(threads);
    #endif

    // Group sizes are identical for every permutation (a permutation
    // reshuffles WHICH samples get which label, not how many of each), so
    // the workspace sizes below are valid for the whole loop, not just the
    // first iteration.
    ivec c_first = C.col(0);
    int n1 = as_scalar(sum(c_first == 1));
    int n2 = as_scalar(sum(c_first == 2));

    #pragma omp parallel
    {
        BicorWorkspace ws1(n1, n_genes);
        BicorWorkspace ws2(n2, n_genes);
        mat diff_sq(n_genes, n_genes);
        vec sum_vals(n_genes);
        ivec counts(n_genes);

        #pragma omp for schedule(dynamic)
        for(int p = 0; p < n_perms; ++p) {
            ivec c_shuffled = C.col(p);

            uvec idx1 = find(c_shuffled == 1);
            uvec idx2 = find(c_shuffled == 2);

            fast_bicor_ws(X.cols(idx1), ws1);
            fast_bicor_ws(X.cols(idx2), ws2);
            mat& r1 = ws1.res;
            mat& r2 = ws2.res;

            // Zero the diagonal (self-correlation) before scoring.
            r1.diag().zeros();
            r2.diag().zeros();

            // Recompute the informative-edge mask from THIS permutation's own
            // (unscaled) correlation matrices - mirrors the observed-data mask
            // construction exactly (abs(cor) > half_threshold in either
            // condition), but re-derived per permutation instead of fixed.
            umat mask_perm = (abs(r1) > half_threshold) || (abs(r2) > half_threshold);

            // Invariant Normalization executed in-place to avoid heap allocations
            r1 /= mad_1;
            r2 /= mad_2;

            diff_sq = square(r1 - r2);

            sum_vals.zeros();
            counts.zeros();

            // Column-Major Loop Interchange Transformation. Preserves absolute hardware locality.
            for(int j = 0; j < n_genes; ++j) {
                for(int i = 0; i < n_genes; ++i) {
                    if(mask_perm(i,j)) {
                        sum_vals(i) += diff_sq(i,j);
                        counts(i)++;
                    }
                }
            }

            for(int i = 0; i < n_genes; ++i) {
                res(i, p) = counts(i) > 0 ? std::sqrt(sum_vals(i) / counts(i)) : 0.0;
            }
        }
    }
    return out_res;
}
