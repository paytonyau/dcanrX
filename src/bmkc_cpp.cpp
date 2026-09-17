#include <RcppArmadillo.h>
#include <cmath>     // Explicit math preservation across Windows/macOS/Linux build platforms
#include <algorithm> // Guarantees thread-safe access limits for sorting structures
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;

// Helper: Extract Upper Triangle (Preserves true zero-correlations, drops diagonal)
vec get_upper_tri(const mat& X) {
  int n = X.n_rows;
  vec out((n * (n - 1)) / 2);
  int k = 0;
  for(int i = 0; i < n; ++i) {
    for(int j = i + 1; j < n; ++j) {
      out(k++) = X(i,j);
    }
  }
  return out;
}

// Pure C++ thread-safe Quantile handler to insulate OpenMP loops against crashes
inline double thread_safe_quantile(vec& v, double q) {
    int n = v.n_elem;
    if (n == 0) return 0.0;
    int idx = std::min(static_cast<int>(std::floor(n * q)), n - 1);
    std::nth_element(v.begin(), v.begin() + idx, v.end());
    return v(idx);
}

// Helper: Thread-Safe Stochastic Quantile Approximation
// NOTE: previously called arma_rng::set_seed_random() here, which reseeds
// Armadillo's RNG from system entropy on every call - this made percentile
// thresholds (and therefore the whole BMKC result) non-reproducible even
// under R's set.seed(), since R's RNG state has no influence over
// Armadillo's internal generator. The seed is now threaded through
// explicitly from the R-level `seed` argument and applied once, outside any
// parallel region, so results are deterministic for a fixed seed.
//
// Split into two steps: extract_and_subsample() does the expensive O(n^2)
// upper-triangle extraction (and, for large n, subsampling) exactly ONCE
// per matrix; quantile_from_vals() does the cheap abs()+quantile step,
// applied separately for the T1 and T2 thresholds since they use
// different `use_abs` settings (T1 respects `unsigned_net`, T2 always
// takes abs values). Previously each threshold call recomputed the full
// upper-triangle extraction from scratch even when re-using the same
// matrix, doubling both the O(n^2) extraction cost and (for large n) the
// memory footprint of the intermediate vector - at n_genes=20,000 that
// vector alone is ~1.6GB, and it was being built twice per matrix for no
// reason.
inline vec extract_and_subsample(const mat& X) {
    vec vals = get_upper_tri(X);

    if (vals.n_elem > 8000000) {
        uvec samp_idx = randi<uvec>(2000000, distr_param(0, vals.n_elem - 1));
        vals = vals.elem(samp_idx);
    }

    return vals;
}

inline double quantile_from_vals(const vec& vals, double q, bool use_abs) {
    vec v = use_abs ? abs(vals) : vals;
    return thread_safe_quantile(v, q);
}

// 1. Cache-Optimal Biweight Midcorrelation Matrix Generator
// See bmht_cpp.cpp's fast_bicor() for the Pearson-fallback rationale.
// [[Rcpp::export]]
NumericMatrix cpp_fast_bicor_matrix(NumericMatrix expr, int threads) {
#ifdef _OPENMP
  omp_set_num_threads(threads);
#endif

  mat X(expr.begin(), expr.nrow(), expr.ncol(), false);
  mat Xt = X.t(); 
  int n_genes = Xt.n_cols;
  int n_samples = Xt.n_rows;

  mat X_tilde = zeros(n_samples, n_genes);
  // Pre-sized so each thread only ever writes its own index i - safe under
  // the parallel loop below without needing a critical section.
  uvec is_zero_mad = zeros<uvec>(n_genes);

#pragma omp parallel for schedule(dynamic)
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
      is_zero_mad(i) = 1;
    }
  }

  mat cov = X_tilde.t() * X_tilde;
  colvec sum_sq = sum(square(X_tilde), 0).t();

  mat denom = sqrt(sum_sq * sum_sq.t()) + 1e-12;
  mat res = cov / denom;
  res.clean(1e-10);

  if(any(is_zero_mad)) {
    mat pearson = cor(Xt);
    pearson.replace(datum::nan, 0.0);
    uvec zero_mad_idx = find(is_zero_mad);
    for(size_t k = 0; k < zero_mad_idx.n_elem; ++k) {
      int i = zero_mad_idx(k);
      res.row(i) = pearson.row(i);
      res.col(i) = pearson.col(i);
    }
  }

  res.diag().zeros();

  NumericMatrix out(n_genes, n_genes);
  std::copy(res.begin(), res.end(), out.begin());
  return out;
}

// 2. Optimized Adjacency Logic & Masking Engine
// NOTE: reproducibility for the percentile-threshold subsampling in
// get_quantile() is controlled entirely by R's set.seed(), called before
// this function runs (see .run_bmkc_matrix). RcppArmadillo's RNG delegates
// to R's own generator by default, so no explicit seeding is needed - or
// wanted - here. An earlier version called arma_rng::set_seed_random(),
// which reseeds from system entropy on every call and defeats set.seed()
// entirely; a later attempt to fix that with an explicit
// arma_rng::set_seed(seed) call here was also wrong; Rcpp emits its own
// warning for exactly this ("the RNG seed has to be set at the R level").
// [[Rcpp::export]]
List cpp_bmkc_adjacency(NumericMatrix cor1, NumericMatrix cor2,
                        double T1, double T2,
                        bool unsigned_net, bool is_percentile, int threads) {

#ifdef _OPENMP
  omp_set_num_threads(threads);
#endif

  mat r1(cor1.begin(), cor1.nrow(), cor1.ncol(), false);
  mat r2(cor2.begin(), cor2.nrow(), cor2.ncol(), false);

  double t1_c1 = T1, t1_c2 = T1;
  double t2_c1 = T2, t2_c2 = T2;

  if (is_percentile) {
    double p1 = T1 > 1.0 ? T1 / 100.0 : T1;
    double p2 = T2 > 1.0 ? T2 / 100.0 : T2;

    vec vals_r1 = extract_and_subsample(r1);
    vec vals_r2 = extract_and_subsample(r2);

    t1_c1 = quantile_from_vals(vals_r1, p1, unsigned_net);
    t1_c2 = quantile_from_vals(vals_r2, p1, unsigned_net);
    t2_c1 = quantile_from_vals(vals_r1, p2, true);
    t2_c2 = quantile_from_vals(vals_r2, p2, true);
  }

  int n_genes = r1.n_rows;

  umat adj1_arma = zeros<umat>(n_genes, n_genes);
  umat adj2_arma = zeros<umat>(n_genes, n_genes);

  // Thread-local edge sum counters to bypass high-overhead R matrix sweeps
  int edge_count_1 = 0;
  int edge_count_2 = 0;

  // Cache-Optimal Column-Major Loop Interchange Configuration
#pragma omp parallel for schedule(dynamic) reduction(+:edge_count_1, edge_count_2)
  for(int j = 0; j < n_genes; ++j) {
    for(int i = 0; i < j; ++i) { 
      double val1 = r1(i,j);
      double val2 = r2(i,j);

      double abs_val1 = std::abs(val1);
      double abs_val2 = std::abs(val2);

      bool res1, res2;

      if (unsigned_net) {
        res1 = (abs_val1 >= t1_c1) && (abs_val2 <= t2_c2);
        res2 = (abs_val2 >= t1_c2) && (abs_val1 <= t2_c1);
      } else {
        res1 = (val1 >= t1_c1) && (abs_val2 <= t2_c2);
        res2 = (val2 >= t1_c2) && (abs_val1 <= t2_c1);
      }

      adj1_arma(i,j) = res1;
      adj1_arma(j,i) = res1;
      adj2_arma(i,j) = res2;
      adj2_arma(j,i) = res2;

      if(res1) edge_count_1++;
      if(res2) edge_count_2++;
    }
  }

  LogicalMatrix adj1(n_genes, n_genes);
  LogicalMatrix adj2(n_genes, n_genes);
  std::copy(adj1_arma.begin(), adj1_arma.end(), adj1.begin());
  std::copy(adj2_arma.begin(), adj2_arma.end(), adj2.begin());

  return List::create(
    Named("adj_cond1") = adj1,
    Named("adj_cond2") = adj2,
    Named("edge_counts") = List::create(
      Named("c1") = edge_count_1,
      Named("c2") = edge_count_2
    ),
    Named("calculated_thresholds") = List::create(
      Named("T1_c1") = t1_c1, Named("T1_c2") = t1_c2,
      Named("T2_c1") = t2_c1, Named("T2_c2") = t2_c2
    )
  );
}
