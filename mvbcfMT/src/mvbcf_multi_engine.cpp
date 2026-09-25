// mvbcf_multi_engine.cpp
// -----------------------------------------------------------------------------
// Unified multivariate Bayesian Causal Forest engine for an ARBITRARY number of
// additive effect components. The response surface is written as
//
//     y_i = sum_{k=0}^{K-1}  D_i^{(k)} (.)  f_k( X_i^{(k)} )  +  eps_i ,
//               eps_i ~ N_q(0, Sigma)
//
// where each f_k is its own sum-of-trees forest with q-variate terminal nodes,
// D^{(k)} is an (n x q) 0/1 (or interaction-product) indicator matrix that
// "switches on" component k, and (.) is the elementwise (Hadamard) product.
//
//   * component k = 0 is the PROGNOSTIC forest  mu(.)  with  D^{(0)} == 1
//   * a main treatment effect t uses           D = Z_t (replicated over q cols)
//   * a pairwise interaction {s,t} uses         D = (Z_s * Z_t)
//   * an order-r interaction uses               D = prod_{t in S} Z_t
//
// This ONE engine therefore implements:
//   - the original single-treatment MVBCF                (K = 2: mu, tau)
//   - the two-treatment factorial MVBCF                  (K = 4: mu, tau1, tau2, tau12)
//   - the general T-treatment ANOVA/Mobius truncation    (K = 1 + sum_{j<=r} C(T,j))
//
// The node conditional posterior and marginal likelihood used here are the
// generic "tau" (indicator-weighted) versions from the original engine; when
// D == 1 they reduce EXACTLY to the prognostic-node updates (see theory doc,
// Prop. 1), so mu is just the D==1 special case and no code is duplicated.
//
// Reuses the Node / Tree / Forest classes of the original fast_bart engine.
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
#include <RcppDist.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppDist)]]

// ============================= Node ==========================================
class Node {
public:
  arma::colvec mu;             // q-variate terminal-node parameter (mu OR tau_k)
  int variable;
  double split_val;
  arma::uvec observations;
  arma::uvec test_observations;
  bool is_terminal;
  bool in_use;

  Node() { variable = -1; split_val = -1; is_terminal = false; in_use = false; }

  Node(const Node& o) {
    mu = o.mu; variable = o.variable; split_val = o.split_val;
    observations = o.observations; test_observations = o.test_observations;
    is_terminal = o.is_terminal; in_use = o.in_use;
  }

  // Generic indicator-weighted node update (Normal-Normal conjugate draw).
  // With D == 1 this reduces to the prognostic (mu) update.
  void update_node(arma::mat sigma, arma::mat sigma_par,
                   arma::mat y_resid, arma::mat D) {
    arma::colvec mu_0(y_resid.n_cols, arma::fill::zeros);
    arma::mat node_d = D.rows(find(observations == 1));
    arma::mat njd = node_d.t() * node_d;                       // q x q
    arma::mat node_resid = y_resid.rows(find(observations == 1));
    arma::mat part1 = arma::inv(arma::inv(sigma_par) + njd % arma::inv(sigma));
    arma::mat tricky(y_resid.n_cols, 1, arma::fill::zeros);
    arma::mat node_d_t = node_d.t();
    arma::mat node_resid_t = node_resid.t();
    for (unsigned int i = 0; i < node_resid.n_rows; i++)
      tricky = tricky + node_d_t.col(i) % (arma::inv(sigma) * node_resid_t.col(i));
    arma::mat part2 = arma::inv(sigma_par) * mu_0 + tricky;
    arma::rowvec temp = rmvnorm(1, part1 * part2, part1);
    mu = temp.t();
  }
};

// ============================= Tree ==========================================
class Tree {
public:
  std::vector<Node> node_vector;

  Tree(int num_nodes = 1, int num_obs = 1, int num_test_obs = 1) {
    node_vector.resize(num_nodes);
    node_vector[0].observations = arma::uvec(num_obs, arma::fill::ones);
    node_vector[0].test_observations = arma::uvec(num_test_obs, arma::fill::ones);
    node_vector[0].in_use = true;
    node_vector[0].is_terminal = true;
  }
  Tree(const Tree& other) {
    for (const Node& node : other.node_vector) node_vector.push_back(Node(node));
  }

  void update_nodes(arma::mat sigma, arma::mat sigma_par,
                    arma::mat y_resid, arma::mat D) {
    for (size_t i = 0; i < node_vector.size(); i++)
      if (node_vector[i].is_terminal && node_vector[i].in_use)
        node_vector[i].update_node(sigma, sigma_par, y_resid, D);
  }

  int get_terminal_node() {
    std::vector<int> v;
    for (size_t i = 0; i < node_vector.size(); i++)
      if (node_vector[i].is_terminal && node_vector[i].in_use) v.push_back(i);
    return v[floor(R::runif(0, v.size()))];
  }
  int get_non_terminal_node() {
    std::vector<int> v;
    for (size_t i = 0; i < node_vector.size(); i++)
      if (!node_vector[i].is_terminal && node_vector[i].in_use) v.push_back(i);
    return v.size() > 0 ? v[floor(R::runif(0, v.size()))] : -1;
  }
  int get_parent_child() {
    std::vector<int> v;
    for (size_t i = 0; i < node_vector.size(); i++)
      if (!node_vector[i].is_terminal && node_vector[i].in_use && (i != 0)) v.push_back(i);
    return v.size() > 0 ? v[floor(R::runif(0, v.size()))] : -1;
  }
  int get_terminal_parent() {
    std::vector<int> v;
    for (size_t i = 0; i < node_vector.size(); i++)
      if (node_vector[i].in_use && !node_vector[i].is_terminal)
        if (node_vector[2*i+1].in_use && node_vector[2*i+1].is_terminal &&
            node_vector[2*i+2].in_use && node_vector[2*i+2].is_terminal) v.push_back(i);
    return v.size() > 0 ? v[floor(R::runif(0, v.size()))] : -1;
  }

  void grow(arma::mat X, arma::mat X_test, int p, int min_nodesize) {
    int gi = get_terminal_node();
    int variable = floor(R::runif(0, p));
    node_vector[gi].variable = variable;
    arma::colvec Xc = X.col(variable), Xtc = X_test.col(variable);
    arma::colvec Xcs = Xc.rows(find(node_vector[gi].observations == 1));
    arma::colvec Xu = arma::unique(Xcs);
    double split_val;
    if (Xu.n_rows > 0) split_val = Xu(floor(R::runif(0, Xu.n_rows)));
    else split_val = -1;
    node_vector[gi].split_val = split_val;
    arma::uvec is_less = Xc <= split_val;
    arma::uvec is_less_test = Xtc <= split_val;
    arma::uvec less_subset = node_vector[gi].observations && is_less;
    arma::uvec more_subset = node_vector[gi].observations && (1 - is_less);
    int sum_less = sum(less_subset), sum_more = sum(more_subset);
    if (sum_more >= min_nodesize && sum_less >= min_nodesize) {
      if ((int)node_vector.size() < 2*gi+2+1) node_vector.resize(2*gi+2+1);
      int cl = 2*gi+1, cr = 2*gi+2;
      node_vector[cl].observations = node_vector[gi].observations && is_less;
      node_vector[cl].test_observations = node_vector[gi].test_observations && is_less_test;
      node_vector[cl].is_terminal = true; node_vector[cl].in_use = true;
      node_vector[cr].observations = node_vector[gi].observations && (1 - is_less);
      node_vector[cr].test_observations = node_vector[gi].test_observations && (1 - is_less_test);
      node_vector[cr].is_terminal = true; node_vector[cr].in_use = true;
      node_vector[gi].is_terminal = false; node_vector[gi].in_use = true;
    }
  }
  void prune() {
    int pi = get_terminal_parent();
    if (pi != -1) {
      node_vector[pi*2+1].in_use = false; node_vector[pi*2+1].is_terminal = false;
      node_vector[pi*2+2].in_use = false; node_vector[pi*2+2].is_terminal = false;
      node_vector[pi].is_terminal = true; node_vector[pi].in_use = true;
    }
  }
  void change_update(arma::mat X, arma::mat X_test) {
    for (size_t i = 0; i < node_vector.size(); i++)
      if (!node_vector[i].is_terminal && node_vector[i].in_use) {
        int cl = 2*i+1, cr = 2*i+2, variable = node_vector[i].variable;
        double sv = node_vector[i].split_val;
        arma::uvec il = X.col(variable) <= sv, im = X.col(variable) > sv;
        arma::uvec ilt = X_test.col(variable) <= sv, imt = X_test.col(variable) > sv;
        node_vector[cl].observations = node_vector[i].observations && il;
        node_vector[cr].observations = node_vector[i].observations && im;
        node_vector[cl].test_observations = node_vector[i].test_observations && ilt;
        node_vector[cr].test_observations = node_vector[i].test_observations && imt;
      }
  }
  void change(arma::mat X, int p) {
    int ci = get_non_terminal_node();
    if (ci != -1) {
      int variable = floor(R::runif(0, p));
      node_vector[ci].variable = variable;
      arma::colvec Xc = X.col(variable);
      Xc = Xc.rows(find(node_vector[ci].observations == 1));
      arma::colvec Xu = arma::unique(Xc);
      if (Xu.size() > 0) node_vector[ci].split_val = Xu(floor(R::runif(0, Xu.n_rows)));
      else node_vector[ci].split_val = -1;
    }
  }
  void swap() {
    int si = get_parent_child();
    if (si != -1) {
      int pi = (si-1)/2;
      int pv = node_vector[pi].variable; double ps = node_vector[pi].split_val;
      int cv = node_vector[si].variable; double cs = node_vector[si].split_val;
      node_vector[pi].variable = cv; node_vector[pi].split_val = cs;
      node_vector[si].variable = pv; node_vector[si].split_val = ps;
    }
  }
  bool has_empty_nodes(int min_nodesize) {
    for (size_t i = 0; i < node_vector.size(); i++)
      if (node_vector[i].in_use && node_vector[i].is_terminal)
        if ((int)sum(node_vector[i].observations) < min_nodesize) return true;
    return false;
  }

  // Generic indicator-weighted integrated log-likelihood + tree prior.
  double log_lik(arma::mat sigma_par, arma::mat sigma, double alpha, double beta,
                 arma::mat y_resid, arma::mat D) {
    double ll = 0.0;
    for (size_t i = 0; i < node_vector.size(); i++) {
      if (node_vector[i].in_use && node_vector[i].is_terminal) {
        double nj = sum(node_vector[i].observations);
        arma::mat nr = y_resid.rows(find(node_vector[i].observations == 1));
        arma::mat nd = D.rows(find(node_vector[i].observations == 1));
        arma::mat tricky(y_resid.n_cols, 1, arma::fill::zeros);
        arma::mat nd_t = nd.t(); arma::mat nr_t = nr.t();
        for (unsigned int t = 0; t < nr.n_rows; t++)
          tricky = tricky + nd_t.col(t) % (arma::inv(sigma) * nr_t.col(t));
        arma::mat s_inv = (nd.t()*nd) % arma::inv(sigma) + arma::inv(sigma_par);
        arma::mat p_j0 = arma::inv(s_inv) * tricky;
        double e1 = (-1.0*nj/2.0)*log(arma::det(sigma));
        double e2 = (-1.0/2.0)*log(arma::det(sigma_par));
        double e3 = (-1.0/2.0)*log(arma::det(arma::inv(sigma_par)+(nd.t()*nd)%arma::inv(sigma)));
        double e4 = (-1.0/2.0)*arma::accu((nr.t()*nr)%arma::inv(sigma));
        double e5 = arma::accu((1.0/2.0)*(p_j0.t())*(s_inv)*(p_j0));
        double e6 = log(1.0 - alpha*pow(1+floor(log2(i+1)), (-1*beta)));
        ll += e1+e2+e3+e4+e5+e6;
      } else if (node_vector[i].in_use && !node_vector[i].is_terminal) {
        ll += log(alpha) - beta*log(1+floor(log2(i+1)));
      }
    }
    return ll;
  }

  arma::mat get_predictions(int q) {
    int nobs = node_vector[0].observations.n_elem;
    arma::mat preds(nobs, q, arma::fill::zeros);
    for (size_t i = 0; i < node_vector.size(); i++)
      for (int j = 0; j < nobs; j++)
        if (node_vector[i].is_terminal && node_vector[i].in_use && node_vector[i].observations[j]==1)
          for (int k = 0; k < q; k++) preds(j,k) = node_vector[i].mu(k);
    return preds;
  }
  arma::mat get_test_predictions(int q) {
    int nobs = node_vector[0].test_observations.n_elem;
    arma::mat preds(nobs, q, arma::fill::zeros);
    for (size_t i = 0; i < node_vector.size(); i++)
      for (int j = 0; j < nobs; j++)
        if (node_vector[i].is_terminal && node_vector[i].in_use && node_vector[i].test_observations[j]==1)
          for (int k = 0; k < q; k++) preds(j,k) = node_vector[i].mu(k);
    return preds;
  }
};

// ============================= Forest ========================================
class Forest {
public:
  std::vector<Tree> tree_vector;
  Forest(int num_trees=1, int num_nodes=1, int num_obs=1, int num_test_obs=1) {
    tree_vector.resize(num_trees);
    for (int i = 0; i < num_trees; i++)
      tree_vector[i] = Tree(num_nodes, num_obs, num_test_obs);
  }
};

// ============================ helpers ========================================
arma::mat sum_over_cube_without_slice(const arma::cube& c, int removed) {
  arma::mat res(c.n_rows, c.n_cols, arma::fill::zeros);
  for (unsigned int i = 0; i < c.n_slices; i++) if ((int)i != removed) res += c.slice(i);
  return res;
}
arma::mat sample_sigma(double n, int v_0, arma::mat y, arma::mat preds, arma::mat sigma_0) {
  arma::mat resid = y - preds;
  return riwish(v_0 + n, sigma_0 + resid.t()*resid);
}
void collect_tree_nodes(const Tree& t, int iter, int tree_num, int comp, std::vector<double>& out) {
  int n = t.node_vector.size(); if (n == 0) return;
  std::vector<int> stack; stack.push_back(0);
  while (!stack.empty()) {
    int ni = stack.back(); stack.pop_back();
    double var=0, val=0, it=0, tr=0, cp=0;
    if (!t.node_vector[ni].is_terminal) {
      var = t.node_vector[ni].variable + 1; val = t.node_vector[ni].split_val;
      it = iter+1; tr = tree_num+1; cp = comp+1;
    } else {
      var = 0; val = t.node_vector[ni].mu[0]; it = iter+1; tr = tree_num+1; cp = comp+1;
    }
    out.push_back(var); out.push_back(val); out.push_back(it); out.push_back(tr); out.push_back(cp);
    int lc = 2*ni+1, rc = 2*ni+2;
    if (rc < n) { if (t.node_vector[rc].in_use) stack.push_back(rc); }
    if (lc < n) { if (t.node_vector[lc].in_use) stack.push_back(lc); }
  }
}

// ====================== main sampler: fast_bart_multi ========================
// [[Rcpp::export]]
List fast_bart_multi(arma::mat y,
                     List X_list,           // K training moderator matrices (n x p_k)
                     List X_test_list,      // K test moderator matrices    (n_test x p_k)
                     List D_list,           // K indicator matrices         (n x q)
                     NumericVector alpha,   // K tree-prior alphas
                     NumericVector beta,    // K tree-prior betas
                     List sigma_par_list,   // K prior covariances          (q x q)
                     NumericVector n_tree,  // K forest sizes
                     LogicalVector add_mean,// K: TRUE for the prognostic comp only
                     int v_0, arma::mat sigma_0,
                     int n_iter, int min_nodesize, int n_burn, int keep_every,
                     NumericVector tree_iters)
{
  int K = X_list.size();
  int n = y.n_rows, q = y.n_cols;

  // materialise inputs
  std::vector<arma::mat> X(K), Xt(K), D(K), Sp(K);
  std::vector<int> m(K);
  int n_test = 0;
  for (int k = 0; k < K; k++) {
    X[k]  = as<arma::mat>(X_list[k]);
    Xt[k] = as<arma::mat>(X_test_list[k]);
    D[k]  = as<arma::mat>(D_list[k]);
    Sp[k] = as<arma::mat>(sigma_par_list[k]);
    m[k]  = n_tree[k];
    n_test = Xt[k].n_rows;
  }

  // how many posterior draws to keep
  int num_keep = 0;
  for (int i = 0; i < n_iter; i++)
    if (i >= (n_burn-1) && (i-n_burn) % keep_every == 0) num_keep++;

  arma::mat sigma(q, q, arma::fill::eye);

  // scale y (unit sd per column)
  arma::rowvec col_means = mean(y, 0);
  arma::rowvec col_stdev = stddev(y, 0);
  arma::mat y_scaled = y.each_row() - col_means;
  y_scaled.each_row() /= col_stdev;

  // per-component storage
  std::vector<arma::cube> tp(K), tpt(K);          // per-tree preds (train / test)
  std::vector<arma::cube> keep(K), keep_t(K);     // kept posterior preds
  std::vector<Forest> forests(K);
  std::vector<std::vector<double>> flat(K);
  for (int k = 0; k < K; k++) {
    tp[k]  = arma::cube(n, q, m[k], arma::fill::zeros);
    tpt[k] = arma::cube(n_test, q, m[k], arma::fill::zeros);
    keep[k]   = arma::cube(n, q, num_keep, arma::fill::zeros);
    keep_t[k] = arma::cube(n_test, q, num_keep, arma::fill::zeros);
    forests[k] = Forest(m[k], 1, n, n_test);
  }
  arma::cube sigmas(q, q, num_keep, arma::fill::zeros);

  StringVector choices = {"Grow", "Prune", "Change", "Swap"};
  int num_kept = 0;

  for (int iter = 0; iter < n_iter; iter++) {

    // total fit from all components (D (.) F) -- recomputed once per sweep
    for (int k = 0; k < K; k++) {
      for (int tree_num = 0; tree_num < m[k]; tree_num++) {
        // residual excluding this tree's contribution to component k
        arma::mat y_resid = y_scaled;
        for (int kk = 0; kk < K; kk++) {
          int rm = (kk == k) ? tree_num : -1;
          y_resid -= D[kk] % sum_over_cube_without_slice(tp[kk], rm);
        }
        String choice = sample(choices, 1)[0];
        Tree prop = Tree(forests[k].tree_vector[tree_num]);
        if (choice == "Grow")  prop.grow(X[k], Xt[k], X[k].n_cols, min_nodesize);
        if (choice == "Prune") prop.prune();
        if (choice == "Change"){ prop.change(X[k], X[k].n_cols); prop.change_update(X[k], Xt[k]); }
        if (choice == "Swap")  { prop.swap(); prop.change_update(X[k], Xt[k]); }

        if (!prop.has_empty_nodes(min_nodesize)) {
          double lnew = prop.log_lik(Sp[k], sigma, alpha[k], beta[k], y_resid, D[k]);
          double lold = forests[k].tree_vector[tree_num].log_lik(Sp[k], sigma, alpha[k], beta[k], y_resid, D[k]);
          if (exp(lnew - lold) > R::runif(0,1))
            forests[k].tree_vector[tree_num] = Tree(prop);
        }
        forests[k].tree_vector[tree_num].update_nodes(sigma, Sp[k], y_resid, D[k]);
        tp[k].slice(tree_num)  = forests[k].tree_vector[tree_num].get_predictions(q);
        tpt[k].slice(tree_num) = forests[k].tree_vector[tree_num].get_test_predictions(q);
      }
    }

    // store tree structures at requested iterations
    for (int c = 0; c < tree_iters.size(); c++)
      if (tree_iters[c] == iter+1)
        for (int k = 0; k < K; k++)
          for (int tree_num = 0; tree_num < m[k]; tree_num++)
            collect_tree_nodes(forests[k].tree_vector[tree_num], iter, tree_num, k, flat[k]);

    Rcpp::Rcout << "Iteration " << iter+1 << " of " << n_iter << " ("
                << (float)(iter+1)/(float)n_iter*100 << "%)        \r";
    Rcpp::Rcout.flush();

    // total fitted mean (scaled space)
    arma::mat total_fit(n, q, arma::fill::zeros);
    for (int k = 0; k < K; k++) total_fit += D[k] % sum_over_cube_without_slice(tp[k], -1);

    bool do_keep = (iter >= (n_burn-1)) && ((iter-n_burn) % keep_every == 0);
    if (do_keep) {
      for (int k = 0; k < K; k++) {
        arma::mat fk  = sum_over_cube_without_slice(tp[k], -1);
        arma::mat fkt = sum_over_cube_without_slice(tpt[k], -1);
        for (unsigned int i = 0; i < q; i++) {
          keep[k].slice(num_kept).col(i)   = fk.col(i)  * col_stdev(i);
          keep_t[k].slice(num_kept).col(i) = fkt.col(i) * col_stdev(i);
          if (add_mean[k]) {
            keep[k].slice(num_kept).col(i)   += col_means(i);
            keep_t[k].slice(num_kept).col(i) += col_means(i);
          }
        }
      }
    }

    sigma = sample_sigma(n, v_0, y_scaled, total_fit, sigma_0);

    if (do_keep) {
      // Sigma on the original (unscaled) scale for reporting
      arma::mat tot_unscaled(n, q, arma::fill::zeros);
      for (int k = 0; k < K; k++) tot_unscaled += D[k] % keep[k].slice(num_kept);
      sigmas.slice(num_kept) = sample_sigma(n, v_0, y, tot_unscaled, sigma_0);
      num_kept++;
    }
  }
  Rcpp::Rcout << "\n";

  // assemble outputs
  List preds(K), preds_test(K), tree_df(K);
  for (int k = 0; k < K; k++) {
    preds[k] = keep[k];
    preds_test[k] = keep_t[k];
    int ncols = flat[k].size() / 5;
    NumericMatrix tdf(5, ncols + 1);
    for (int r = 0; r < 5; r++) tdf(r, 0) = -2;
    for (int j = 0; j < ncols; j++)
      for (int r = 0; r < 5; r++) tdf(r, j+1) = flat[k][5*j + r];
    tree_df[k] = tdf;
  }

  return List::create(
    Named("predictions") = preds,          // K-list of (n x q x draws) cubes
    Named("predictions_test") = preds_test,// K-list of (n_test x q x draws) cubes
    Named("sigmas") = sigmas,
    Named("tree_df") = tree_df,            // K-list of 5 x (nodes+1): var,val,iter,tree,comp
    Named("K") = K);
}
