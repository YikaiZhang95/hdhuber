// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
extern "C" {
  #include <R_ext/BLAS.h>   // for F77_CALL(daxpy)
}

using namespace Rcpp;
using namespace std;

// --------------------------- small utilities ---------------------------

inline double clip_huber(double r, double h) {
  if (r >  h) return  h;
  if (r < -h) return -h;
  return r;
}

// return a with the sign of b (Fortran SIGN(a,b))
inline double sign_copy(double a, double b) {
  return (b >= 0.0) ? std::abs(a) : -std::abs(a);
}

// elementwise max with 0 for matrix pf
static void clamp_nonneg(NumericMatrix pf) {
  const int p = pf.nrow(), q = pf.ncol();
  for (int j = 0; j < p; ++j)
    for (int l = 0; l < q; ++l)
      if (pf(j, l) < 0.0) pf(j, l) = 0.0;
}

// elementwise max with 0 for vector pf2
static void clamp_nonneg_vec(NumericVector v) {
  for (int i = 0; i < v.size(); ++i)
    if (v[i] < 0.0) v[i] = 0.0;
}

// gradient helper used for strong rule / KKT checks
static NumericVector huber_drv_internal(const NumericMatrix& X,
                                        const NumericVector& r,
                                        const double hval) {
  const int n = X.nrow(), p = X.ncol();
  const double ninv = 1.0 / static_cast<double>(n);
  NumericVector vl(p);
  // precompute clipped residuals
  std::vector<double> dl(n);
  for (int i = 0; i < n; ++i) dl[i] = clip_huber(r[i], hval);
  // vl_j = -(1/n) * <dl, X[,j]>
  for (int j = 0; j < p; ++j) {
    double s = 0.0;
    const double* colj = &X(0, j);
    for (int i = 0; i < n; ++i) s += colj[i] * dl[i];
    vl[j] = -s * ninv;
  }
  return vl;
}

// --------------------------- chkvars ---------------------------
// Fortran: CALL chkvars(nobs, nvars, x, ju)
// Here: return ju: 1 if column j has non-zero variance, else 0
static IntegerVector chkvars_cpp(const NumericMatrix& X) {
  const int n = X.nrow(), p = X.ncol();
  IntegerVector ju(p);
  for (int j = 0; j < p; ++j) {
    const double* col = &X(0, j);
    double mu = 0.0;
    for (int i = 0; i < n; ++i) mu += col[i];
    mu /= static_cast<double>(n);
    double ss = 0.0;
    for (int i = 0; i < n; ++i) {
      const double z = col[i] - mu;
      ss += z * z;
    }
    ju[j] = (ss > 0.0) ? 1 : 0;
  }
  return ju;
}

// --------------------------- Standard ---------------------------
// Fortran: CALL Standard(nobs, nvars, x, ju, isd, xmean, xnorm, maj)
//
// Behavior implemented:
// - Center X columns by mean always (intercept handled separately).
// - If isd==1, divide by SD; maj[j] = 1 (since mean-square becomes 1).
// - If isd==0, do not scale; maj[j] = mean of squared centered column.
//
// NOTE: We do NOT modify ju here; zero-variance columns should already be 0.
//       xmean = column means (pre-centering), xnorm = SD used if isd==1 (else SD
//       is returned for post-unscaling), maj = per-feature majorization const.
static void standardize_cpp(NumericMatrix& X,
                            const IntegerVector& ju,
                            const int isd,
                            NumericVector& xmean,
                            NumericVector& xnorm,
                            NumericVector& maj) {
  const int n = X.nrow(), p = X.ncol();
  xmean = NumericVector(p);
  xnorm = NumericVector(p);
  maj   = NumericVector(p);

  for (int j = 0; j < p; ++j) {
    if (ju[j] == 0) { // not used
      xmean[j] = 0.0;
      xnorm[j] = 1.0;
      maj[j]   = 0.0;
      continue;
    }
    double mu = 0.0;
    for (int i = 0; i < n; ++i) mu += X(i, j);
    mu /= static_cast<double>(n);
    xmean[j] = mu;

    // center
    for (int i = 0; i < n; ++i) X(i, j) -= mu;

    // variance / maj
    double ss = 0.0;
    for (int i = 0; i < n; ++i) ss += X(i, j) * X(i, j);
    const double msq = ss / static_cast<double>(n);
    const double sd  = (msq > 0.0) ? std::sqrt(msq) : 0.0;

    xnorm[j] = (sd > 0.0) ? sd : 1.0; // store sd for later unscale

    if (isd == 1 && sd > 0.0) {
      // scale to unit sd
      const double invsd = 1.0 / sd;
      for (int i = 0; i < n; ++i) X(i, j) *= invsd;
      maj[j] = 1.0; // mean-square after scaling is 1
    } else {
      maj[j] = msq; // mean-square (centered) if not scaling
    }
  }
}

// --------------------------- huber_path core ---------------------------

static List huber_path_core(double alpha,
                            double lam2_in,
                            double hval,
                            NumericVector maj, // will be scaled by mval (like Fortran)
                            double mval,
                            const NumericMatrix& X,
                            const NumericVector& y,
                            const IntegerVector& ju,
                            int pfncol,
                            const NumericMatrix& pf,
                            const NumericVector& pf2,
                            int dfmax,
                            int pmax,
                            int nlam,
                            double flmin_in,
                            const NumericVector& ulam,
                            double eps,
                            int maxit,
                            int istrong) {

  // constants
  const double BIG  = 9.9e30;
  const double MFL  = 1.0e-6;
  const int    MNLAM = 6;

  const int n = X.nrow();
  const int p = X.ncol();

  // working copies
  NumericVector r = clone(y);            // residual/margin
  NumericVector ga(p), vl(p);
  NumericVector b0_out(nlam);            // per-lambda intercept
  NumericMatrix beta_out(pmax, nlam);    // compressed coefficients
  IntegerVector m(pmax);                 // index of active variables (1-based when returned)
  IntegerVector mm(p);                   // position map: 0 if not in m; else 1..ni
  IntegerVector nbeta_out(nlam);         // number of compressed coefs per lambda
  NumericVector alam(nlam);              // lambda values
  IntegerVector jxx(p);                  // strong set flags

  // initializations
  double al = 0.0;
  double al0 = 0.0;
  double alf = 0.01; // will be updated when flmin<1
  int npass = 0;
  int nalam = 0;
  int jerr  = 0;

  const double* Xbase = X.begin();
  auto colptr = [&](int j)->const double* { return Xbase + (size_t)j*n; };
  // effective maj scaling
  for (int j = 0; j < p; ++j) maj[j] *= mval;

  // strong rule switches
  if (istrong == 1) {
    std::fill(jxx.begin(), jxx.end(), 0);
  } else {
    std::fill(jxx.begin(), jxx.end(), 1);
  }

  // path grid rule
  int mnl = std::min(MNLAM, nlam);
  double flmin = flmin_in;
  if (flmin < 1.0) {
    flmin = std::max(MFL, flmin);
    if (nlam > 1) alf = std::pow(flmin, 1.0 / (static_cast<double>(nlam) - 1.0));
  }

  // coefficients (excluding intercept) and their old copies
  std::vector<double> b(p, 0.0), oldb_all(p, 0.0);
  double b0 = 0.0, oldb0 = 0.0;

  int ni = 0; // active count

  // lambdas loop
  for (int l = 0; l < nlam; ++l) {
    const int pfl = (pfncol == 1 ? 0 : l); // column in pf to use

    al0 = al;
    if (flmin >= 1.0) {
      al = ulam[l];
    } else {
      if (l > 1) {
        al = al * alf;
      } else if (l == 0) {
        al = BIG;
      } else if (l == 1) {
        al0 = 0.0;
        vl  = huber_drv_internal(X, r, hval);
        for (int j = 0; j < p; ++j) ga[j] = std::abs(vl[j]);
        double maxv = 0.0;
        for (int j = 0; j < p; ++j) {
          if (ju[j] == 0) continue;
          const double pf_j = pf(j, pfl);
          if (pf_j > 0.0) {
            maxv = std::max(maxv, ga[j] / pf_j);
          }
        }
        al = maxv * alf;
      }
    }

    // lam2 effective (depends on current lambda if alpha provided)
    double lam2 = lam2_in;
    if (alpha != -1.0) {
      // Fortran: lam2 = al * (1 - alpha) * 0.5D0 / alpha
      lam2 = al * (1.0 - alpha) * 0.5 / alpha;
    }

    // strong rule update
    const double tlam = 2.0 * al - al0;
    for (int j = 0; j < p; ++j) {
      if (jxx[j] == 1) continue; // already in
      if (ga[j] > pf(j, pfl) * tlam) jxx[j] = 1;
    }

    // ---------------- outer loop ----------------
    for (;;) {
      // save old
      oldb0 = b0;
      if (ni > 0) {
        for (int t = 0; t < ni; ++t) {
          const int k = m[t] - 1; // m stores 1-based, keep 0-based here
          oldb_all[k] = b[k];
        }
      }

      // ------------- middle loop -------------
      for (;;) {
        npass += 1;
        double dif = 0.0;

        // (1) update all beta_j in the current strong set
        std::vector<double> dl(n); // clipped residual
        for (int k = 0; k < p; ++k) {
          if (jxx[k] == 0) continue;
          if (ju[k] == 0) continue;

          const double oldbk = b[k];
          // compute u = maj(k)*b(k) - (1/n) * sum_i V'(r_i) * x_ik
          double u = 0.0;
          for (int i = 0; i < n; ++i) {
            dl[i] = clip_huber(r[i], hval);
            u -= dl[i] * X(i, k);
          }
          u = maj[k] * b[k] - u / static_cast<double>(n);
          double v = al * pf(k, pfl);
          v = std::abs(u) - v;
          
          if (v > 0.0) {
            b[k] = sign_copy(v, u) / (maj[k] + pf2[k] * lam2);
          } else {
            b[k] = 0.0;
          }

          const double d = b[k] - oldbk;
          if (std::abs(d) > 0.0) {
            dif = std::max(dif, d * d);
            for (int i = 0; i < n; ++i) r[i] -= X(i, k) * d;

            if (mm[k] == 0) {
              ni += 1;
              if (ni > pmax) goto exceed_pmax_middle;
              mm[k] = ni;
              m[ni - 1] = k + 1; // store 1-based index
            }
          }
        }

        // (2) update intercept
        {
          double d = 0.0;
          for (int i = 0; i < n; ++i) {
            const double dl_i = clip_huber(r[i], hval);
            d -= dl_i;
          }
          d = -d / static_cast<double>(n) / mval;
          if (d != 0.0) {
            b0  += d;
            for (int i = 0; i < n; ++i) r[i] -= d;
            dif = std::max(dif, d * d);
          }
        }

        if (dif < eps) break;
        if (npass > maxit) { jerr = -1; goto done_path; }

        // ------------- inner loop -------------
        for (;;) {
          npass += 1;
          double dif = 0.0;

          // (1) update beta_j in active set only
          for (int t = 0; t < ni; ++t) {
            const int k = m[t] - 1;
            const double* xk = colptr(k);
            const double oldbk = b[k];

            double u = 0.0;
            for (int i = 0; i < n; ++i) {
              const double dl_i = clip_huber(r[i], hval);
              u -= dl_i * X(i, k);
            }
            u = maj[k] * b[k] - u / static_cast<double>(n);

            double v = al * pf(k, pfl);
            v = std::abs(u) - v;

            if (v > 0.0) {
              b[k] = sign_copy(v, u) / (maj[k] + pf2[k] * lam2);
            } else {
              b[k] = 0.0;
            }

            const double d = b[k] - oldbk;
            if (std::abs(d) > 0.0) {
              dif = std::max(dif, d * d);
              // for (int i = 0; i < n; ++i) r[i] -= X(i, k) * d;
              int nn = n, inc = 1; double nd = -d;
              F77_CALL(daxpy)(&nn, &nd, xk, &inc, REAL(r), &inc);
            }
          }

          // (2) update intercept
          {
            double d = 0.0;
            for (int i = 0; i < n; ++i) {
              const double dl_i = clip_huber(r[i], hval);
              d -= dl_i;
            }
            d = -d / static_cast<double>(n) / mval;
            if (d != 0.0) {
              b0  += d;
              for (int i = 0; i < n; ++i) r[i] -= d;
              dif = std::max(dif, d * d);
            }
            
          }

          if (dif < eps) break;
          if (npass > maxit) { jerr = -1; goto done_path; }
        } // inner
      } // middle

      // ----- final check: did any active coef/intercept change beyond eps?
      {
        int jx = 0;
        if ((b0 - oldb0) * (b0 - oldb0) >= eps) jx = 1;
        if (!jx) {
          for (int t = 0; t < ni; ++t) {
            const int k = m[t] - 1;
            const double diff = b[k] - oldb_all[k];
            if (diff * diff >= eps) { jx = 1; break; }
          }
        }
        if (jx) continue;

        // KKT check for discarded variables
        vl = huber_drv_internal(X, r, hval);
        for (int j = 0; j < p; ++j) ga[j] = std::abs(vl[j]);
        for (int j = 0; j < p; ++j) {
          if (jxx[j] == 1) continue;
          if (ga[j] > al * pf(j, pfl)) { jxx[j] = 1; jx = 1; }
        }
        if (jx) continue; // need another sweep
      }

      break; // outer converged
    } // outer

    // ---- save results for this lambda ----
    if (ni > pmax) { jerr = -10000 - (l + 1); goto done_path; }
    if (ni > 0) {
      for (int t = 0; t < ni; ++t) {
        const int k = m[t] - 1;
        beta_out(t, l) = b[k];
      }
    }
    nbeta_out[l] = ni;
    b0_out[l]    = b0;
    alam[l]      = al;
    nalam        = l + 1;
    // early stopping by dfmax (only after a few lambdas and when flmin<1)
    if (l >= (mnl - 1) && flmin < 1.0) {
      int me = 0;
      for (int t = 0; t < ni; ++t) if (beta_out(t, l) != 0.0) ++me;
      if (me > dfmax) break;
    }

    continue;

  exceed_pmax_middle:
    jerr = -10000 - (l + 1);
    break;
  } // for each lambda

done_path:
  // Return 1-based m like Fortran/R
  IntegerVector ibeta = clone(m);

  return List::create(
    _["nalam"] = nalam,
    _["b0"]    = b0_out,
    _["beta"]  = beta_out,
    _["ibeta"] = ibeta,
    _["nbeta"] = nbeta_out,
    _["alam"]  = alam,
    _["npass"] = npass,
    _["jerr"]  = jerr
  );
}

// --------------------------- Exported: huber_path_cpp ---------------------------

/**
 * R interface mirroring the Fortran huber_path interface (but nobs/nvars
 * are inferred from X). Expects:
 *  - X, y already standardized/centered consistently with maj & mval,
 *  - ju indicates usable variables (1/0),
 *  - pf: (nvars x pfncol), pf2: length nvars,
 *  - ulam length nlam.
 *
 * Returns a list: nalam, b0, beta, ibeta, nbeta, alam, npass, jerr.
 */
// [[Rcpp::export]]
List huber_path_cpp(double alpha,
                    double lam2,
                    double hval,
                    NumericVector maj,
                    double mval,
                    NumericMatrix X,
                    NumericVector y,
                    IntegerVector ju,
                    int pfncol,
                    NumericMatrix pf,
                    NumericVector pf2,
                    int dfmax,
                    int pmax,
                    int nlam,
                    double flmin,
                    NumericVector ulam,
                    double eps,
                    int maxit,
                    int istrong) {
  // simple shape checks
  const int n = X.nrow(), p = X.ncol();
  if (y.size() != n) stop("Length of y must equal nrow(X).");
  if (maj.size() != p) stop("Length of maj must equal ncol(X).");
  if (pf.nrow() != p) stop("pf must have nrow = ncol(X).");
  if (!(pfncol == 1 || pf.ncol() == nlam))
    stop("pfncol must be 1 or equal to nlam and match pf.ncol().");
  if (pf2.size() != p) stop("pf2 must have length ncol(X).");
  if (ulam.size() != nlam && flmin >= 1.0)
    stop("When flmin >= 1, length(ulam) must equal nlam.");

  return huber_path_core(alpha, lam2, hval, maj, mval,
                         X, y, ju, pfncol, pf, pf2,
                         dfmax, pmax, nlam, flmin, ulam,
                         eps, maxit, istrong);
}

// --------------------------- Exported: huber_cd_cpp ---------------------------

/**
 * High-level routine that matches Fortran huber_cd:
 *
 * Inputs (scalar args as in Fortran):
 *  alpha, lam2, hval, nobs, nvars, X, y, jd, pfncol, pf, pf2,
 *  dfmax, pmax, nlam, flmin, ulam, eps, isd, maxit, istrong.
 *
 * Returns a list containing:
 *  nalam, b0, beta, ibeta, nbeta, alam, npass, jerr
 *  (beta is already unscaled/uncentered to original X scale).
 */
// [[Rcpp::export]]
List huber_cd_cpp(double alpha,
                  double lam2,
                  double hval,
                  int nobs,
                  int nvars,
                  NumericMatrix X,
                  NumericVector y,
                  IntegerVector jd,        // deletion flags; 1st element is count, followed by 1-based indices
                  int pfncol,
                  NumericMatrix pf,
                  NumericVector pf2,
                  int dfmax,
                  int pmax,
                  int nlam,
                  double flmin,
                  NumericVector ulam,
                  double eps,
                  int isd,
                  int maxit,
                  int istrong) {

  // shape checks
  if (X.nrow() != nobs || X.ncol() != nvars) stop("X must be nobs x nvars.");
  if (y.size() != nobs) stop("y must have length nobs.");
  if (pf.nrow() != nvars) stop("pf must have nrow=nvars.");
  if (!(pfncol == 1 || pf.ncol() == nlam))
    stop("pfncol must be 1 or equal to nlam and match pf.ncol().");
  if (pf2.size() != nvars) stop("pf2 must have length nvars.");

  // 1) chkvars -> ju
  IntegerVector ju = chkvars_cpp(X);

  // delete variables specified in jd (1-based indices)
  if (jd.size() > 0 && jd[0] > 0) {
    const int deln = jd[0];
    if (jd.size() < deln + 1) stop("jd length inconsistent with jd[1].");
    for (int k = 0; k < deln; ++k) {
      const int idx1 = jd[k + 1]; // 1-based
      if (idx1 >= 1 && idx1 <= nvars) ju[idx1 - 1] = 0;
    }
  }

  // error checks like Fortran
  {
    int anyju = 0;
    for (int j = 0; j < nvars; ++j) if (ju[j] != 0) { anyju = 1; break; }
    if (!anyju) {
      return List::create(_["nalam"]=0, _["b0"]=NumericVector(0), _["beta"]=NumericMatrix(pmax, 0),
                          _["ibeta"]=IntegerVector(pmax), _["nbeta"]=IntegerVector(0),
                          _["alam"]=NumericVector(0), _["npass"]=0, _["jerr"]=7777);
    }
  }

  // positivity checks for pf/pf2
  {
    double maxpf = R_NegInf, maxpf2 = R_NegInf;
    for (int j = 0; j < pf.nrow(); ++j)
      for (int c = 0; c < pf.ncol(); ++c)
        if (pf(j, c) > maxpf) maxpf = pf(j, c);
    for (int j = 0; j < pf2.size(); ++j)
      if (pf2[j] > maxpf2) maxpf2 = pf2[j];

    if (!(maxpf > 0.0)) {
      return List::create(_["nalam"]=0, _["b0"]=NumericVector(0), _["beta"]=NumericMatrix(pmax, 0),
                          _["ibeta"]=IntegerVector(pmax), _["nbeta"]=IntegerVector(0),
                          _["alam"]=NumericVector(0), _["npass"]=0, _["jerr"]=10000);
    }
    if (!(maxpf2 > 0.0)) {
      return List::create(_["nalam"]=0, _["b0"]=NumericVector(0), _["beta"]=NumericMatrix(pmax, 0),
                          _["ibeta"]=IntegerVector(pmax), _["nbeta"]=IntegerVector(0),
                          _["alam"]=NumericVector(0), _["npass"]=0, _["jerr"]=10000);
    }
  }

  // clamp pf and pf2 to be nonnegative (like Fortran: pf = max(0, pf))
  clamp_nonneg(pf);
  clamp_nonneg_vec(pf2);

  // 2) Standardize X (center always; scale if isd==1). Collect xmean, xnorm, maj.
  NumericMatrix Xstd = clone(X);
  NumericVector xmean, xnorm, maj;
  standardize_cpp(Xstd, ju, isd, xmean, xnorm, maj);

  // 3) Path solver on standardized X
  List res = huber_path_core(alpha, lam2, hval, clone(maj), 1.0,
                             Xstd, y, ju, pfncol, pf, pf2,
                             dfmax, pmax, nlam, flmin, ulam,
                             eps, maxit, istrong);

  int jerr = as<int>(res["jerr"]);
  if (jerr > 0) return res; // fatal; pass through

  // 4) Post-process to original scale (like Fortran block after huber_path)
  int nalam = as<int>(res["nalam"]);
  IntegerVector ibeta = res["ibeta"];            // 1-based
  IntegerVector nbeta = res["nbeta"];
  NumericMatrix beta  = res["beta"];             // compressed rows
  NumericVector b0    = res["b0"];

  if (nalam > 0) {
    for (int l = 0; l < nalam; ++l) {
      int nk = nbeta[l];
      if (isd == 1) {
        for (int j = 0; j < nk; ++j) {
          const int col = ibeta[j] - 1; // 0-based
          if (xnorm[col] != 0.0) beta(j, l) /= xnorm[col];
        }
      }
      // b0(l) = b0(l) - dot(beta(1:nk,l), xmean(ibeta(1:nk)))
      double adj = 0.0;
      for (int j = 0; j < nk; ++j) {
        const int col = ibeta[j] - 1; // 0-based
        adj += beta(j, l) * xmean[col];
      }
      b0[l] -= adj;
    }
  }
  NumericMatrix beta_full(nvars, nalam);
  std::fill(beta_full.begin(), beta_full.end(), 0.0);
  for (int l = 0; l < nalam; ++l) {
    const int nk = nbeta[l];
    for (int j = 0; j < nk; ++j) {
      const int var = ibeta[j] - 1;
      beta_full(var, l) = beta(j, l);
    }
  }

  res["beta_full"] = beta_full;
  res["beta"] = beta;
  res["b0"]   = b0;
  return res;
}

// --------------------------- Exported: huber_drv_cpp ---------------------------

/**
 * Gradient used throughout:
 *   vl_j = -(1/n) * sum_i X_{ij} * clip(r_i, -hval, +hval)
 */
// [[Rcpp::export]]
NumericVector huber_drv_cpp(const NumericMatrix& X,
                            const NumericVector& r,
                            const double hval) {
  const int n = X.nrow();
  if (r.size() != n) stop("Length of r must match nrow(X).");
  return huber_drv_internal(X, r, hval);
}
