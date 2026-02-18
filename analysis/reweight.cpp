// reweight_module.cpp
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <vector>
#include <cmath>
#include <random>
#include <stdexcept>

#include "ranlxs.h"

namespace py = pybind11;
using ld = long double;

// -------------------------
// Autocorrelation function
// -------------------------
void AutoCorr(const std::vector<ld>& data, ld &autocorr, ld &err) {
    int len = data.size();
    if(len == 0) {
        autocorr = err = 0.0L;
        return;
    }

    ld C0 = 0.0L, Ct = 0.0L, rho = 0.0L, tint = 0.5L;
    ld avr = 0.0L, f2 = 0.0L, f = 0.0L;
    // t=0: compute C0
    for (int i = 0; i < len; i++){
        avr += data[i];
        f2 += data[i] * data[i];
    }
    f = 2.0L * avr;
    avr /= (ld) len;
    C0 = (f2/(ld)len) - avr * avr;
    
    bool valid = false;
    int M;
    for (M = 1; M < len; M++) {
        f2 = 0.0L; f = 0.0L; avr = 0.0L;
        for (int i = 0; i < len - M; i++){
            f2 += data[i] * data[i+M];
            f += data[i] + data[i+M];
            avr += data[i];
        }
        int n_M = len - M;
        avr /= (ld)n_M;
        Ct = (f2/(ld)n_M) + avr*(avr - (f/(ld)n_M));
        rho = Ct / C0;
        tint += rho;
        if(M > 4.0L * tint) { valid = true; break; }
    }
    if(valid) {
        autocorr = tint;
        err = std::sqrtl(2.0L*(2*M+1)/(ld)len) * tint;
    } else {
        autocorr = err = 0.0L; // not valid: no proper determination
    }
}

// -------------------------
// Helper: log-sum-exp for two numbers
// -------------------------
ld logsumexp(ld a, ld b) {
    // returns log(exp(a) + exp(b)) in a numerically stable way
    if (a > b) {
        ld e_val = std::expl(b - a);
        if (e_val == 0.0L) {
            py::print("Underflow in logsumexp: exponent = ", (b - a),
                      py::arg("flush")=true);
        }
        return a + std::log1pl(e_val);
    } else {
        ld e_val = std::expl(a - b);
        if (e_val == 0.0L) {
            py::print("Underflow in logsumexp: exponent = ", (a - b),
                      py::arg("flush")=true);
        }
        return b + std::log1pl(e_val);
    }
}


// -------------------------
// LnZ: Computes the log normalization (free energy) 
// -------------------------
ld LnZ(ld beta, int nrun, const std::vector<ld>& Beta, 
       const std::vector< std::vector<ld> >& d, 
       const std::vector<ld>& LogZ) {
    std::vector<ld> LogN(nrun);
    for (int i = 0; i < nrun; i++) {
        LogN[i] = std::logl((ld)d[i].size());
    }
    bool init = true;
    ld res = 0.0L;
    // Loop over all simulation runs and their samples
    for (int i = 0; i < nrun; i++) {
        for (size_t idx = 0; idx < d[i].size(); idx++) {
            ld E = d[i][idx];
            // LogDen is log of the denominator in self-consistency equation
            ld LogDen = LogN[0] - LogZ[0] + (beta - Beta[0]) * E;
            for (int j = 1; j < nrun; j++) {
                ld term = LogN[j] - LogZ[j] + (beta - Beta[j]) * E;
                LogDen = logsumexp(LogDen, term);
            }
            if(init) {
                res = -LogDen;
                init = false;
            } 
            else {
                ld diff = res + LogDen;
                if (diff > 0) {
                    ld tmp = std::expl(-diff);
                    if (tmp == 0.0L) {
                        py::print("Underflow in std::expl(-diff): diff =", -diff,
                                  py::arg("flush")=true);
                    }
                    res += std::log1pl(tmp);
                } else {
                    ld tmp = std::expl(diff);
                    if (tmp == 0.0L) {
                        py::print("Underflow in std::expl(diff): diff =", diff,
                                  py::arg("flush")=true);
                    }
                    res = -LogDen + std::log1pl(tmp);
                }
            }            
        }
    }
    return res;
}

// -------------------------
// LnO_n: Similar to LnZ but with observable moments.
// -------------------------
ld LnO_n(const std::vector< std::vector<ld> >& O, ld n, ld beta, int nrun, 
           const std::vector<ld>& Beta, 
           const std::vector< std::vector<ld> >& d, 
           const std::vector<ld>& LogZ) {
    std::vector<ld> LogN(nrun);
    for (int i = 0; i < nrun; i++) {
        if(d[i].size() != O[i].size())
            throw std::runtime_error("Mismatch between energy and observable sample sizes.");
        LogN[i] = std::logl((ld)d[i].size());
    }
    bool init = true;
    ld res = 0.0L;
    for (int i = 0; i < nrun; i++) {
        for (size_t idx = 0; idx < d[i].size(); idx++) {
            ld E = d[i][idx];
            ld LogDen = LogN[0] - LogZ[0] + (beta - Beta[0]) * E;
            for (int j = 1; j < nrun; j++) {
                ld term = LogN[j] - LogZ[j] + (beta - Beta[j]) * E;
                LogDen = logsumexp(LogDen, term);
            }
            if (O[i][idx] <= 0.0L) {
                throw std::runtime_error("Observable has nonpositive value, log() is invalid.");
            }
            ld LogO = n * std::logl(O[i][idx]);
            LogDen -= LogO;
            if(init) {
                res = -LogDen;
                init = false;
            } 
            else {
                ld diff = res + LogDen;
                if (diff > 0) {
                    ld tmp = std::expl(-diff);
                    if (tmp == 0.0L) {
                        py::print("Underflow in std::expl(-diff): diff =", -diff,
                                  py::arg("flush")=true);
                    }
                    res += std::log1pl(tmp);
                } else {
                    ld tmp = std::expl(diff);
                    if (tmp == 0.0L) {
                        py::print("Underflow in std::expl(diff): diff =", diff,
                                  py::arg("flush")=true);
                    }
                    res = -LogDen + std::log1pl(tmp);
                }
            }
        }
    }
    return res;
}

// -------------------------
// MultiHistRw: Iteratively update LogZ until convergence.
// -------------------------
void MultiHistRw(int nrun, const std::vector<ld>& beta, 
                 const std::vector< std::vector<ld> >& d, 
                 std::vector<ld>& LogZ, ld tol) {
    std::vector<ld> new_LogZ(nrun, 0.0L);
    ld Delta2;
    do {
        // Shift LogZ to remove an arbitrary additive constant
        // ld A=(LogZ[0]+LogZ[nrun-1])/2.0;
        // note that this is the middle simulation
        ld A = LogZ[nrun/2];
        for (int k = 0; k < nrun; k++)
            LogZ[k] -= A;
        // Update new_LogZ for each simulation run's beta
        for (int k = 0; k < nrun; k++)
            new_LogZ[k] = LnZ(beta[k], nrun, beta, d, LogZ);
        // Compute convergence measure Delta2
        Delta2 = 0.0L;
        for (int k = 0; k < nrun; k++) {
            ld term = std::expm1l(new_LogZ[k] - LogZ[k]);
            Delta2 += term * term;
        }
        // Update LogZ
        for (int k = 0; k < nrun; k++)
            LogZ[k] = new_LogZ[k];
    } while (Delta2 > tol);
}

// -------------------------
// reweight: main function exposed to Python
// -------------------------
//
// reweight(O, E, betas, rw_betas, N_bootstraps, tol)
//   O: (K,N) numpy array of observables (float64)
//   E: (K,N) numpy array of energies (float64)
//   betas: (K,) numpy array of simulation betas (float64)
//   rw_betas: (M,) numpy array of target betas (float64)
//   N_bootstraps: number of bootstrap samples (default 500)
//   tol: tolerance for convergence in MultiHistRw (default 1e-14)
// 
// Returns a tuple of numpy arrays (O_means, O_means_err, O_vars, O_vars_err, 
// O_vars_argmax_bs)
// Each of shape (M,) for the first four outputs, and the last output shape (N_bootstrap,).
py::tuple reweight(py::array_t<double, py::array::c_style | py::array::forcecast> O_array,
                py::array_t<double, py::array::c_style | py::array::forcecast> E_array,
                py::array_t<double, py::array::c_style | py::array::forcecast> betas_array,
                py::array_t<double, py::array::c_style | py::array::forcecast> rw_betas_array,
                int N_bootstraps = 500,
                double tol = 1e-14) {

    // --- Convert inputs from numpy arrays to C++ vectors ---
    py::buffer_info O_buf = O_array.request();
    py::buffer_info E_buf = E_array.request();
    py::buffer_info betas_buf = betas_array.request();
    py::buffer_info rw_betas_buf = rw_betas_array.request();
    
    if (O_buf.ndim != 2 || E_buf.ndim != 2)
        throw std::runtime_error("O and E must be 2D arrays");
    if (betas_buf.ndim != 1)
        throw std::runtime_error("betas must be a 1D array");
    if (rw_betas_buf.ndim != 1)
        throw std::runtime_error("rw_betas must be a 1D array");
    
    int K = O_buf.shape[0]; // number of simulation runs
    int N = O_buf.shape[1]; // number of samples per run
    int M = rw_betas_buf.shape[0]; // number of target beta values
    
    if(E_buf.shape[0] != K || E_buf.shape[1] != N)
        throw std::runtime_error("O and E must have the same shape");
    if(betas_buf.shape[0] != K)
        throw std::runtime_error("Length of betas must equal first dimension of O and E");
    
    // Simulation betas
    std::vector<ld> sim_betas(K);
    double* betas_ptr = static_cast<double*>(betas_buf.ptr);
    for (int i = 0; i < K; i++)
        sim_betas[i] = betas_ptr[i];
    
    // Target reweighting betas
    std::vector<ld> target_betas(M);
    double* rw_betas_ptr = static_cast<double*>(rw_betas_buf.ptr);
    for (int i = 0; i < M; i++)
        target_betas[i] = rw_betas_ptr[i];
    
    // Build per-run data vectors (convert from double to long double)
    std::vector< std::vector<ld> > O_runs(K), E_runs(K);
    double* O_ptr = static_cast<double*>(O_buf.ptr);
    double* E_ptr = static_cast<double*>(E_buf.ptr);
    for (int i = 0; i < K; i++) {
        O_runs[i].resize(N);
        E_runs[i].resize(N);
        for (int j = 0; j < N; j++) {
            O_runs[i][j] = O_ptr[i * N + j];
            E_runs[i][j] = E_ptr[i * N + j];
        }
    }
    
    // Storage for bootstrap results:
    // For each bootstrap replicate (0...N_bootstraps-1) and each target beta,
    // we will compute a reweighted observable average and variance.
    std::vector< std::vector<ld> > O_mean_bootstrap(N_bootstraps, std::vector<ld>(M, 0.0L));
    std::vector< std::vector<ld> > O_var_bootstrap(N_bootstraps, std::vector<ld>(M, 0.0L));
    // Storage for the beta values (from rw_betas) corresponding to the maximum O_vars in each bootstrap replicate.
    std::vector<ld> argmax_beta_vals;

    
    // // Setup a random number generator RANDLXS
    rlxs_init(0,872348673);
    // Define a helper lambda to generate one random number in [0,1)
    auto get_rand = []() -> ld {
        double r;
        ranlxs(&r, 1); // ranlxs fills r with one random number in [0,1)
        return (ld) r;
    };
    
    // ----- Begin bootstrap loop -----
    for (int rep = 0; rep < N_bootstraps; rep++) {
        // For each simulation run, create a bootstrap sample (with effective decorrelation)
        std::vector< std::vector<ld> > O_bs(K), E_bs(K);
        for (int i = 0; i < K; i++) {
            int len = N;
            ld tau, tau_err;
            AutoCorr(O_runs[i], tau, tau_err);
            if (tau == 0.0L) {  
                py::print("Warning: Autocorrelation time is cannot be determined.", py::arg("flush") = true);
            }
            if (tau < 1.0L)
                tau = 1.0L;
            // Determine number of bootstrap samples for run i: use samples at intervals ~tau
            int n_bs = (int) std::ceil((ld)len / tau);
            O_bs[i].resize(n_bs);
            E_bs[i].resize(n_bs);
            for (int j = 0; j < n_bs; j++) {
                int idx = (int) std::floor(len * get_rand());
                if (idx >= len) idx = len - 1;
                O_bs[i][j] = O_runs[i][idx];
                E_bs[i][j] = E_runs[i][idx];
            }
        }
        
        // Reweighting: initialize LogZ for each simulation run to 0.
        std::vector<ld> LogZ(K, 0.0L);
        MultiHistRw(K, sim_betas, E_bs, LogZ, tol);

        // For each target beta, compute the reweighted observable average and variance.
        for (int m = 0; m < M; m++) {
            ld b = target_betas[m];
            ld FreeNRG = LnZ(b, K, sim_betas, E_bs, LogZ);
            ld LnO1 = LnO_n(O_bs, 1.0L, b, K, sim_betas, E_bs, LogZ);
            ld O_mean = std::expl(LnO1 - FreeNRG);
            if (O_mean == 0.0L) {
                py::print("Underflow detected in std::expl: exponent =", LnO1 - FreeNRG, py::arg("flush") = true);
            } else if (std::isinf(O_mean)) {
                py::print("Overflow detected in std::expl: exponent =", LnO1 - FreeNRG, py::arg("flush") = true);
            }
            ld LnO2 = LnO_n(O_bs, 2.0L, b, K, sim_betas, E_bs, LogZ);
            ld O2 = std::expl(LnO2 - FreeNRG);
            if (O_mean == 0.0L) {
                py::print("Underflow detected in std::expl: exponent =", LnO2 - FreeNRG, py::arg("flush") = true);
            } else if (std::isinf(O_mean)) {
                py::print("Overflow detected in std::expl: exponent =", LnO2 - FreeNRG, py::arg("flush") = true);
            }
            ld variance = O2 - O_mean * O_mean;
            O_mean_bootstrap[rep][m] = O_mean;
            O_var_bootstrap[rep][m] = variance;
        }

        // Compute the argmax: find the target beta (from rw_betas) that yields the maximum O_var for this bootstrap replicate.
        ld max_var = O_var_bootstrap[rep][0];
        int max_index = 0;
        for (int m = 1; m < M; m++) {
            if (O_var_bootstrap[rep][m] > max_var) {
                max_var = O_var_bootstrap[rep][m];
                max_index = m;
            }
        }
        argmax_beta_vals.push_back(target_betas[max_index]);
    }
    // ----- End bootstrap loop -----
    
    // ----- Compute bootstrap averages and errors for each target beta -----
    std::vector<double> O_means_rw(M, 0.0), O_means_rw_err(M, 0.0);
    std::vector<double> O_vars_rw(M, 0.0), O_vars_rw_err(M, 0.0);
    for (int m = 0; m < M; m++) {
        ld sum_mean = 0.0L, sum_mean2 = 0.0L;
        ld sum_var = 0.0L, sum_var2 = 0.0L;
        for (int rep = 0; rep < N_bootstraps; rep++) {
            ld om = O_mean_bootstrap[rep][m];
            ld ov = O_var_bootstrap[rep][m];
            sum_mean += om;
            sum_mean2 += om * om;
            sum_var += ov;
            sum_var2 += ov * ov;
        }
        ld mean_mean = sum_mean / N_bootstraps;
        ld mean_var = sum_var / N_bootstraps;
        ld std_mean = std::sqrtl((sum_mean2 / N_bootstraps) - (mean_mean * mean_mean));
        ld std_var = std::sqrtl((sum_var2 / N_bootstraps) - (mean_var * mean_var));
        O_means_rw[m] = (double) mean_mean;
        O_means_rw_err[m] = (double) std_mean;
        O_vars_rw[m] = (double) mean_var;
        O_vars_rw_err[m] = (double) std_var;
    }

    // Create output numpy arrays (each of shape (M,))
    auto result_O_means = py::array_t<double>(M);
    auto result_O_means_err = py::array_t<double>(M);
    auto result_O_vars = py::array_t<double>(M);
    auto result_O_vars_err = py::array_t<double>(M);

    py::buffer_info r1 = result_O_means.request();
    py::buffer_info r2 = result_O_means_err.request();
    py::buffer_info r3 = result_O_vars.request();
    py::buffer_info r4 = result_O_vars_err.request();

    double* r1_ptr = static_cast<double*>(r1.ptr);
    double* r2_ptr = static_cast<double*>(r2.ptr);
    double* r3_ptr = static_cast<double*>(r3.ptr);
    double* r4_ptr = static_cast<double*>(r4.ptr);
    for (int m = 0; m < M; m++) {
        r1_ptr[m] = O_means_rw[m];
        r2_ptr[m] = O_means_rw_err[m];
        r3_ptr[m] = O_vars_rw[m];
        r4_ptr[m] = O_vars_rw_err[m];
    }

    // Create a numpy array for the argmax beta values across bootstraps
    auto result_argmax = py::array_t<double>(N_bootstraps);
    py::buffer_info r5 = result_argmax.request();
    double* r5_ptr = static_cast<double*>(r5.ptr);
    for (int rep = 0; rep < N_bootstraps; rep++) {
        r5_ptr[rep] = argmax_beta_vals[rep];
    }

    // Return python tuple
    return py::make_tuple(result_O_means, result_O_means_err, result_O_vars, result_O_vars_err, result_argmax);
                }

// -------------------------
// Module definition
// -------------------------
PYBIND11_MODULE(reweight, m) {
m.doc() = "Multihistogram reweighting module implemented in C++ with pybind11";
m.def("reweight", &reweight,
      py::arg("O"), py::arg("E"), py::arg("betas"), py::arg("rw_betas"),
      py::arg("N_bootstraps") = 500, py::arg("tol") = 1e-14,
      "Perform multihistogram reweighting and return (O_means, O_means_err, O_vars, O_vars_err, O_vars_argmax_bs)");
}