#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <filesystem>

struct MCStats {
    double mean{};
    double se{};
    double rel_se{};
    double ci_low{};
    double ci_high{};
};

struct ValidationRow {
    std::string section;
    std::string instrument;
    double value{};
    double se{};
    double rel_se{};
    double ci_low{};
    double ci_high{};
    double target_rel_se{};
    bool pass{};
    bool has_exact{};
    double exact{};
    double abs_diff{};
};

struct SensitivityRow {
    double a{};
    double sigma{};
    double rcc5{};
    double f5{};
    double zcb5_exact{};
    double zcb5_mc{};
    double zcb5_rel_se{};
    double short_rate_call{};
    double short_rate_call_rel_se{};
    double bond_call{};
    double bond_call_rel_se{};
};

static bool approx_equal(double x, double y, double eps = 1e-12) {
    return std::fabs(x - y) <= eps;
}

static MCStats mean_and_se(const std::vector<double>& x) {
    if (x.empty()) throw std::runtime_error("Empty sample");
    const double n = static_cast<double>(x.size());
    double sum = 0.0;
    for (double v : x) sum += v;
    const double mean = sum / n;
    double ss = 0.0;
    for (double v : x) {
        const double d = v - mean;
        ss += d * d;
    }
    const double var = (x.size() > 1) ? (ss / static_cast<double>(x.size() - 1)) : 0.0;
    const double se = std::sqrt(var / n);
    const double rel_se = (std::fabs(mean) > 0.0) ? (se / std::fabs(mean)) : std::numeric_limits<double>::infinity();
    const double z = 1.96;

    return {mean, se, rel_se, mean - z * se, mean + z * se};
}

class DiscountCurve {
    private:
        std::map<double, double> points_;

    public:
        DiscountCurve() {points_[0.0] = 1.0;}

        void add_point(double t, double df) {
            if (t < 0.0) throw std::runtime_error("Negative maturity");
            if (!(df > 0.0 && df <= 1.0 + 1e-10)) throw std::runtime_error("Discount factor out of range");
            points_[t] = df;
        }

        double discount_factor(double t) const {
            if (t < 0.0) throw std::runtime_error("Negative maturity");
            const auto it = points_.find(t);
            if (it != points_.end()) return it->second;
            const auto upper = points_.lower_bound(t);
            if (upper == points_.begin() || upper == points_.end()) {
                throw std::runtime_error("Interpolation point out of range");
            }

            const auto lower = std::prev(upper);
            const double t0 = lower->first;
            const double t1 = upper->first;
            const double d0 = lower->second;
            const double d1 = upper->second;
            const double w = (t - t0) / (t1 - t0);
            return d0 + w * (d1 - d0);
        }

        double zero_rate_cc(double t) const {
            if (t <= 0.0) return 0.0;
            return -std::log(discount_factor(t)) / t;
        }

        double simple_forward_rate(double t0, double t1) const {
            if (!(t1 > t0)) throw std::runtime_error("Require t1 > t0");
            const double d0 = discount_factor(t0);
            const double d1 = discount_factor(t1);
            return (d0 / d1 - 1.0) / (t1 - t0);
        }

        const std::map<double, double>& points() const {
            return points_;
        }
};

class Swap {
    private:
        double maturity_{};
        double rate_{};
        int f_{};

    public:
        Swap(double maturity_years, double par_rate, int cashflows_per_year) : maturity_(maturity_years), rate_(par_rate), f_(cashflows_per_year) {
            if (maturity_ <= 0.0) throw std::runtime_error("Swap maturity must be positive");
            if (f_ <= 0) throw std::runtime_error("Cashflow frequency must be positive");
        }

        std::vector<double> cashflow_times() const {
            const int n = static_cast<int>(std::round(maturity_ * static_cast<double>(f_)));
            std::vector<double> times;
            times.reserve(static_cast<std::size_t>(n));
            for (int i = 1; i <= n; ++i) {
                times.push_back(static_cast<double>(i) / static_cast<double>(f_));
            }
            return times;
        }

        double price(const DiscountCurve& curve) const {
            double fixed_leg = 0.0;
            for (double t : cashflow_times()) {
                fixed_leg += curve.discount_factor(t);
            }
            fixed_leg *= (rate_ / static_cast<double>(f_));
            return fixed_leg + curve.discount_factor(maturity_) - 1.0;
        }

        double par_rate_from_curve(const DiscountCurve& curve) const {
            double annuity = 0.0;
            for (double t : cashflow_times()) {
                annuity += curve.discount_factor(t);
            }
            return static_cast<double>(f_) * (1.0 - curve.discount_factor(maturity_)) / annuity;
        }
};

class Bootstrapper {
    private:
        int f_{};

    public:
        explicit Bootstrapper(int cashflows_per_year) : f_(cashflows_per_year) {
            if (f_ <= 0) throw std::runtime_error("Cashflow frequency must be positive");
        }

        DiscountCurve bootstrap(const std::vector<std::pair<double, double>>& market_swaps) const {
            if (market_swaps.empty()) throw std::runtime_error("No market swap data provided");
            double prev_T = 0.0;
            for (const auto& [T, _] : market_swaps) {
                if (T <= prev_T) throw std::runtime_error("Swap maturities must be strictly increasing");
                prev_T = T;
            }

            DiscountCurve curve;
            prev_T = 0.0;
            double prev_D = 1.0;

            for (const auto& [T, S] : market_swaps) {
                const int n = static_cast<int>(std::round(T * static_cast<double>(f_)));
                double constant_part = -1.0;
                double coeff_on_DT = 1.0;

                for (int i = 1; i <= n; ++i) {
                    const double t = static_cast<double>(i) / static_cast<double>(f_);
                    double alpha = 0.0;
                    double beta = 0.0;
                    if (approx_equal(t, T)) {
                        beta = 1.0;
                    } else if (t <= prev_T + 1e-12) {
                        alpha = curve.discount_factor(t);
                    } else {
                        const double w = (t - prev_T) / (T - prev_T);
                        alpha = (1.0 - w) * prev_D;
                        beta = w;
                    }
                    constant_part += (S / static_cast<double>(f_)) * alpha;
                    coeff_on_DT += (S / static_cast<double>(f_)) * beta;
                }
                const double DT = -constant_part / coeff_on_DT;
                curve.add_point(T, DT);
                prev_T = T;
                prev_D = DT;
            }
            return curve;
        }
};

class Vasicek {
    private:
        double a_{};
        double b_{};
        double sigma_{};
        double r0_{};

    public:
        Vasicek(double a, double b, double sigma, double r0) : a_(a), b_(b), sigma_(sigma), r0_(r0) {
            if (a_ <= 0.0) throw std::runtime_error("a must be positive");
            if (sigma_ < 0.0) throw std::runtime_error("sigma must be non-negative");
        }
        double a() const { return a_; }
        double b() const { return b_; }
        double sigma() const { return sigma_; }
        double r0() const { return r0_; }

        double B(double T) const {
            if (T == 0.0) return 0.0;
            return (1.0 - std::exp(-a_ * T)) / a_;
        }

        double A(double T) const {
            const double Bv = B(T);
            const double term1 = (Bv - T) * ((a_ * a_ * b_ - 0.5 * sigma_ * sigma_) / (a_ * a_));
            const double term2 = (sigma_ * sigma_ * Bv * Bv) / (4.0 * a_);
            return term1 - term2;
        }

        double discount_factor(double T, double r = std::numeric_limits<double>::quiet_NaN()) const {
            const double r_use = std::isnan(r) ? r0_ : r;
            return std::exp(A(T) - B(T) * r_use);
        }

        double spot_rate_cc(double T, double r = std::numeric_limits<double>::quiet_NaN()) const {
            const double r_use = std::isnan(r) ? r0_ : r;
            if (T == 0.0) return r_use;
            return -std::log(discount_factor(T, r_use)) / T;
        }

        double forward_rate_cc(double T, double r = std::numeric_limits<double>::quiet_NaN()) const {
            const double r_use = std::isnan(r) ? r0_ : r;
            const double exp_at = std::exp(-a_ * T);
            const double Bv = B(T);
            const double Bp = exp_at;
            const double Ap = ((Bp - 1.0) * ((a_ * a_ * b_ - 0.5 * sigma_ * sigma_) / (a_ * a_))) - (sigma_ * sigma_ * Bv * Bp) / (2.0 * a_);
            return -(Ap - Bp * r_use);
        }
};

struct SimResult {
    std::vector<double> t;
    std::vector<std::vector<double>> r;
};

static SimResult simulate_short_rate_paths_exact(const Vasicek& model, double T, double dt_target, std::size_t n_paths, std::uint64_t seed) {
    if (T <= 0.0) throw std::runtime_error("Simulation horizon must be positive");
    if (dt_target <= 0.0) throw std::runtime_error("dt must be positive");
    if (n_paths == 0) throw std::runtime_error("n_paths must be positive");

    const std::size_t steps = static_cast<std::size_t>(std::ceil(T / dt_target));
    const double dt = T / static_cast<double>(steps);
    const double m = std::exp(-model.a() * dt);
    const double s = model.sigma() * std::sqrt((1.0 - std::exp(-2.0 * model.a() * dt)) / (2.0 * model.a()));

    std::mt19937_64 rng(seed);
    std::normal_distribution<double> nd(0.0, 1.0);

    SimResult out;
    out.t.resize(steps + 1);
    for (std::size_t k = 0; k <= steps; ++k) {
        out.t[k] = T * static_cast<double>(k) / static_cast<double>(steps);
    }
    out.r.assign(n_paths, std::vector<double>(steps + 1, 0.0));
    for (std::size_t i = 0; i < n_paths; ++i) out.r[i][0] = model.r0();
    for (std::size_t k = 0; k < steps; ++k) {
        for (std::size_t i = 0; i < n_paths; ++i) {
            const double z = nd(rng);
            out.r[i][k + 1] = model.b() + (out.r[i][k] - model.b()) * m + s * z;
        }
    }

    return out;
}

static std::vector<double> path_discount_factors_trapezoid(const std::vector<std::vector<double>>& r_paths, const std::vector<double>& t_grid) {
    const std::size_t n_paths = r_paths.size();
    const std::size_t steps = t_grid.size() - 1;
    std::vector<double> disc(n_paths, 1.0);

    for (std::size_t i = 0; i < n_paths; ++i) {
        double integral = 0.0;
        for (std::size_t k = 0; k < steps; ++k) {
            const double dt = t_grid[k + 1] - t_grid[k];
            integral += 0.5 * (r_paths[i][k] + r_paths[i][k + 1]) * dt;
        }
        disc[i] = std::exp(-integral);
    }

    return disc;
}

static void write_histogram_csv(const std::vector<double>& samples, std::size_t bins, const std::string& filename) {
    if (samples.empty() || bins == 0) return;
    const auto [mn_it, mx_it] = std::minmax_element(samples.begin(), samples.end());
    double mn = *mn_it;
    double mx = *mx_it;
    if (mx == mn) mx = mn + 1e-12;
    const double width = (mx - mn) / static_cast<double>(bins);
    std::vector<std::size_t> counts(bins, 0);

    for (double x : samples) {
        std::size_t idx = static_cast<std::size_t>((x - mn) / width);
        if (idx >= bins) idx = bins - 1;
        counts[idx]++;
    }

    std::ofstream out(filename);
    if (!out) throw std::runtime_error("Failed to open: " + filename);
    out << "bin_left,bin_right,count\n";
    for (std::size_t i = 0; i < bins; ++i) {
        const double left = mn + width * static_cast<double>(i);
        const double right = left + width;
        out << std::setprecision(12) << left << "," << right << "," << counts[i] << "\n";
    }
}

struct ProjectConfig {
    int swap_frequency = 2;
    std::vector<std::pair<double, double>> market_swaps = {
        {1.0, 0.0264},
        {2.0, 0.0302},
        {3.0, 0.0342},
        {5.0, 0.0411},
        {7.0, 0.0456},
        {10.0, 0.0497}
    };

    double seven_year_test_rate = 0.0442;
    double q1_curve_grid_step = 0.5;
    double a = 0.25;
    double b = 0.10;
    double sigma = 0.02;
    double r0 = 0.07;

    double dt = 1.0 / 252.0;
    std::uint64_t seed = 1234;

    std::size_t n_paths_hist = 30000;        
    std::size_t n_paths_bond = 30000;     
    std::size_t n_paths_opt = 30000;         
    std::size_t n_paths_sensitivity = 10000; 
    std::size_t hist_bins = 50;          

    double q2_curve_T_max = 15.0;
    std::size_t q2_curve_grid_n = 600;
    std::vector<int> histogram_years = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

    double bond_T = 5.0;
    double short_rate_call_K = 0.11;
    double short_rate_call_M = 1000000.0;
    double short_rate_call_T = 5.0;
    double bond_call_T1 = 1.0;
    double bond_call_T2 = 5.0;
    double bond_call_K = 0.69;

    double bond_target_rel_se = 0.001;
    double option_target_rel_se = 0.01;

    std::vector<double> sensitivity_a = {0.10, 0.25, 0.50};
    std::vector<double> sensitivity_sigma = {0.01, 0.02, 0.04};
};

static void write_q1_curve_csv(const DiscountCurve& curve, double step, double T_max, const std::string& filename) {
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("Failed to open: " + filename);
    out << "T,discount_factor,zero_rate_cc,forward_6m_simple\n";
    for (double t = 0.0; t <= T_max + 1e-12; t += step) {
        const double df = curve.discount_factor(t);
        const double zr = curve.zero_rate_cc(t);
        double fwd = std::numeric_limits<double>::quiet_NaN();
        if (t >= step && t <= T_max - step + 1e-12) {
            fwd = curve.simple_forward_rate(t, t + step);
        }
        out << std::setprecision(12) << t << "," << df << "," << zr << ",";
        if (std::isnan(fwd)) out << "";
        else out << fwd;
        out << "\n";
    }
}

static void write_q1_summary_csv(double zcb_1y_price, double deposit_rate_2y_cc, double seven_year_swap_price, double nine_year_par_rate, const std::string& filename) {
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("Failed to open: " + filename);
    out << "metric,value\n";
    out << "zcb_1y_price," << std::setprecision(12) << zcb_1y_price << "\n";
    out << "deposit_rate_2y_cc," << deposit_rate_2y_cc << "\n";
    out << "swap_7y_price_4.42pct," << seven_year_swap_price << "\n";
    out << "swap_9y_par_rate," << nine_year_par_rate << "\n";
}

static void write_q2_term_structure_csv(const Vasicek& model, double T_max, std::size_t n_grid, const std::string& filename) {
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("Failed to open: " + filename);
    out << "T,discount_factor,spot_rate_cc,forward_rate_cc\n";
    for (std::size_t i = 0; i < n_grid; ++i) {
        const double T = T_max * static_cast<double>(i) / static_cast<double>(n_grid - 1);
        out << std::setprecision(12) << T << ","
            << model.discount_factor(T) << ","
            << model.spot_rate_cc(T) << ","
            << model.forward_rate_cc(T) << "\n";
    }
}

static void write_hist_summary_csv(const std::vector<std::tuple<int, double, double>>& rows, const std::string& filename) {
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("Failed to open: " + filename);
    out << "year,mean_short_rate,se\n";
    for (const auto& [year, mean, se] : rows) {
        out << year << "," << std::setprecision(12) << mean << "," << se << "\n";
    }
}

static void write_validation_csv(const std::vector<ValidationRow>& rows, const std::string& filename) {
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("Failed to open: " + filename);
    out << "section,instrument,value,se,rel_se,ci_low,ci_high,target_rel_se,pass,exact,abs_diff\n";
    for (const auto& row : rows) {
        out << row.section << ","
            << row.instrument << ","
            << std::setprecision(12) << row.value << ","
            << row.se << ","
            << row.rel_se << ","
            << row.ci_low << ","
            << row.ci_high << ","
            << row.target_rel_se << ","
            << (row.pass ? 1 : 0) << ",";
        if (row.has_exact) out << row.exact; 
        out << ",";
        if (row.has_exact) out << row.abs_diff;
        out << "\n";
    }
}

static void write_sensitivity_csv(const std::vector<SensitivityRow>& rows, const std::string& filename) {
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("Failed to open: " + filename);
    out << "a,sigma,rcc5,f5,zcb5_exact,zcb5_mc,zcb5_rel_se,short_rate_call,short_rate_call_rel_se,bond_call,bond_call_rel_se\n";
    for (const auto& row : rows) {
        out << std::setprecision(12)
            << row.a << ","
            << row.sigma << ","
            << row.rcc5 << ","
            << row.f5 << ","
            << row.zcb5_exact << ","
            << row.zcb5_mc << ","
            << row.zcb5_rel_se << ","
            << row.short_rate_call << ","
            << row.short_rate_call_rel_se << ","
            << row.bond_call << ","
            << row.bond_call_rel_se << "\n";
    }
}

static ValidationRow price_zcb_mc(const Vasicek& model, double T, double dt, std::size_t n_paths, std::uint64_t seed, double target_rel_se) {
    const auto sim = simulate_short_rate_paths_exact(model, T, dt, n_paths, seed);
    const auto disc = path_discount_factors_trapezoid(sim.r, sim.t);
    const auto stats = mean_and_se(disc);
    const double exact = model.discount_factor(T);

    return {
        "Q2(c)",
        "5Y zero coupon bond",
        stats.mean,
        stats.se,
        stats.rel_se,
        stats.ci_low,
        stats.ci_high,
        target_rel_se,
        stats.rel_se < target_rel_se,
        true,
        exact,
        std::fabs(stats.mean - exact)
    };
}

static ValidationRow price_short_rate_call_mc(const Vasicek& model, double K, double M, double T, double dt, std::size_t n_paths, std::uint64_t seed, double target_rel_se) {
    const auto sim = simulate_short_rate_paths_exact(model, T, dt, n_paths, seed);
    const auto disc = path_discount_factors_trapezoid(sim.r, sim.t);
    std::vector<double> pv(n_paths, 0.0);
    for (std::size_t i = 0; i < n_paths; ++i) {
        const double payoff = M * std::max(sim.r[i].back() - K, 0.0);
        pv[i] = disc[i] * payoff;
    }
    const auto stats = mean_and_se(pv);
    return {
        "Q2(d)",
        "Short-rate call",
        stats.mean,
        stats.se,
        stats.rel_se,
        stats.ci_low,
        stats.ci_high,
        target_rel_se,
        stats.rel_se < target_rel_se,
        false,
        0.0,
        0.0
    };
}

static ValidationRow price_bond_call_mc(const Vasicek& model, double T1, double T2, double K, double dt, std::size_t n_paths, std::uint64_t seed, double target_rel_se) {
    if (!(T2 > T1)) throw std::runtime_error("Require T2 > T1.");
    const auto sim = simulate_short_rate_paths_exact(model, T1, dt, n_paths, seed);
    const auto disc_0_T1 = path_discount_factors_trapezoid(sim.r, sim.t);
    const double tau = T2 - T1;
    std::vector<double> pv(n_paths, 0.0);
    for (std::size_t i = 0; i < n_paths; ++i) {
        const double P_T1_T2 = std::exp(model.A(tau) - model.B(tau) * sim.r[i].back());
        const double payoff = std::max(P_T1_T2 - K, 0.0);
        pv[i] = disc_0_T1[i] * payoff;
    }
    const auto stats = mean_and_se(pv);
    return {
        "Q2(e)",
        "Bond call",
        stats.mean,
        stats.se,
        stats.rel_se,
        stats.ci_low,
        stats.ci_high,
        target_rel_se,
        stats.rel_se < target_rel_se,
        false,
        0.0,
        0.0
    };
}

static std::vector<std::tuple<int, double, double>> run_histograms(const Vasicek& model, const ProjectConfig& cfg, const std::string& filename_prefix){
    std::vector<std::tuple<int, double, double>> rows;
    rows.reserve(cfg.histogram_years.size());
    for (int year : cfg.histogram_years) {
        const auto sim = simulate_short_rate_paths_exact(model, static_cast<double>(year), cfg.dt, cfg.n_paths_hist, cfg.seed + static_cast<std::uint64_t>(year));
        std::vector<double> rT(cfg.n_paths_hist);
        for (std::size_t i = 0; i < cfg.n_paths_hist; ++i) rT[i] = sim.r[i].back();
        const auto stats = mean_and_se(rT);
        rows.emplace_back(year, stats.mean, stats.se);
        const std::string filename = filename_prefix + std::to_string(year) + "y.csv";
        write_histogram_csv(rT, cfg.hist_bins, filename);
    }
    return rows;
}

static void print_validation_row(const ValidationRow& row) {
    std::cout << row.section << " - " << row.instrument << "\n";
    std::cout << "  value      = " << std::setprecision(12) << row.value << "\n";
    std::cout << "  SE         = " << row.se << "\n";
    std::cout << "  rel SE     = " << (100.0 * row.rel_se) << "%\n";
    std::cout << "  95% CI     = [" << row.ci_low << ", " << row.ci_high << "]\n";
    std::cout << "  target     = " << (100.0 * row.target_rel_se) << "%\n";
    std::cout << "  pass       = " << (row.pass ? "YES" : "NO") << "\n";
    if (row.has_exact) {
        std::cout << "  exact      = " << row.exact << "\n";
        std::cout << "  abs diff   = " << row.abs_diff << "\n";
    }
    std::cout << "\n";
}

int main() {
    try {
        std::cout << "Working directory: ";     
        std::cout << std::setprecision(12);
        std::cout << "FIXED INCOME VASICEK MODELLING & MONTE CARLO PRICING\n\n";
        const ProjectConfig cfg;
        const Bootstrapper bootstrapper(cfg.swap_frequency);
        const DiscountCurve curve = bootstrapper.bootstrap(cfg.market_swaps);

        const double zcb_1y_price = curve.discount_factor(1.0);
        const double deposit_rate_2y_cc = curve.zero_rate_cc(2.0);
        const Swap seven_year_swap(7.0, cfg.seven_year_test_rate, cfg.swap_frequency);
        const double seven_year_swap_price = seven_year_swap.price(curve);
        const Swap nine_year_swap(9.0, 0.0, cfg.swap_frequency);
        const double nine_year_par_rate = nine_year_swap.par_rate_from_curve(curve);

        write_q1_curve_csv(curve, cfg.q1_curve_grid_step, 10.0, "q1_discount_curve.csv");
        write_q1_summary_csv(zcb_1y_price, deposit_rate_2y_cc, seven_year_swap_price, nine_year_par_rate, "q1_summary.csv");

        std::cout << "Q1 summary\n";
        std::cout << "  1Y ZCB price            = " << zcb_1y_price << "\n";
        std::cout << "  2Y zero rate (cc)       = " << deposit_rate_2y_cc << "  (" << 100.0 * deposit_rate_2y_cc << "%)\n";
        std::cout << "  7Y swap price @ 4.42%   = " << seven_year_swap_price << "\n";
        std::cout << "  9Y par swap rate        = " << nine_year_par_rate << "  (" << 100.0 * nine_year_par_rate << "%)\n";

        const Vasicek base_model(cfg.a, cfg.b, cfg.sigma, cfg.r0);
        write_q2_term_structure_csv(base_model, cfg.q2_curve_T_max, cfg.q2_curve_grid_n, "q2_term_structure.csv");

        std::cout << "Q2(a)\n";
        std::cout << "  D(5)   = " << base_model.discount_factor(5.0) << "\n";
        std::cout << "  Rcc(5) = " << base_model.spot_rate_cc(5.0) << "\n";
        std::cout << "  f(5)   = " << base_model.forward_rate_cc(5.0) << "\n";

        const auto hist_rows = run_histograms(base_model, cfg, "q2_hist_");
        write_hist_summary_csv(hist_rows, "q2_hist_summary.csv");

        std::cout << "Q2(b)\n";
        for (const auto& [year, mean, se] : hist_rows) {
            std::cout << "  year=" << year << "  mean(r_t)=" << mean << "  SE=" << se << "\n";
        }
        std::cout << "  wrote q2_hist_summary.csv and q2_hist_{1y-15y}.csv\n\n";

        std::vector<ValidationRow> validation_rows;
        validation_rows.push_back(price_zcb_mc(base_model, cfg.bond_T, cfg.dt, cfg.n_paths_bond, cfg.seed, cfg.bond_target_rel_se));
        validation_rows.push_back(price_short_rate_call_mc(base_model, cfg.short_rate_call_K, cfg.short_rate_call_M, cfg.short_rate_call_T, cfg.dt, cfg.n_paths_opt, cfg.seed + 101, cfg.option_target_rel_se));
        validation_rows.push_back(price_bond_call_mc(base_model, cfg.bond_call_T1, cfg.bond_call_T2, cfg.bond_call_K, cfg.dt, cfg.n_paths_opt, cfg.seed + 202, cfg.option_target_rel_se));
        write_validation_csv(validation_rows, "q2_validation.csv");

        std::cout << "Validation summary\n";
        for (const auto& row : validation_rows) {
            print_validation_row(row);
        }
        std::cout << "  wrote q2_validation.csv\n\n";

        std::vector<SensitivityRow> sensitivity_rows;
        for (double a : cfg.sensitivity_a) {
            for (double sigma : cfg.sensitivity_sigma) {
                const Vasicek model(a, cfg.b, sigma, cfg.r0);
                const auto zcb_row = price_zcb_mc(model, cfg.bond_T, cfg.dt, cfg.n_paths_sensitivity, cfg.seed + static_cast<std::uint64_t>(1000 * a + 10000 * sigma), cfg.bond_target_rel_se);
                const auto sr_row = price_short_rate_call_mc(model, cfg.short_rate_call_K, cfg.short_rate_call_M, cfg.short_rate_call_T, cfg.dt, cfg.n_paths_sensitivity, cfg.seed + static_cast<std::uint64_t>(2000 * a + 20000 * sigma), cfg.option_target_rel_se);
                const auto bc_row = price_bond_call_mc(model, cfg.bond_call_T1, cfg.bond_call_T2, cfg.bond_call_K, cfg.dt, cfg.n_paths_sensitivity, cfg.seed + static_cast<std::uint64_t>(3000 * a + 30000 * sigma), cfg.option_target_rel_se);
                sensitivity_rows.push_back({a, sigma, model.spot_rate_cc(5.0), model.forward_rate_cc(5.0), model.discount_factor(5.0), zcb_row.value, zcb_row.rel_se, sr_row.value, sr_row.rel_se, bc_row.value, bc_row.rel_se});
            }
        }

        write_sensitivity_csv(sensitivity_rows, "q2_sensitivity.csv");
        std::cout << "Sensitivity analysis\n";
        std::cout << "  scenarios = " << sensitivity_rows.size() << "\n";
        std::cout << "  wrote q2_sensitivity.csv\n\n";
        std::cout << "Run complete.\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }
}