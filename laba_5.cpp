#include <iostream>
#include <cmath>
#include <functional>
#include <iomanip>
#include <stdexcept>


class NumericIntegrator {
public:
    // ---- Комбинированный метод Симпсона ----
    static double compositeSimpson(const std::function<double(double)>& f, double a, double b, int n) {
        if (n % 2 != 0) throw std::runtime_error("Simpson rule requires even n");

        double h = (b - a) / n;
        double sum = f(a) + f(b);

        for (int i = 1; i < n; ++i) {
            double x = a + i * h;
            sum += (i % 2 != 0) ? 4 * f(x) : 2 * f(x);
        }
        return sum * h / 3.0;
    }

    // ---- Комбинированная формула Гаусса (2 точки) ----
    static double compositeGauss2(const std::function<double(double)>& f, double a, double b, int n) {
        double h = (b - a) / n;
        const double t = 1.0 / std::sqrt(3.0);
        double sum = 0.0;

        for (int i = 0; i < n; i++) {
            double x0 = a + i * h;
            double x1 = x0 + h;
            double mid = (x0 + x1) / 2.0;
            double r = h / 2.0;

            double x_plus = mid + r * t;
            double x_minus = mid - r * t;

            sum += f(x_plus) + f(x_minus);
        }
        return sum * (h / 2.0);
    }
};

int main() {
    std::cout << std::scientific << std::setprecision(12);


    double a = 0.05, b = 0.56;
    auto phi = [](double x) { return std::exp(x); };


    double Itrue = std::exp(b) - std::exp(a);
    std::cout << "Analytical solution I* = e^b - e^a = " << Itrue << "\n\n";


    int n_h = 10;
    int n_h2 = 20;
    int n_h4 = 40;
    int n_h8 = 80;

    double h = (b - a) / n_h;
    double h2 = (b - a) / n_h2;
    double h4 = (b - a) / n_h4;
    double h8 = (b - a) / n_h8;

    std::cout << "Uniform grids:\n";
    std::cout << "h   = " << h << "   (n = " << n_h << ")\n";
    std::cout << "h/2 = " << h2 << " (n = " << n_h2 << ")\n";
    std::cout << "h/4 = " << h4 << " (n = " << n_h4 << ")\n";
    std::cout << "h/8 = " << h8 << " (n = " << n_h8 << ")\n\n";


    double S_h = NumericIntegrator::compositeSimpson(phi, a, b, n_h);
    double S_h2 = NumericIntegrator::compositeSimpson(phi, a, b, n_h2);
    double S_h4 = NumericIntegrator::compositeSimpson(phi, a, b, n_h4);
    double S_h8 = NumericIntegrator::compositeSimpson(phi, a, b, n_h8);

    double G_h = NumericIntegrator::compositeGauss2(phi, a, b, n_h);
    double G_h2 = NumericIntegrator::compositeGauss2(phi, a, b, n_h2);
    double G_h4 = NumericIntegrator::compositeGauss2(phi, a, b, n_h4);
    double G_h8 = NumericIntegrator::compositeGauss2(phi, a, b, n_h8);

    int k = 4;  // порядок формулы (Симпсон и Гаусс-2: O(h^4))


    auto print_table = [&](std::string method,
        double Ih, double Ih2, double Ih4, double Ih8,
        double h, double h2, double h4, double h8)
        {
            std::cout << std::scientific;
            double err_h = std::abs(Itrue - Ih);
            double err_h2 = std::abs(Itrue - Ih2);
            double p_est = std::log(err_h / err_h2) / std::log(2.0);
            double IR = Ih2 + (Ih2 - Ih) / (std::pow(2.0, k) - 1.0);

            std::cout << "============================\n";
            std::cout << method << "\n";
            std::cout << "h   = " << h << "\n";
            std::cout << "h/2 = " << h2 << "\n";
            std::cout << "h/4 = " << h4 << "\n";
            std::cout << "h/8 = " << h8 << "\n\n";

            std::cout << "I(h)   = " << Ih << "\n";
            std::cout << "I(h/2) = " << Ih2 << "\n";
            std::cout << "I(h/4) = " << Ih4 << "\n";
            std::cout << "I(h/8) = " << Ih8 << "\n\n";

            std::cout << "Error(h)   = " << err_h << "\n";
            std::cout << "Error(h/2) = " << err_h2 << "\n";
            std::cout << "Error(h/4) = " << std::abs(Itrue - Ih4) << "\n";
            std::cout << "Error(h/8) = " << std::abs(Itrue - Ih8) << "\n";

            std::cout << "Estimated order p ≈ log(err_h/err_h2)/log(2) = " << p_est << "\n";
            std::cout << "Richardson I^R = " << IR << "\n";
            std::cout << "============================\n\n";
        };


    print_table("Simpson composite", S_h, S_h2, S_h4, S_h8, h, h2, h4, h8);
    print_table("Gauss 2-point composite", G_h, G_h2, G_h4, G_h8, h, h2, h4, h8);


    auto print_detailed_table = [&](const std::string& method,
        double Ih, double Ih2, double Ih4, double Ih8,
        double Itrue, double h, int k)
        {
            std::cout << "\n" << method << " - detailed table:\n";
            std::cout << std::setw(10) << "h"
                << std::setw(20) << "I*-I^h"
                << std::setw(20) << "(I*-I^h)/(I*-I^{h/2})"
                << std::setw(20) << "(I^{h/2}-I^h)/(2^k-1)"
                << std::setw(20) << "I^R"
                << std::setw(20) << "I*-I^R"
                << "\n";

            auto calc_row = [&](double Ih_local, double Ih_next) {
                double diff = Itrue - Ih_local;
                double ratio = (Ih_next != 0) ? diff / (Itrue - Ih_next) : 0.0;
                double correction = (Ih_next != 0) ? (Ih_next - Ih_local) / (std::pow(2.0, k) - 1.0) : 0.0;
                double IR = Ih_next + correction;
                double IR_err = Itrue - IR;
                return std::make_tuple(diff, ratio, correction, IR, IR_err);
                };

            // строки для h, h/2, h/4, h/8
            double diff, ratio, correction, IR, IR_err;

            std::tie(diff, ratio, correction, IR, IR_err) = calc_row(Ih, Ih2);
            std::cout << std::setw(10) << h
                << std::setw(20) << diff
                << std::setw(20) << ratio
                << std::setw(20) << correction
                << std::setw(20) << IR
                << std::setw(20) << IR_err
                << "\n";

            std::tie(diff, ratio, correction, IR, IR_err) = calc_row(Ih2, Ih4);
            std::cout << std::setw(10) << h / 2
                << std::setw(20) << diff
                << std::setw(20) << ratio
                << std::setw(20) << correction
                << std::setw(20) << IR
                << std::setw(20) << IR_err
                << "\n";

            std::tie(diff, ratio, correction, IR, IR_err) = calc_row(Ih4, Ih8);
            std::cout << std::setw(10) << h / 4
                << std::setw(20) << diff
                << std::setw(20) << ratio
                << std::setw(20) << correction
                << std::setw(20) << IR
                << std::setw(20) << IR_err
                << "\n";

            // для h/8 выводим только I*-I^h, остальные значения не вычисляются
            diff = Itrue - Ih8;
            std::cout << std::setw(10) << h / 8
                << std::setw(20) << diff
                << std::setw(20) << " - "
                << std::setw(20) << " - "
                << std::setw(20) << " - "
                << std::setw(20) << " - "
                << "\n";
        };



    print_detailed_table("Simpson composite", S_h, S_h2, S_h4, S_h8, Itrue, h, k);
    print_detailed_table("Gauss 2-point composite", G_h, G_h2, G_h4, G_h8, Itrue, h, k);

    return 0;
}
