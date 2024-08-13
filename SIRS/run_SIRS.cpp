/**
 * @brief  ODE system for SIRS dynamics in a coupled master-equation + mean-field system
 * @author  LHD
 * @since   2023-10-16
 */

#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>
#include <sstream>

#include <boost/multi_array.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include "dyn_diff.hpp"

double binomial(int n, int k, double p);

using namespace std;

double Poisson(int k, const double average)
{

    if (average == 0)
        return 1.0;
    else
    {

        boost::math::poisson_distribution<> poiss(average);

        double answer = boost::math::pdf(poiss, k);

        return answer;
    }
}

int main(int argc, const char *argv[])
{
    if (argc < 6)
    {
        cerr << "Requires 6 parameters:\n"
             << "transmission rate (beta)\n"
             << "recovery rate (alpha)\n"
             << "waning rate (gamma)\n"
             << "side length of location grid (n_l)\n"
             << "number of master infectious states\n"
             << "number of master recovered states\n"
             << endl;
        return 0;
    }

    // Model parameters
    double beta = atof(argv[1]);  // transmission rate
    double alpha = atof(argv[2]); // recovery rate
    double gamma = atof(argv[3]); // waning rate
    int n_l = atoi(argv[4]);        // side length of location grid
    int n_me1 = atoi(argv[5]) + 2;  // number of states in master equation of dimension 1
    int n_me2 = atoi(argv[6]) + 2;  // number of states in master equation of dimension 2

    Sparam param = {beta, alpha, gamma, n_l, n_me1, n_me2};

    // Integrator parameters
    double t = 0.0;
    double dt = 0.1;
    double t_step = 1.0;
    const double eps_abs = 1e-4;
    const double eps_rel = 1e-5;

    // Setting initial conditions
    typedef boost::multi_array<double, 5> mat_type;
    typedef mat_type::index index;
    mat_type y(boost::extents[n_l][n_l][1][n_me1][n_me2]);
    fill(y.data(), y.data() + y.num_elements(), 0.0);

    // Initial conditions
    if (n_me1 == 2 && n_me2 == 2)
    {
        for (int d1 = 0; d1 < n_l; ++d1)
            for (int d2 = 0; d2 < n_l; ++d2)
            { // loop over all sites
                if (d1 < 5 && d2 < 5)
                {
                    y[d1][d2][0][0][0] = 1.0;
                    y[d1][d2][0][1][0] = 1.0;
                    y[d1][d2][0][0][1] = 0.0;
                }
                else
                {
                    y[d1][d2][0][0][0] = 1.0;
                    y[d1][d2][0][1][0] = 0.0;
                    y[d1][d2][0][0][1] = 0.0;
                }
            }
    }
    else
    {
        for (int d1 = 0; d1 < n_l; ++d1)
            for (int d2 = 0; d2 < n_l; ++d2)
            { // loop over all sites
                if (d1 < 5 && d2 < 5)
                    y[d1][d2][0][1][0] = 1.0;
                else
                    y[d1][d2][0][0][0] = 1.0;

                for (int d3 = 0; d3 < n_me1 - 1; ++d3)
                    y[d1][d2][0][d3][n_me2 - 1] = n_me2 - 1;
                for (int d4 = 0; d4 < n_me2 - 1; ++d4)
                    y[d1][d2][0][n_me1 - 1][d4] = n_me1 - 1;
            }
    }

    // Define GSL odeiv parameters
    long unsigned int sys_size = n_l * n_l * 1 * n_me1 * n_me2;
    const gsl_odeiv_step_type *step_type = gsl_odeiv_step_rkf45;
    gsl_odeiv_step *step = gsl_odeiv_step_alloc(step_type, sys_size);
    gsl_odeiv_control *control = gsl_odeiv_control_y_new(eps_abs, eps_rel);
    gsl_odeiv_evolve *evolve = gsl_odeiv_evolve_alloc(sys_size);
    gsl_odeiv_system sys = {dydt, NULL, sys_size, &param};

    // Integration
    int status(GSL_SUCCESS);
    for (double t_target = t + t_step; t_target < 10.0; t_target += t_step)
    { // stop by time

        while (t < t_target)
        {
            status = gsl_odeiv_evolve_apply(evolve, control, step, &sys, &t, t_target, &dt, y.data());
            if (status != GSL_SUCCESS)
            {
                cout << "ODE solver failed!" << endl;
                break;
            }
        } // end while

        // distribution output
        for (int d1 = 0; d1 < n_l; ++d1)
            for (int d2 = 0; d2 < n_l; ++d2)
            {
                for (int d3 = 0; d3 < n_me1 - 2; ++d3)
                {
                    double sum = 0.0;
                    for (int d4 = 0; d4 < n_me2 - 1; ++d4)
                        sum += y[d1][d2][0][d3][d4];
                    cout << "t=" << t
                         << ", d1=" << d1
                         << ", d2=" << d2
                         << ", d3=" << d3
                         << ", sum=" << sum << "\n";
                }
            }

    } // end while

    // final state output

    cout.flush();

    // Free memory
    gsl_odeiv_evolve_free(evolve);
    gsl_odeiv_control_free(control);
    gsl_odeiv_step_free(step);

    return 0;
}
