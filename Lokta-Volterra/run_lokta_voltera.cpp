/**
 * @file   tevol_source_diff.hpp
 * @brief  ODE system for Lotka-Volterra dynamics in a coupled master-equation + mean-field system
 * @author  LHD
 * @since   2024-01-02
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

int main(int argc, const char *argv[])
{
    if (argc < 7)
    {
        cerr << "Requires 7 parameters:\n"
             << "transmission rate (beta)\n"
             << "reproduction rate (mu)\n"
             << "carrying capacity (K)\n"
             << "death rate (nu)\n"
             << "side length of location grid (n_l)\n"
             << "number of master infectious states\n"
             << "number of master recovered states\n"
             << endl;
        return 0;
    }

    cout << "Starting simulation..." << endl;
    
    // Model parameters
    double beta = atof(argv[1]); // transmission rate
    double mu = atof(argv[2]);   // reproduction rate
    double K = atof(argv[3]);    // carrying capacity
    double nu = atof(argv[4]);   // death rate
    int n_l = atoi(argv[5]);       // side length of location grid
    int n_me1 = atoi(argv[6]) + 2; // number of states in master equation of dimension 1
    int n_me2 = atoi(argv[7]) + 2; // number of states in master equation of dimension 2

    Sparam param = {beta, mu, K, nu, n_l, n_me1, n_me2};

    // Integrator parameters
    double t = 0.0;
    double dt = 0.1;
    double t_step = 0.25;
    const double eps_abs = 1e-1;
    const double eps_rel = 1e-2;

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
                    y[d1][d2][0][1][0] = 5.0;
                    y[d1][d2][0][0][1] = 5.0;
                }
                else
                {
                    y[d1][d2][0][0][0] = 1.0;
                    y[d1][d2][0][1][0] = 5.0;
                    y[d1][d2][0][0][1] = 5.0;
                }
            }
    }
    else
    {
        for (int d1 = 0; d1 < n_l; ++d1)
            for (int d2 = 0; d2 < n_l; ++d2)
            { // loop over all sites
                if (d1 < 5 && d2 < 5)
                    y[d1][d2][0][1][1] = 1.0;
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
    for (double t_target = t + t_step; t_target < 100.0; t_target += t_step)
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

        // extinction timeseries and average mature output
        cout << "t=" << t
             << ", y=" << y[0][0][0][0][0]
             << endl;

    } // end while

    // final state output

    cout.flush();

    // Free memory
    gsl_odeiv_evolve_free(evolve);
    gsl_odeiv_control_free(control);
    gsl_odeiv_step_free(step);

    return 0;
}
