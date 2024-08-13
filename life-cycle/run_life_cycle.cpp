/**
 * @file   tevol_source_diff.hpp
 * @brief  ODE system for species dispersion in a coupled master-equation + mean-field system
 *
 * Source code. All parameters passed as arguments, but specify and compile to change precision or output format.
 * g++ -std=c++11 -O3 -o tevol_source_diff ./tevol_source_diff.cpp $(gsl-config --cflags) $(gsl-config --libs)
 *
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

int main(int argc, const char *argv[])
{
    if (argc < 11)
    {
        cerr << "Requires 11 parameters:\n"
             << "basic diffusion rate (beta)\n"
             << "rate of seed to seedling (seed2seedling)\n"
             << "seed death rate (seed_death)\n"
             << "seedling to sapling (seedling2sapling)\n"
             << "seedling death rate (seedling_death)\n"
             << "sapling to tree (sapling2tree)\n"
             << "sapling death rate (saplingdeath)\n"
             << "tree death rate (tree_death)\n"
             << "side length of location grid (n_l)\n"
             << "number of master sapling states\n"
             << "number of master adult states\n"
             << endl;
        return 0;
    }

    // Model parameters
    double beta = atof(argv[1]);             // basic diffusion rate
    double seed2seedling = atof(argv[2]);    // rate of seed 2 seedling transition
    double seed_death = atof(argv[3]);       // death rate for seeds
    double seedling2sapling = atof(argv[4]); // rate of seedling to sapling transition
    double seedling_death = atof(argv[5]);   // death rate for seedlings
    double sapling2tree = atof(argv[6]);     // rate of sapling to tree transition
    double sapling_death = atof(argv[7]);    // death rate for saplings
    double tree_death = atof(argv[8]);       // death rate for trees
    int n_l = atoi(argv[9]);                 // side length of location grid
    int n_me1 = atoi(argv[10]) + 2;          // number of states in master equation of dimension 1
    int n_me2 = atoi(argv[11]) + 2;          // number of states in master equation of dimension 2

    Sparam param = {
        beta,
        seed2seedling,
        seed_death,
        seedling2sapling,
        seedling_death,
        sapling2tree,
        sapling_death,
        tree_death,
        n_l,
        n_me1,
        n_me2};

    // Integrator parameters
    double t = 0.0;
    double dt = 1.0;
    double t_step = 0.1;
    const double eps_abs = 1e-1;
    const double eps_rel = 1e-3;

    // Setting initial conditions
    typedef boost::multi_array<double, 5> mat_type;
    typedef mat_type::index index;
    mat_type y(boost::extents[n_l][n_l][3][n_me1][n_me2]);
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
                    y[d1][d2][0][0][1] = 1.0;
                    y[d1][d2][1][0][0] = 10.0;
                    y[d1][d2][2][0][0] = 20.0;
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
                {
                    y[d1][d2][0][1][1] = 1.0;
                    y[d1][d2][1][1][1] = 10.0;
                    y[d1][d2][2][1][1] = 20.0;
                }
                else
                    y[d1][d2][0][0][0] = 1.0;
            }
    }

    // Define GSL odeiv parameters
    long unsigned int sys_size = n_l * n_l * 3 * n_me1 * n_me2;
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

        for (int d1 = 0; d1 < n_l; ++d1)
            for (int d2 = 0; d2 < n_l; ++d2)
            {
                double ext = 0.0;
                double avg = 0.0;
                double sum = 0.0;
                for (int d3 = 0; d3 < n_me1 - 1; ++d3)
                {
                    ext += y[d1][d2][0][d3][0];
                    for (int d4 = 0; d4 < n_me2 - 2; ++d4)
                    {
                        avg += 1.0 * d4 * y[d1][d2][0][d3][d4];
                        sum += y[d1][d2][0][d3][d4];
                    }
                    avg += y[d1][d2][0][d3][n_me2 - 2] * y[d1][d2][0][d3][n_me2 - 1];
                    sum += y[d1][d2][0][d3][n_me2 - 2];
                }
                cout << "t=" << t
                     << ", d1=" << d1
                     << ", d2=" << d2
                     << ", ext=" << ext
                     << ", avg=" << avg
                     << ", sum=" << sum
                     << endl;
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
