#ifndef DYN_DIFF_HPP_INCLUDED
#define DYN_DIFF_HPP_INCLUDED

#include <boost/multi_array.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <math.h>

using namespace std;

/**
 * @brief  ODE system for Lotka-Voleterra dynamics in a coupled master-equation + mean-field system
 * @author  LHD
 * @since   2024-01-02
 */

struct Sparam
{
    const double beta;
    const double mu;
    const double K;
    const double nu;
    const int n_l;
    const int n_me1;
    const int n_me2;
}; // parameter structure

double truncated_Poisson_coupling(int k, const double average)
{

    if (average <= 0 || k <= 0)
        return 0.0;
    else
    {

        boost::math::poisson_distribution<> poiss(average);

        double answer = boost::math::pdf(poiss, k) / (1.0 - boost::math::cdf(poiss, k - 1));

        if (isinf(answer) || answer > 1.0)
            answer = 1.0;

        return answer;
    }
}

//********** function dydt definition **************************************************************
int dydt(double t, const double y[], double f[], void *param)
{
    // ODE system for diffusion

    // std::cout << t << std::endl;

    // Cast parameters
    Sparam &p = *static_cast<Sparam *>(param);

    // Create multi_array reference to y and f
    typedef boost::multi_array_ref<const double, 5> CSTmatref_type;
    typedef boost::multi_array_ref<double, 5> matref_type;
    typedef CSTmatref_type::index indexref;
    CSTmatref_type yref(y, boost::extents[p.n_l][p.n_l][1][p.n_me1][p.n_me2]);
    matref_type fref(f, boost::extents[p.n_l][p.n_l][1][p.n_me1][p.n_me2]);

    // Process all dynamics
    for (int d1 = 0; d1 < p.n_l; ++d1)
    {
        for (int d2 = 0; d2 < p.n_l; ++d2)
        {
            for (int d3 = 0; d3 < p.n_me1 - 1; ++d3)
            {
                for (int d4 = 0; d4 < p.n_me2 - 1; ++d4)
                {

                    if (d3 < p.n_me1 - 3 && d4 < p.n_me2 - 3)
                    { // general master equation
                        fref[d1][d2][0][d3][d4] = -p.mu * d3 * ((p.K - d3) / p.K) * yref[d1][d2][0][d3][d4] - p.nu * d4 * yref[d1][d2][0][d3][d4] - p.beta * d3 * d4 * yref[d1][d2][0][d3][d4];
                        fref[d1][d2][0][d3][d4] += p.nu * (d4 + 1) * yref[d1][d2][0][d3][d4 + 1];
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] += p.mu * (d3 - 1) * ((p.K - d3 + 1) / p.K) * yref[d1][d2][0][d3 - 1][d4];
                        if (d4 > 0)
                            fref[d1][d2][0][d3][d4] += p.beta * (d3 + 1) * (d4 - 1) * yref[d1][d2][0][d3 + 1][d4 - 1];
                    }

                    if (d3 == p.n_me1 - 3 && d4 < p.n_me2 - 3)
                    { // final states of master equation for prey
                        fref[d1][d2][0][d3][d4] = -p.mu * d3 * ((p.K - d3) / p.K) * yref[d1][d2][0][d3][d4] - p.nu * d4 * yref[d1][d2][0][d3][d4] - p.beta * d3 * d4 * yref[d1][d2][0][d3][d4];
                        fref[d1][d2][0][d3][d4] += p.nu * (d4 + 1) * yref[d1][d2][0][d3][d4 + 1];
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] += p.mu * (d3 - 1) * ((p.K - d3 + 1) / p.K) * yref[d1][d2][0][d3 - 1][d4];
                        if (d4 > 0)
                            fref[d1][d2][0][d3][d4] += p.beta * (d3 + 1) * truncated_Poisson_coupling(d3 + 1, yref[d1][d2][0][d3 + 2][d4 - 1]) * (d4 - 1) * yref[d1][d2][0][d3 + 1][d4 - 1];
                    }

                    if (d3 < p.n_me1 - 3 && d4 == p.n_me2 - 3)
                    { // final states of master equation for predator
                        fref[d1][d2][0][d3][d4] = -p.mu * d3 * ((p.K - d3) / p.K) * yref[d1][d2][0][d3][d4] - p.nu * d4 * yref[d1][d2][0][d3][d4] - p.beta * d3 * d4 * yref[d1][d2][0][d3][d4];
                        fref[d1][d2][0][d3][d4] += p.nu * (d4 + 1) * truncated_Poisson_coupling(d4 + 1, yref[d1][d2][0][d3][d4 + 2]) * yref[d1][d2][0][d3][d4 + 1];
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] += p.mu * (d3 - 1) * ((p.K - d3 + 1) / p.K) * yref[d1][d2][0][d3 - 1][d4];
                        if (d4 > 0)
                            fref[d1][d2][0][d3][d4] += p.beta * (d3 + 1) * (d4 - 1) * yref[d1][d2][0][d3 + 1][d4 - 1];
                    }

                    if (d3 == p.n_me1 - 3 && d4 == p.n_me2 - 3)
                    { // final states of master equations for prey + predator
                        fref[d1][d2][0][d3][d4] = -p.mu * d3 * ((p.K - d3) / p.K) * yref[d1][d2][0][d3][d4] - p.nu * d4 * yref[d1][d2][0][d3][d4] - p.beta * d3 * d4 * yref[d1][d2][0][d3][d4];
                        fref[d1][d2][0][d3][d4] += p.nu * (d4 + 1) * truncated_Poisson_coupling(d4 + 1, yref[d1][d2][0][d3][d4 + 2]) * yref[d1][d2][0][d3][d4 + 1];
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] += p.mu * (d3 - 1) * ((p.K - d3 + 1) / p.K) * yref[d1][d2][0][d3 - 1][d4];
                        if (d4 > 0)
                            fref[d1][d2][0][d3][d4] += p.beta * (d3 + 1) * truncated_Poisson_coupling(d3 + 1, yref[d1][d2][0][d3 + 2][d4 - 1]) * (d4 - 1) * yref[d1][d2][0][d3 + 1][d4 - 1];
                    }

                    if (d3 == p.n_me1 - 2 && d4 < p.n_me2 - 3)
                    { // prey mean-field equation, general predator master equation
                        fref[d1][d2][0][d3][d4] = -p.nu * d4 * yref[d1][d2][0][d3][d4];
                        fref[d1][d2][0][d3][d4] -= p.beta * ((1.0 - truncated_Poisson_coupling(d3, yref[d1][d2][0][d3 + 1][d4])) * yref[d1][d2][0][d3 + 1][d4] + truncated_Poisson_coupling(d3, yref[d1][d2][0][d3 + 1][d4]) * d3) * d4 * yref[d1][d2][0][d3][d4];
                        fref[d1][d2][0][d3][d4] += p.nu * (d4 + 1) * yref[d1][d2][0][d3][d4 + 1];
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] += p.mu * (d3 - 1) * ((p.K - d3 + 1) / p.K) * yref[d1][d2][0][d3 - 1][d4];
                        if (d4 > 0)
                            fref[d1][d2][0][d3][d4] += p.beta * yref[d1][d2][0][d3 + 1][d4 - 1] * (1.0 - truncated_Poisson_coupling(d3, yref[d1][d2][0][d3 + 1][d4 - 1])) * (d4 - 1) * yref[d1][d2][0][d3][d4 - 1];

                        fref[d1][d2][0][d3 + 1][d4] = p.mu * yref[d1][d2][0][d3 + 1][d4] * (p.K - yref[d1][d2][0][d3 + 1][d4]) / p.K - p.beta * d4 * yref[d1][d2][0][d3 + 1][d4]; // mean-field state
                    }

                    if (d3 == p.n_me1 - 2 && d4 == p.n_me2 - 3)
                    { // prey mean-field equation, final predator master equation state
                        fref[d1][d2][0][d3][d4] = -p.nu * d4 * yref[d1][d2][0][d3][d4];
                        fref[d1][d2][0][d3][d4] -= p.beta * ((1.0 - truncated_Poisson_coupling(d3, yref[d1][d2][0][d3 + 1][d4])) * yref[d1][d2][0][d3 + 1][d4] + truncated_Poisson_coupling(d3, yref[d1][d2][0][d3 + 1][d4]) * d3) * d4 * yref[d1][d2][0][d3][d4];
                        fref[d1][d2][0][d3][d4] += p.nu * (d4 + 1) * truncated_Poisson_coupling(d4 + 1, yref[d1][d2][0][d3][d4 + 2]) * yref[d1][d2][0][d3][d4 + 1];
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] += p.mu * (d3 - 1) * ((p.K - d3 + 1) / p.K) * yref[d1][d2][0][d3 - 1][d4];
                        if (d4 > 0)
                            fref[d1][d2][0][d3][d4] += p.beta * yref[d1][d2][0][d3 + 1][d4 - 1] * (1.0 - truncated_Poisson_coupling(d3, yref[d1][d2][0][d3 + 1][d4 - 1])) * (d4 - 1) * yref[d1][d2][0][d3][d4 - 1];

                        fref[d1][d2][0][d3 + 1][d4] = p.mu * yref[d1][d2][0][d3 + 1][d4] * (p.K - yref[d1][d2][0][d3 + 1][d4]) / p.K - p.beta * d4 * yref[d1][d2][0][d3 + 1][d4]; // mean-field state
                    }

                    if (d3 < p.n_me1 - 3 && d4 == p.n_me2 - 2)
                    { // general prey master equation, predator mean-field equation
                        fref[d1][d2][0][d3][d4] = -p.mu * d3 * ((p.K - d3) / p.K) * yref[d1][d2][0][d3][d4] - p.nu * d4 * truncated_Poisson_coupling(d4, yref[d1][d2][0][d3][d4 + 1]) * yref[d1][d2][0][d3][d4] - p.beta * d3 * yref[d1][d2][0][d3][d4 + 1] * yref[d1][d2][0][d3][d4];
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] += p.mu * (d3 - 1) * ((p.K - d3 + 1) / p.K) * yref[d1][d2][0][d3 - 1][d4];
                        if (d4 > 0)
                            fref[d1][d2][0][d3][d4] += p.beta * (d3 + 1) * (d4 - 1) * yref[d1][d2][0][d3 + 1][d4 - 1];
                        if (d4 > 0)
                            fref[d1][d2][0][d3][d4] += p.beta * (d3 + 1) * yref[d1][d2][0][d3 + 1][d4 + 1] * yref[d1][d2][0][d3 + 1][d4];

                        fref[d1][d2][0][d3][d4 + 1] = p.beta * d3 * yref[d1][d2][0][d3][d4 + 1] - p.nu * yref[d1][d2][0][d3][d4 + 1]; // mean-field state
                    }

                    if (d3 == p.n_me1 - 3 && d4 == p.n_me2 - 2)
                    { // final prey master equation state, predator mean-field equation
                        fref[d1][d2][0][d3][d4] = -p.mu * d3 * ((p.K - d3) / p.K) * yref[d1][d2][0][d3][d4] - p.nu * d4 * truncated_Poisson_coupling(d4, yref[d1][d2][0][d3][d4 + 1]) * yref[d1][d2][0][d3][d4] - p.beta * d3 * yref[d1][d2][0][d3][d4 + 1] * yref[d1][d2][0][d3][d4];
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] += p.mu * (d3 - 1) * ((p.K - d3 + 1) / p.K) * yref[d1][d2][0][d3 - 1][d4];
                        if (d4 > 0)
                            fref[d1][d2][0][d3][d4] += p.beta * (d3 + 1) * (d4 - 1) * truncated_Poisson_coupling(d3 + 1, yref[d1][d2][0][d3 + 2][d4 - 1]) * yref[d1][d2][0][d3 + 1][d4 - 1];
                        if (d4 > 0)
                            fref[d1][d2][0][d3][d4] += p.beta * (d3 + 1) * truncated_Poisson_coupling(d3 + 1, yref[d1][d2][0][d3 + 2][d4]) * yref[d1][d2][0][d3 + 1][d4 + 1] * yref[d1][d2][0][d3 + 1][d4];

                        fref[d1][d2][0][d3][d4 + 1] = p.beta * d3 * yref[d1][d2][0][d3][d4 + 1] - p.nu * yref[d1][d2][0][d3][d4 + 1]; // mean-field state
                    }

                    if (d3 == p.n_me1 - 2 && d4 == p.n_me2 - 2)
                    { // double mean-regime
                        fref[d1][d2][0][d3][d4] = -p.nu * d4 * truncated_Poisson_coupling(d4, yref[d1][d2][0][d3][d4 + 1]) * yref[d1][d2][0][d3][d4] - p.beta * d3 * truncated_Poisson_coupling(d3, yref[d1][d2][0][d3 + 1][d4]) * yref[d1][d2][0][d3][d4 + 1] * yref[d1][d2][0][d3][d4];
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] += p.mu * (d3 - 1) * ((p.K - d3 + 1) / p.K) * yref[d1][d2][0][d3 - 1][d4];
                        if (d4 > 0)
                            fref[d1][d2][0][d3][d4] += p.beta * (d4 - 1) * yref[d1][d2][0][d3 + 1][d4 - 1] * (1.0 - truncated_Poisson_coupling(d3, yref[d1][d2][0][d3 + 1][d4 - 1])) * yref[d1][d2][0][d3][d4 - 1];

                        fref[d1][d2][0][d3][d4 + 1] = p.beta * yref[d1][d2][0][d3 + 1][d4] * yref[d1][d2][0][d3][d4 + 1] - p.nu * yref[d1][d2][0][d3][d4 + 1];                                             // mean-field state
                        fref[d1][d2][0][d3 + 1][d4] = p.mu * yref[d1][d2][0][d3 + 1][d4] * (p.K - yref[d1][d2][0][d3 + 1][d4]) / p.K - p.beta * yref[d1][d2][0][d3][d4 + 1] * yref[d1][d2][0][d3 + 1][d4]; // mean-field state
                    }
                }
            }
        }
    }

    return GSL_SUCCESS;

} //********** end function dydt definition ********************************************************

#endif // DYN_DIFF_HPP_INCLUDED
