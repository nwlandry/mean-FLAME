#ifndef DYN_DIFF_HPP_INCLUDED
#define DYN_DIFF_HPP_INCLUDED

#include <boost/multi_array.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <math.h>

using namespace std;

/**
 * @file   dyn_diff.hpp
 * @brief  ODE system for Lotka-Voleterra dynamics in a coupled master-equation + mean-field system
 *
 * @author  LHD
 * @since   2024-01-02
 */

struct Sparam
{
    const double beta;
    const double mu;
    const double K;
    const double nu;
    const int N;
    const int MF3;
    const int MF4;
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
    CSTmatref_type yref(y, boost::extents[p.N][p.N][1][p.MF3][p.MF4]);
    // fill(yref.data(),yref.data()+yref.num_elements(),0.0);
    matref_type fref(f, boost::extents[p.N][p.N][1][p.MF3][p.MF4]);
    // fill(fref.data(),fref.data()+fref.num_elements(),0.0);

    // Process all dynamics
    for (int d01 = 0; d01 < p.N; ++d01)
    {
        for (int d02 = 0; d02 < p.N; ++d02)
        {
            for (int d3 = 0; d3 < p.MF3 - 1; ++d3)
            {
                for (int d4 = 0; d4 < p.MF4 - 1; ++d4)
                {

                    if (d3 < p.MF3 - 3 && d4 < p.MF4 - 3)
                    { // general master equation
                        fref[d01][d02][0][d3][d4] = -p.mu * d3 * ((p.K - d3) / p.K) * yref[d01][d02][0][d3][d4] - p.nu * d4 * yref[d01][d02][0][d3][d4] - p.beta * d3 * d4 * yref[d01][d02][0][d3][d4];
                        fref[d01][d02][0][d3][d4] += p.nu * (d4 + 1) * yref[d01][d02][0][d3][d4 + 1];
                        if (d3 > 0)
                            fref[d01][d02][0][d3][d4] += p.mu * (d3 - 1) * ((p.K - d3 + 1) / p.K) * yref[d01][d02][0][d3 - 1][d4];
                        if (d4 > 0)
                            fref[d01][d02][0][d3][d4] += p.beta * (d3 + 1) * (d4 - 1) * yref[d01][d02][0][d3 + 1][d4 - 1];
                    }

                    if (d3 == p.MF3 - 3 && d4 < p.MF4 - 3)
                    { // final states of master equation for prey
                        fref[d01][d02][0][d3][d4] = -p.mu * d3 * ((p.K - d3) / p.K) * yref[d01][d02][0][d3][d4] - p.nu * d4 * yref[d01][d02][0][d3][d4] - p.beta * d3 * d4 * yref[d01][d02][0][d3][d4];
                        fref[d01][d02][0][d3][d4] += p.nu * (d4 + 1) * yref[d01][d02][0][d3][d4 + 1];
                        if (d3 > 0)
                            fref[d01][d02][0][d3][d4] += p.mu * (d3 - 1) * ((p.K - d3 + 1) / p.K) * yref[d01][d02][0][d3 - 1][d4];
                        if (d4 > 0)
                            fref[d01][d02][0][d3][d4] += p.beta * (d3 + 1) * truncated_Poisson_coupling(d3 + 1, yref[d01][d02][0][d3 + 2][d4 - 1]) * (d4 - 1) * yref[d01][d02][0][d3 + 1][d4 - 1];
                    }

                    if (d3 < p.MF3 - 3 && d4 == p.MF4 - 3)
                    { // final states of master equation for predator
                        fref[d01][d02][0][d3][d4] = -p.mu * d3 * ((p.K - d3) / p.K) * yref[d01][d02][0][d3][d4] - p.nu * d4 * yref[d01][d02][0][d3][d4] - p.beta * d3 * d4 * yref[d01][d02][0][d3][d4];
                        fref[d01][d02][0][d3][d4] += p.nu * (d4 + 1) * truncated_Poisson_coupling(d4 + 1, yref[d01][d02][0][d3][d4 + 2]) * yref[d01][d02][0][d3][d4 + 1];
                        if (d3 > 0)
                            fref[d01][d02][0][d3][d4] += p.mu * (d3 - 1) * ((p.K - d3 + 1) / p.K) * yref[d01][d02][0][d3 - 1][d4];
                        if (d4 > 0)
                            fref[d01][d02][0][d3][d4] += p.beta * (d3 + 1) * (d4 - 1) * yref[d01][d02][0][d3 + 1][d4 - 1];
                    }

                    if (d3 == p.MF3 - 3 && d4 == p.MF4 - 3)
                    { // final states of master equations for prey + predator
                        fref[d01][d02][0][d3][d4] = -p.mu * d3 * ((p.K - d3) / p.K) * yref[d01][d02][0][d3][d4] - p.nu * d4 * yref[d01][d02][0][d3][d4] - p.beta * d3 * d4 * yref[d01][d02][0][d3][d4];
                        fref[d01][d02][0][d3][d4] += p.nu * (d4 + 1) * truncated_Poisson_coupling(d4 + 1, yref[d01][d02][0][d3][d4 + 2]) * yref[d01][d02][0][d3][d4 + 1];
                        if (d3 > 0)
                            fref[d01][d02][0][d3][d4] += p.mu * (d3 - 1) * ((p.K - d3 + 1) / p.K) * yref[d01][d02][0][d3 - 1][d4];
                        if (d4 > 0)
                            fref[d01][d02][0][d3][d4] += p.beta * (d3 + 1) * truncated_Poisson_coupling(d3 + 1, yref[d01][d02][0][d3 + 2][d4 - 1]) * (d4 - 1) * yref[d01][d02][0][d3 + 1][d4 - 1];
                    }

                    if (d3 == p.MF3 - 2 && d4 < p.MF4 - 3)
                    { // prey mean-field equation, general predator master equation
                        fref[d01][d02][0][d3][d4] = -p.nu * d4 * yref[d01][d02][0][d3][d4];
                        fref[d01][d02][0][d3][d4] -= p.beta * ((1.0 - truncated_Poisson_coupling(d3, yref[d01][d02][0][d3 + 1][d4])) * yref[d01][d02][0][d3 + 1][d4] + truncated_Poisson_coupling(d3, yref[d01][d02][0][d3 + 1][d4]) * d3) * d4 * yref[d01][d02][0][d3][d4];
                        fref[d01][d02][0][d3][d4] += p.nu * (d4 + 1) * yref[d01][d02][0][d3][d4 + 1];
                        if (d3 > 0)
                            fref[d01][d02][0][d3][d4] += p.mu * (d3 - 1) * ((p.K - d3 + 1) / p.K) * yref[d01][d02][0][d3 - 1][d4];
                        if (d4 > 0)
                            fref[d01][d02][0][d3][d4] += p.beta * yref[d01][d02][0][d3 + 1][d4 - 1] * (1.0 - truncated_Poisson_coupling(d3, yref[d01][d02][0][d3 + 1][d4 - 1])) * (d4 - 1) * yref[d01][d02][0][d3][d4 - 1];

                        fref[d01][d02][0][d3 + 1][d4] = p.mu * yref[d01][d02][0][d3 + 1][d4] * (p.K - yref[d01][d02][0][d3 + 1][d4]) / p.K - p.beta * d4 * yref[d01][d02][0][d3 + 1][d4]; // mean-field state
                    }

                    if (d3 == p.MF3 - 2 && d4 == p.MF4 - 3)
                    { // prey mean-field equation, final predator master equation state
                        fref[d01][d02][0][d3][d4] = -p.nu * d4 * yref[d01][d02][0][d3][d4];
                        fref[d01][d02][0][d3][d4] -= p.beta * ((1.0 - truncated_Poisson_coupling(d3, yref[d01][d02][0][d3 + 1][d4])) * yref[d01][d02][0][d3 + 1][d4] + truncated_Poisson_coupling(d3, yref[d01][d02][0][d3 + 1][d4]) * d3) * d4 * yref[d01][d02][0][d3][d4];
                        fref[d01][d02][0][d3][d4] += p.nu * (d4 + 1) * truncated_Poisson_coupling(d4 + 1, yref[d01][d02][0][d3][d4 + 2]) * yref[d01][d02][0][d3][d4 + 1];
                        if (d3 > 0)
                            fref[d01][d02][0][d3][d4] += p.mu * (d3 - 1) * ((p.K - d3 + 1) / p.K) * yref[d01][d02][0][d3 - 1][d4];
                        if (d4 > 0)
                            fref[d01][d02][0][d3][d4] += p.beta * yref[d01][d02][0][d3 + 1][d4 - 1] * (1.0 - truncated_Poisson_coupling(d3, yref[d01][d02][0][d3 + 1][d4 - 1])) * (d4 - 1) * yref[d01][d02][0][d3][d4 - 1];

                        fref[d01][d02][0][d3 + 1][d4] = p.mu * yref[d01][d02][0][d3 + 1][d4] * (p.K - yref[d01][d02][0][d3 + 1][d4]) / p.K - p.beta * d4 * yref[d01][d02][0][d3 + 1][d4]; // mean-field state
                    }

                    if (d3 < p.MF3 - 3 && d4 == p.MF4 - 2)
                    { // general prey master equation, predator mean-field equation
                        fref[d01][d02][0][d3][d4] = -p.mu * d3 * ((p.K - d3) / p.K) * yref[d01][d02][0][d3][d4] - p.nu * d4 * truncated_Poisson_coupling(d4, yref[d01][d02][0][d3][d4 + 1]) * yref[d01][d02][0][d3][d4] - p.beta * d3 * yref[d01][d02][0][d3][d4 + 1] * yref[d01][d02][0][d3][d4];
                        if (d3 > 0)
                            fref[d01][d02][0][d3][d4] += p.mu * (d3 - 1) * ((p.K - d3 + 1) / p.K) * yref[d01][d02][0][d3 - 1][d4];
                        if (d4 > 0)
                            fref[d01][d02][0][d3][d4] += p.beta * (d3 + 1) * (d4 - 1) * yref[d01][d02][0][d3 + 1][d4 - 1];
                        if (d4 > 0)
                            fref[d01][d02][0][d3][d4] += p.beta * (d3 + 1) * yref[d01][d02][0][d3 + 1][d4 + 1] * yref[d01][d02][0][d3 + 1][d4];

                        fref[d01][d02][0][d3][d4 + 1] = p.beta * d3 * yref[d01][d02][0][d3][d4 + 1] - p.nu * yref[d01][d02][0][d3][d4 + 1]; // mean-field state
                    }

                    if (d3 == p.MF3 - 3 && d4 == p.MF4 - 2)
                    { // final prey master equation state, predator mean-field equation
                        fref[d01][d02][0][d3][d4] = -p.mu * d3 * ((p.K - d3) / p.K) * yref[d01][d02][0][d3][d4] - p.nu * d4 * truncated_Poisson_coupling(d4, yref[d01][d02][0][d3][d4 + 1]) * yref[d01][d02][0][d3][d4] - p.beta * d3 * yref[d01][d02][0][d3][d4 + 1] * yref[d01][d02][0][d3][d4];
                        if (d3 > 0)
                            fref[d01][d02][0][d3][d4] += p.mu * (d3 - 1) * ((p.K - d3 + 1) / p.K) * yref[d01][d02][0][d3 - 1][d4];
                        if (d4 > 0)
                            fref[d01][d02][0][d3][d4] += p.beta * (d3 + 1) * (d4 - 1) * truncated_Poisson_coupling(d3 + 1, yref[d01][d02][0][d3 + 2][d4 - 1]) * yref[d01][d02][0][d3 + 1][d4 - 1];
                        if (d4 > 0)
                            fref[d01][d02][0][d3][d4] += p.beta * (d3 + 1) * truncated_Poisson_coupling(d3 + 1, yref[d01][d02][0][d3 + 2][d4]) * yref[d01][d02][0][d3 + 1][d4 + 1] * yref[d01][d02][0][d3 + 1][d4];

                        fref[d01][d02][0][d3][d4 + 1] = p.beta * d3 * yref[d01][d02][0][d3][d4 + 1] - p.nu * yref[d01][d02][0][d3][d4 + 1]; // mean-field state
                    }

                    if (d3 == p.MF3 - 2 && d4 == p.MF4 - 2)
                    { // double mean-regime
                        fref[d01][d02][0][d3][d4] = -p.nu * d4 * truncated_Poisson_coupling(d4, yref[d01][d02][0][d3][d4 + 1]) * yref[d01][d02][0][d3][d4] - p.beta * d3 * truncated_Poisson_coupling(d3, yref[d01][d02][0][d3 + 1][d4]) * yref[d01][d02][0][d3][d4 + 1] * yref[d01][d02][0][d3][d4];
                        if (d3 > 0)
                            fref[d01][d02][0][d3][d4] += p.mu * (d3 - 1) * ((p.K - d3 + 1) / p.K) * yref[d01][d02][0][d3 - 1][d4];
                        if (d4 > 0)
                            fref[d01][d02][0][d3][d4] += p.beta * (d4 - 1) * yref[d01][d02][0][d3 + 1][d4 - 1] * (1.0 - truncated_Poisson_coupling(d3, yref[d01][d02][0][d3 + 1][d4 - 1])) * yref[d01][d02][0][d3][d4 - 1];

                        fref[d01][d02][0][d3][d4 + 1] = p.beta * yref[d01][d02][0][d3 + 1][d4] * yref[d01][d02][0][d3][d4 + 1] - p.nu * yref[d01][d02][0][d3][d4 + 1];                                               // mean-field state
                        fref[d01][d02][0][d3 + 1][d4] = p.mu * yref[d01][d02][0][d3 + 1][d4] * (p.K - yref[d01][d02][0][d3 + 1][d4]) / p.K - p.beta * yref[d01][d02][0][d3][d4 + 1] * yref[d01][d02][0][d3 + 1][d4]; // mean-field state
                    }
                }
            }
        }
    }

    return GSL_SUCCESS;

} //********** end function dydt definition ********************************************************

#endif // DYN_DIFF_HPP_INCLUDED
