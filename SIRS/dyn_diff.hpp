#ifndef DYN_DIFF_HPP_INCLUDED
#define DYN_DIFF_HPP_INCLUDED

#include <boost/multi_array.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <math.h>

using namespace std;

/**
 * @file   dyn_diff.hpp
 * @brief  ODE system for SIRS dynamics in a coupled master-equation + mean-field system
 *
 * @author  LHD
 * @since   2023-10-16
 */

struct Sparam
{
    const double beta;
    const double alpha;
    const double gamma;
    const int N;
    const int MF3;
    const int MF4;
}; // parameter structure

double truncated_Poisson_coupling(int k, const double average)
{

    if (average <= 0 || k <= 0 || isinf(average) || isnan(average))
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

int Population(double time, double x, double y)
{

    // double local_factor = 100.0*abs(cos(x/5)*sin(sqrt(y/2)));
    double local_factor = 100.0 * abs(cos(x / 5) * sin(sqrt(y / 2)) * cos(x / 5) * sin(sqrt(y / 2)));

    return 10 + ceil(local_factor);
}

//********** function dydt definition **************************************************************
int dydt(double t, const double y[], double f[], void *param)
{
    // ODE system for diffusion

    // Cast parameters
    Sparam &p = *static_cast<Sparam *>(param);

    // Create multi_array reference to y and f
    typedef boost::multi_array_ref<const double, 5> CSTmatref_type;
    typedef boost::multi_array_ref<double, 5> matref_type;
    typedef CSTmatref_type::index indexref;
    CSTmatref_type yref(y, boost::extents[p.N][p.N][1][p.MF3][p.MF4]);
    matref_type fref(f, boost::extents[p.N][p.N][1][p.MF3][p.MF4]);

    // Calculate mean-field diffusion
    std::vector<std::vector<double>> diffusion_matrix(p.N, std::vector<double>(p.N, 0.0));
    for (int d1 = 0; d1 < p.N; ++d1)
        for (int d2 = 0; d2 < p.N; ++d2)
        { // loop over sites
            // Draw location-dependent diffusion parameter
            const double diff_rate = 0.1 * p.beta;
            for (int d4 = 0; d4 < p.MF4 - 1; ++d4)
            { // loop over recovered individuals
                for (int d3 = 0; d3 < p.MF3 - 2; ++d3)
                { // loop over infectious individuals
                    // Loop over diffusion kernel (currently von neumann neighborhood)
                    if (d1 > 0)
                        diffusion_matrix[d1][d2] += diff_rate * d3 * yref[d1 - 1][d2][0][d3][d4];
                    if (d1 < p.N - 1)
                        diffusion_matrix[d1][d2] += diff_rate * d3 * yref[d1 + 1][d2][0][d3][d4];
                    if (d2 > 0)
                        diffusion_matrix[d1][d2] += diff_rate * d3 * yref[d1][d2 - 1][0][d3][d4];
                    if (d2 < p.N - 1)
                        diffusion_matrix[d1][d2] += diff_rate * d3 * yref[d1][d2 + 1][0][d3][d4];
                }
                // Diffusion from mean-field limit of infectious individuals
                // Loop over diffusion kernel (currently von neumann neighborhood)
                if (d1 > 0)
                    diffusion_matrix[d1][d2] += diff_rate * yref[d1 - 1][d2][0][p.MF3 - 1][d4] * yref[d1 - 1][d2][0][p.MF3 - 2][d4];
                if (d1 < p.N - 1)
                    diffusion_matrix[d1][d2] += diff_rate * yref[d1 + 1][d2][0][p.MF3 - 1][d4] * yref[d1 + 1][d2][0][p.MF3 - 2][d4];
                if (d2 > 0)
                    diffusion_matrix[d1][d2] += diff_rate * yref[d1][d2 - 1][0][p.MF3 - 1][d4] * yref[d1][d2 - 1][0][p.MF3 - 2][d4];
                if (d2 < p.N - 1)
                    diffusion_matrix[d1][d2] += diff_rate * yref[d1][d2 + 1][0][p.MF3 - 1][d4] * yref[d1][d2 + 1][0][p.MF3 - 2][d4];
            }
        }

    // Process all dynamics
    for (int d1 = 0; d1 < p.N; ++d1)
    {
        for (int d2 = 0; d2 < p.N; ++d2)
        {
            for (int d3 = 0; d3 < p.MF3 - 1; ++d3)
            {
                for (int d4 = 0; d4 < p.MF4 - 1; ++d4)
                {
                    if (d3 < p.MF3 - 3 && d4 < p.MF4 - 3)
                    { // general master equation
                        double S = Population(t, d1, d2) - d3 - d4;
                        fref[d1][d2][0][d3][d4] = -(p.beta * d3 + diffusion_matrix[d1][d2]) * S * yref[d1][d2][0][d3][d4] - 1.0 * d3 * p.alpha * yref[d1][d2][0][d3][d4];
                        fref[d1][d2][0][d3][d4] += p.gamma * (d4 + 1) * yref[d1][d2][0][d3][d4 + 1] - p.gamma * d4 * yref[d1][d2][0][d3][d4];
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] += (p.beta * (d3 - 1) + diffusion_matrix[d1][d2]) * (S + 1) * yref[d1][d2][0][d3 - 1][d4];
                        if (d4 > 0)
                            fref[d1][d2][0][d3][d4] += 1.0 * (d3 + 1) * p.alpha * yref[d1][d2][0][d3 + 1][d4 - 1];
                    }

                    if (d3 == p.MF3 - 3 && d4 < p.MF4 - 3)
                    { // final states of master equation for infectious
                        double S = Population(t, d1, d2) - d3 - d4;
                        fref[d1][d2][0][d3][d4] = -(p.beta * d3 + diffusion_matrix[d1][d2]) * S * yref[d1][d2][0][d3][d4] - 1.0 * d3 * p.alpha * yref[d1][d2][0][d3][d4];
                        fref[d1][d2][0][d3][d4] += p.gamma * (d4 + 1) * yref[d1][d2][0][d3][d4 + 1] - p.gamma * d4 * yref[d1][d2][0][d3][d4];
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] += (p.beta * (d3 - 1) + diffusion_matrix[d1][d2]) * (S + 1) * yref[d1][d2][0][d3 - 1][d4];
                        if (d4 > 0)
                            fref[d1][d2][0][d3][d4] += truncated_Poisson_coupling(d3 + 1, yref[d1][d2][0][d3 + 2][d4 - 1]) * (d3 + 1) * p.alpha * yref[d1][d2][0][d3 + 1][d4 - 1];
                    }

                    if (d3 < p.MF3 - 3 && d4 == p.MF4 - 3)
                    { // final states of master equation for recovered
                        double S = Population(t, d1, d2) - d3 - d4;
                        fref[d1][d2][0][d3][d4] = -(p.beta * d3 + diffusion_matrix[d1][d2]) * S * yref[d1][d2][0][d3][d4] - 1.0 * d3 * p.alpha * yref[d1][d2][0][d3][d4];
                        fref[d1][d2][0][d3][d4] += p.gamma * (d4 + 1) * yref[d1][d2][0][d3][d4 + 1] * truncated_Poisson_coupling(d4 + 1, yref[d1][d2][0][d3][d4 + 2]) - p.gamma * d4 * yref[d1][d2][0][d3][d4];
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] += (p.beta * (d3 - 1) + diffusion_matrix[d1][d2]) * (S + 1) * yref[d1][d2][0][d3 - 1][d4];
                        if (d4 > 0)
                            fref[d1][d2][0][d3][d4] += 1.0 * (d3 + 1) * p.alpha * yref[d1][d2][0][d3 + 1][d4 - 1];
                    }

                    if (d3 == p.MF3 - 3 && d4 == p.MF4 - 3)
                    { // final states of master equations for infectious + recovered
                        double S = Population(t, d1, d2) - d3 - d4;
                        fref[d1][d2][0][d3][d4] = -(p.beta * d3 + diffusion_matrix[d1][d2]) * S * yref[d1][d2][0][d3][d4] - 1.0 * d3 * p.alpha * yref[d1][d2][0][d3][d4];
                        fref[d1][d2][0][d3][d4] += p.gamma * (d4 + 1) * yref[d1][d2][0][d3][d4 + 1] * truncated_Poisson_coupling(d4 + 1, yref[d1][d2][0][d3][d4 + 2]) - p.gamma * d4 * yref[d1][d2][0][d3][d4];
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] += (p.beta * (d3 - 1) + diffusion_matrix[d1][d2]) * (S + 1) * yref[d1][d2][0][d3 - 1][d4];
                        if (d4 > 0)
                            fref[d1][d2][0][d3][d4] += truncated_Poisson_coupling(d3 + 1, yref[d1][d2][0][d3 + 2][d4 - 1]) * (d3 + 1) * p.alpha * yref[d1][d2][0][d3 + 1][d4 - 1];
                    }

                    if (d3 == p.MF3 - 2 && d4 < p.MF4 - 3)
                    { // infectious mean-field equation, general recovered master equation
                        double S = Population(t, d1, d2) - yref[d1][d2][0][d3 + 1][d4] - d4;
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] = (p.beta * (d3 - 1) + diffusion_matrix[d1][d2]) * (Population(t, d1, d2) - d3 + 1 - d4) * yref[d1][d2][0][d3 - 1][d4];
                        fref[d1][d2][0][d3][d4] += p.gamma * (d4 + 1) * yref[d1][d2][0][d3][d4 + 1] - p.gamma * d4 * yref[d1][d2][0][d3][d4];
                        fref[d1][d2][0][d3][d4] -= p.alpha * yref[d1][d2][0][d3 + 1][d4] * (1.0 - truncated_Poisson_coupling(d3, yref[d1][d2][0][d3 + 1][d4])) * yref[d1][d2][0][d3][d4];
                        fref[d1][d2][0][d3][d4] -= p.alpha * d3 * truncated_Poisson_coupling(d3, yref[d1][d2][0][d3 + 1][d4]) * yref[d1][d2][0][d3][d4];
                        if (d4 > 0)
                            fref[d1][d2][0][d3][d4] += (1.0 - truncated_Poisson_coupling(d3, yref[d1][d2][0][d3 + 1][d4 - 1])) * yref[d1][d2][0][d3 + 1][d4 - 1] * p.alpha * yref[d1][d2][0][d3][d4 - 1];

                        fref[d1][d2][0][d3 + 1][d4] = (p.beta * yref[d1][d2][0][d3 + 1][d4] + diffusion_matrix[d1][d2]) * S - yref[d1][d2][0][d3 + 1][d4] * p.alpha; // mean-field state
                    }

                    if (d3 == p.MF3 - 2 && d4 == p.MF4 - 3)
                    { // infectious mean-field equation, final recovered master equation state
                        double S = Population(t, d1, d2) - yref[d1][d2][0][d3 + 1][d4] - d4;
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] = (p.beta * (d3 - 1) + diffusion_matrix[d1][d2]) * (Population(t, d1, d2) - d3 + 1 - d4) * yref[d1][d2][0][d3 - 1][d4];
                        fref[d1][d2][0][d3][d4] += p.gamma * (d4 + 1) * yref[d1][d2][0][d3][d4 + 1] * truncated_Poisson_coupling(d4 + 1, yref[d1][d2][0][d3][d4 + 2]) - p.gamma * d4 * yref[d1][d2][0][d3][d4];
                        fref[d1][d2][0][d3][d4] -= p.alpha * yref[d1][d2][0][d3 + 1][d4] * (1.0 - truncated_Poisson_coupling(d3, yref[d1][d2][0][d3 + 1][d4])) * yref[d1][d2][0][d3][d4];
                        fref[d1][d2][0][d3][d4] -= p.alpha * d3 * truncated_Poisson_coupling(d3, yref[d1][d2][0][d3 + 1][d4]) * yref[d1][d2][0][d3][d4];
                        if (d4 > 0)
                            fref[d1][d2][0][d3][d4] += (1.0 - truncated_Poisson_coupling(d3, yref[d1][d2][0][d3 + 1][d4 - 1])) * yref[d1][d2][0][d3 + 1][d4 - 1] * p.alpha * yref[d1][d2][0][d3][d4 - 1];

                        fref[d1][d2][0][d3 + 1][d4] = (p.beta * yref[d1][d2][0][d3 + 1][d4] + diffusion_matrix[d1][d2]) * S - yref[d1][d2][0][d3 + 1][d4] * p.alpha; // mean-field state
                    }

                    if (d3 < p.MF3 - 3 && d4 == p.MF4 - 2)
                    { // general infectious master equation, recovered mean-field equation
                        double S = Population(t, d1, d2) - d3 - yref[d1][d2][0][d3][d4 + 1];
                        fref[d1][d2][0][d3][d4] = -(p.beta * d3 + diffusion_matrix[d1][d2]) * S * yref[d1][d2][0][d3][d4] - d3 * p.alpha * yref[d1][d2][0][d3][d4];
                        fref[d1][d2][0][d3][d4] -= p.gamma * d4 * yref[d1][d2][0][d3][d4] * truncated_Poisson_coupling(d4, yref[d1][d2][0][d3][d4 + 1]);
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] += (p.beta * (d3 - 1) + diffusion_matrix[d1][d2]) * (Population(t, d1, d2) - d3 + 1 - yref[d1][d2][0][d3 - 1][d4 + 1]) * yref[d1][d2][0][d3 - 1][d4];
                        if (d4 > 0)
                            fref[d1][d2][0][d3][d4] += 1.0 * (d3 + 1) * p.alpha * yref[d1][d2][0][d3 + 1][d4 - 1];
                        fref[d1][d2][0][d3][d4] += 1.0 * (d3 + 1) * p.alpha * yref[d1][d2][0][d3 + 1][d4];

                        fref[d1][d2][0][d3][d4 + 1] = 1.0 * (d3 + 1) * p.alpha * yref[d1][d2][0][d3 + 1][d4];
                        if (d4 > 0)
                            fref[d1][d2][0][d3][d4 + 1] += 1.0 * (d3 + 1) * p.alpha * yref[d1][d2][0][d3 + 1][d4 - 1];
                        fref[d1][d2][0][d3][d4 + 1] -= p.gamma * yref[d1][d2][0][d3][d4 + 1] * yref[d1][d2][0][d3][d4]; // mean-field state
                    }

                    if (d3 == p.MF3 - 3 && d4 == p.MF4 - 2)
                    { // final infectious master equation state, recovered mean-field equation
                        double S = Population(t, d1, d2) - d3 - yref[d1][d2][0][d3][d4 + 1];
                        fref[d1][d2][0][d3][d4] = -(p.beta * d3 + diffusion_matrix[d1][d2]) * S * yref[d1][d2][0][d3][d4] - d3 * p.alpha * yref[d1][d2][0][d3][d4];
                        fref[d1][d2][0][d3][d4] -= p.gamma * d4 * yref[d1][d2][0][d3][d4] * truncated_Poisson_coupling(d4, yref[d1][d2][0][d3][d4 + 1]);
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] += (p.beta * (d3 - 1) + diffusion_matrix[d1][d2]) * (Population(t, d1, d2) - d3 + 1 - yref[d1][d2][0][d3 - 1][d4 + 1]) * yref[d1][d2][0][d3 - 1][d4];
                        if (d4 > 0)
                            fref[d1][d2][0][d3][d4] += 1.0 * (d3 + 1) * truncated_Poisson_coupling(d3 + 1, yref[d1][d2][0][d3 + 2][d4 - 1]) * p.alpha * yref[d1][d2][0][d3 + 1][d4 - 1];
                        fref[d1][d2][0][d3][d4] += 1.0 * (d3 + 1) * truncated_Poisson_coupling(d3 + 1, yref[d1][d2][0][d3 + 2][d4]) * p.alpha * yref[d1][d2][0][d3 + 1][d4];

                        fref[d1][d2][0][d3][d4 + 1] = 1.0 * (d3 + 1) * truncated_Poisson_coupling(d3 + 1, yref[d1][d2][0][d3 + 2][d4]) * p.alpha * yref[d1][d2][0][d3 + 1][d4];
                        if (d4 > 0)
                            fref[d1][d2][0][d3][d4 + 1] += 1.0 * (d3 + 1) * truncated_Poisson_coupling(d3 + 1, yref[d1][d2][0][d3 + 2][d4 - 1]) * p.alpha * yref[d1][d2][0][d3 + 1][d4 - 1];
                        fref[d1][d2][0][d3][d4 + 1] -= p.gamma * yref[d1][d2][0][d3][d4 + 1] * yref[d1][d2][0][d3][d4]; // mean-field state
                    }

                    if (d3 == p.MF3 - 2 && d4 == p.MF4 - 2)
                    { // double mean-regime
                        double S = Population(t, d1, d2) - yref[d1][d2][0][d3 + 1][d4] - yref[d1][d2][0][d3][d4 + 1];
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] = -1.0 * d3 * p.alpha * yref[d1][d2][0][d3][d4] * truncated_Poisson_coupling(d3, yref[d1][d2][0][d3 + 1][d4]);
                        if (d4 > 0)
                            fref[d1][d2][0][d3][d4] -= p.gamma * d4 * yref[d1][d2][0][d3][d4] * truncated_Poisson_coupling(d4, yref[d1][d2][0][d3][d4 + 1]);
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] += (p.beta * (d3 - 1) + diffusion_matrix[d1][d2]) * (Population(t, d1, d2) - d3 + 1 - yref[d1][d2][0][d3 - 1][d4 + 1]) * yref[d1][d2][0][d3 - 1][d4];
                        if (d4 > 0)
                            fref[d1][d2][0][d3][d4] += (1.0 - truncated_Poisson_coupling(d3, yref[d1][d2][0][d3 + 1][d4 - 1])) * yref[d1][d2][0][d3 + 1][d4 - 1] * p.alpha * yref[d1][d2][0][d3][d4 - 1];

                        fref[d1][d2][0][d3 + 1][d4] = (p.beta * yref[d1][d2][0][d3 + 1][d4] + diffusion_matrix[d1][d2]) * S - p.alpha * yref[d1][d2][0][d3 + 1][d4]; // mean-field state
                        fref[d1][d2][0][d3][d4 + 1] = p.alpha * yref[d1][d2][0][d3 + 1][d4] - p.gamma * yref[d1][d2][0][d3][d4 + 1];                                   // mean-field state
                    }
                }
            }
        }
    }

    return GSL_SUCCESS;

} //********** end function dydt definition ********************************************************

#endif // DYN_DIFF_HPP_INCLUDED
