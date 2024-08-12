#ifndef DYN_DIFF_HPP_INCLUDED
#define DYN_DIFF_HPP_INCLUDED

#include <boost/multi_array.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <math.h>

using namespace std;

/**
 * @file   dyn_diff.hpp
 * @brief  ODE system for species dispersion in a coupled master-equation + mean-field system
 *
 * @author  LHD
 * @since   2023-10-16
 */

struct Sparam
{
    const double beta;
    const double seed2seedling;
    const double seeddeath;
    const double seedling2sapling;
    const double seedlingdeath;
    const double sapling2tree;
    const double sapplingdeath;
    const double treedeath;
    const int N;
    const int MF3;
    const int MF4;
}; // parameter structure

double Poisson_coupling(int k, const double average)
{

    // This needs something to deal with the case where k is large and average is low?

    if (average <= 0 || k <= 0 || isinf(average) || isnan(average))
        return 0.0;
    else
    {

        boost::math::poisson_distribution<> poiss(average);

        double answer = boost::math::pdf(poiss, k) / (1.0 - boost::math::cdf(poiss, k - 1));

        if (isinf(answer) || answer > 1.0)
            answer = 1.0; // to conserve flow if cdf approx fails
        else if (isnan(answer) || answer < 0.0)
            answer = 0.0; // For Poisson too close to zero

        return answer;
    }
}

double CarryingCapacity(double time, double x, double y)
{

    // double local_factor = cos(abs(x/5)+abs(y/5))*(sqrt(x)+abs(y));
    double local_factor = 20.0 * abs(pow(cos((x + time) / 5), 2.0) * pow(sin(sqrt((y + time) / 2)), 2.0));

    return 0.0 + abs(local_factor);
}

double Growth(double pop, double K)
{
    double rate = 1.0 - pop / K;
    if (rate < 0 || K == 0.0)
        return 0.0;
    else
        return rate;
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

    CSTmatref_type yref(y, boost::extents[p.N][p.N][3][p.MF3][p.MF4]);
    matref_type fref(f, boost::extents[p.N][p.N][3][p.MF3][p.MF4]);

    // Calculate mean-field diffusion
    std::vector<std::vector<double>> diffusion_matrix(p.N, std::vector<double>(p.N, 0.0));
    for (int d01 = 0; d01 < p.N; ++d01)
        for (int d02 = 0; d02 < p.N; ++d02)
        { // loop over sites
            // Draw climate-dependent diffusion parameter
            const double diff_rate = p.beta;
            for (int d3 = 0; d3 < p.MF3 - 1; ++d3)
            { // loop over sappling numbers
                // Diffusion from master equation of adult trees
                for (int d4 = 0; d4 < p.MF4 - 2; ++d4)
                {
                    // Loop over diffusion kernel (currently von neumann neighborhood)
                    if (d01 > 0)
                        diffusion_matrix[d01][d02] += diff_rate * d4 * yref[d01 - 1][d02][0][d3][d4];
                    if (d01 < p.N - 1)
                        diffusion_matrix[d01][d02] += diff_rate * d4 * yref[d01 + 1][d02][0][d3][d4];
                    if (d02 > 0)
                        diffusion_matrix[d01][d02] += diff_rate * d4 * yref[d01][d02 - 1][0][d3][d4];
                    if (d02 < p.N - 1)
                        diffusion_matrix[d01][d02] += diff_rate * d4 * yref[d01][d02 + 1][0][d3][d4];
                }
                // Diffusion from mean-field limit of adult trees
                // Loop over diffusion kernel (currently von neumann neighborhood)
                if (d01 > 0)
                    diffusion_matrix[d01][d02] += diff_rate * yref[d01 - 1][d02][0][d3][p.MF4 - 1] * yref[d01 - 1][d02][0][d3][p.MF4 - 2];
                if (d01 < p.N - 1)
                    diffusion_matrix[d01][d02] += diff_rate * yref[d01 + 1][d02][0][d3][p.MF4 - 1] * yref[d01 + 1][d02][0][d3][p.MF4 - 2];
                if (d02 > 0)
                    diffusion_matrix[d01][d02] += diff_rate * yref[d01][d02 - 1][0][d3][p.MF4 - 1] * yref[d01][d02 - 1][0][d3][p.MF4 - 2];
                if (d02 < p.N - 1)
                    diffusion_matrix[d01][d02] += diff_rate * yref[d01][d02 + 1][0][d3][p.MF4 - 1] * yref[d01][d02 + 1][0][d3][p.MF4 - 2];
            }
        }

    // Add external seeding and compute seed + seedling mean-fields
    for (int d01 = 0; d01 < p.N; ++d01)
    {
        for (int d02 = 0; d02 < p.N; ++d02)
        { // loop over sites
            // Draw climate-dependent diffusion parameter
            const double diff_rate = p.beta;
            // double K = Suitability(t,d01,d02)*100;
            for (int d3 = 0; d3 < p.MF3 - 1; ++d3)
            {
                // Dynamics with master number of adult trees
                for (int d4 = 0; d4 < p.MF4 - 2; ++d4)
                {
                    fref[d01][d02][1][d3][d4] = diffusion_matrix[d01][d02] + diff_rate * d4 * yref[d01][d02][0][d3][d4] - yref[d01][d02][1][d3][d4] * (p.seeddeath + p.seed2seedling);
                    fref[d01][d02][2][d3][d4] = yref[d01][d02][1][d3][d4] * p.seed2seedling - yref[d01][d02][2][d3][d4] * (p.seedlingdeath + p.seedling2sapling);
                }
                // Plus dynamics with mean-field adult trees
                fref[d01][d02][1][d3][p.MF4 - 2] = diffusion_matrix[d01][d02] + diff_rate * yref[d01][d02][0][d3][p.MF4 - 1] * yref[d01][d02][0][d3][p.MF4 - 2] - yref[d01][d02][1][d3][p.MF4 - 2] * (p.seeddeath + p.seed2seedling);
                fref[d01][d02][2][d3][p.MF4 - 2] = yref[d01][d02][1][d3][p.MF4 - 2] * p.seed2seedling - yref[d01][d02][2][d3][p.MF4 - 2] * (p.seedlingdeath + p.seedling2sapling);
            }
        }
    }

    // Process all internal dynamics for sapplings and trees
    for (int d01 = 0; d01 < p.N; ++d01)
    {
        for (int d02 = 0; d02 < p.N; ++d02)
        {

            double K = CarryingCapacity(t, d02, d01);

            for (int d3 = 0; d3 < p.MF3 - 1; ++d3)
            {
                for (int d4 = 0; d4 < p.MF4 - 1; ++d4)
                {

                    if (d3 < p.MF3 - 3 && d4 < p.MF4 - 3)
                    { // general master equation
                        fref[d01][d02][0][d3][d4] = 0.0;
                        fref[d01][d02][0][d3][d4] = yref[d01][d02][0][d3 + 1][d4] * (d3 + 1) * p.sapplingdeath; // sapling death influx
                        if (d3 > 0)
                            fref[d01][d02][0][d3][d4] += p.seedling2sapling * yref[d01][d02][2][d3 - 1][d4] * yref[d01][d02][0][d3 - 1][d4];
                        fref[d01][d02][0][d3][d4] -= p.seedling2sapling * yref[d01][d02][2][d3][d4] * yref[d01][d02][0][d3][d4]; // seedling to sapling mean-field forcing
                        if (d3 > 0)
                            fref[d01][d02][0][d3][d4] -= yref[d01][d02][0][d3][d4] * d3 * p.sapplingdeath + yref[d01][d02][0][d3][d4] * p.sapling2tree * pow(1.0 * d3, 0.5) * Growth(1.0 * d4, K); // sapling loss outflux
                        if (d4 > 0)
                            fref[d01][d02][0][d3][d4] += yref[d01][d02][0][d3 + 1][d4 - 1] * p.sapling2tree * pow(1.0 + d3, 0.5) * Growth(1.0 * d4 - 1.0, K); // sapling growth influx
                        fref[d01][d02][0][d3][d4] -= yref[d01][d02][0][d3][d4] * d4 * p.treedeath;                                                            // tree loss outflux
                        fref[d01][d02][0][d3][d4] += yref[d01][d02][0][d3][d4 + 1] * (d4 + 1) * p.treedeath;                                                  // tree loss influx
                    }

                    if (d3 == p.MF3 - 3 && d4 < p.MF4 - 3)
                    { // final states of master equation for sappling
                        fref[d01][d02][0][d3][d4] = 0.0;
                        if (d3 > 0)
                            fref[d01][d02][0][d3][d4] = p.seedling2sapling * (yref[d01][d02][2][d3 - 1][d4] * yref[d01][d02][0][d3 - 1][d4] - yref[d01][d02][2][d3][d4] * yref[d01][d02][0][d3][d4]); // seedling to sapling mean-field forcing
                        if (d3 > 0)
                            fref[d01][d02][0][d3][d4] -= yref[d01][d02][0][d3][d4] * d3 * p.sapplingdeath + yref[d01][d02][0][d3][d4] * p.sapling2tree * pow(1.0 * d3, 0.5) * Growth(1.0 * d4, K); // sapling loss outflux
                        fref[d01][d02][0][d3][d4] += yref[d01][d02][0][d3 + 1][d4] * (d3 + 1) * p.sapplingdeath * Poisson_coupling(d3 + 1, yref[d01][d02][0][d3 + 2][d4]);                         // sapling death influx
                        if (d4 > 0)
                            fref[d01][d02][0][d3][d4] += yref[d01][d02][0][d3 + 1][d4 - 1] * pow(1.0 + d3, 0.5) * p.sapling2tree * Growth(1.0 * d4 - 1.0, K) * Poisson_coupling(d3 + 1, yref[d01][d02][0][d3 + 2][d4 - 1]); // sapling growth influx
                        fref[d01][d02][0][d3][d4] -= yref[d01][d02][0][d3][d4] * d4 * p.treedeath;                                                                                                                          // tree loss outflux
                        fref[d01][d02][0][d3][d4] += yref[d01][d02][0][d3][d4 + 1] * (d4 + 1) * p.treedeath;                                                                                                                // tree loss influx
                    }

                    if (d3 < p.MF3 - 3 && d4 == p.MF4 - 3)
                    { // final states of master equation for trees
                        fref[d01][d02][0][d3][d4] = 0.0;
                        if (d3 > 0)
                            fref[d01][d02][0][d3][d4] = -yref[d01][d02][0][d3][d4] * d3 * p.sapplingdeath - yref[d01][d02][0][d3][d4] * p.sapling2tree * pow(1.0 * d3, 0.5) * Growth(1.0 * d4, K); // sapling loss outflux
                        if (d3 > 0)
                            fref[d01][d02][0][d3][d4] += p.seedling2sapling * yref[d01][d02][2][d3 - 1][d4] * yref[d01][d02][0][d3 - 1][d4];
                        fref[d01][d02][0][d3][d4] -= p.seedling2sapling * yref[d01][d02][2][d3][d4] * yref[d01][d02][0][d3][d4]; // seedling to sapling mean-field forcing
                        fref[d01][d02][0][d3][d4] += yref[d01][d02][0][d3 + 1][d4] * (d3 + 1) * p.sapplingdeath;                 // sapling death influx
                        if (d4 > 0)
                            fref[d01][d02][0][d3][d4] += yref[d01][d02][0][d3 + 1][d4 - 1] * p.sapling2tree * pow(1.0 + d3, 0.5) * Growth(1.0 * d4 - 1.0, K);          // sapling growth influx
                        fref[d01][d02][0][d3][d4] -= yref[d01][d02][0][d3][d4] * d4 * p.treedeath;                                                                     // tree loss outflux
                        fref[d01][d02][0][d3][d4] += yref[d01][d02][0][d3][d4 + 1] * (d4 + 1) * p.treedeath * Poisson_coupling(d4 + 1, yref[d01][d02][0][d3][d4 + 2]); // tree loss influx
                    }

                    if (d3 == p.MF3 - 3 && d4 == p.MF4 - 3)
                    { // final states of master equations for sappling + trees
                        fref[d01][d02][0][d3][d4] = 0.0;
                        if (d3 > 0)
                            fref[d01][d02][0][d3][d4] = -yref[d01][d02][0][d3][d4] * d3 * p.sapplingdeath - yref[d01][d02][0][d3][d4] * p.sapling2tree * pow(1.0 * d3, 0.5) * Growth(1.0 * d4, K); // sapling loss outflux
                        if (d3 > 0)
                            fref[d01][d02][0][d3][d4] += p.seedling2sapling * yref[d01][d02][2][d3 - 1][d4] * yref[d01][d02][0][d3 - 1][d4];
                        fref[d01][d02][0][d3][d4] -= p.seedling2sapling * yref[d01][d02][2][d3][d4] * yref[d01][d02][0][d3][d4];                                           // seedling to sapling mean-field forcing
                        fref[d01][d02][0][d3][d4] += yref[d01][d02][0][d3 + 1][d4] * (d3 + 1) * p.sapplingdeath * Poisson_coupling(d3 + 1, yref[d01][d02][0][d3 + 2][d4]); // sapling death influx
                        if (d4 > 0)
                            fref[d01][d02][0][d3][d4] += yref[d01][d02][0][d3 + 1][d4 - 1] * p.sapling2tree * pow(1.0 + d3, 0.5) * Growth(1.0 * d4 - 1.0, K) * Poisson_coupling(d3 + 1, yref[d01][d02][0][d3 + 2][d4 - 1]); // sapling growth influx
                        fref[d01][d02][0][d3][d4] -= yref[d01][d02][0][d3][d4] * d4 * p.treedeath;                                                                                                                          // tree loss outflux
                        fref[d01][d02][0][d3][d4] += yref[d01][d02][0][d3][d4 + 1] * (d4 + 1) * p.treedeath * Poisson_coupling(d4 + 1, yref[d01][d02][0][d3][d4 + 2]);                                                      // tree loss influx
                    }

                    if (d3 == p.MF3 - 2 && d4 < p.MF4 - 3)
                    { // sapplings mean-field equation, general adult master equation
                        fref[d01][d02][0][d3][d4] = 0.0;
                        if (d3 > 0)
                            fref[d01][d02][0][d3][d4] = yref[d01][d02][2][d3 - 1][d4] * p.seedling2sapling * yref[d01][d02][0][d3 - 1][d4];                  // occupation number influx
                        fref[d01][d02][0][d3][d4] -= yref[d01][d02][0][d3][d4] * d3 * p.sapplingdeath * Poisson_coupling(d3, yref[d01][d02][0][d3 + 1][d4]); // occupation number outflux
                        if (d3 > 0)
                            fref[d01][d02][0][d3][d4] -= yref[d01][d02][0][d3][d4] * p.sapling2tree * pow(1.0 * d3, 0.5) * Growth(1.0 * d4, K) * Poisson_coupling(d3, yref[d01][d02][0][d3 + 1][d4]); // occupation number outflux (diagonal move)
                        if (yref[d01][d02][0][d3 + 1][d4] > 0)
                            fref[d01][d02][0][d3][d4] -= yref[d01][d02][0][d3][d4] * yref[d01][d02][0][d3 + 1][d4] * p.sapling2tree * pow(yref[d01][d02][0][d3 + 1][d4], -0.5) * Growth(1.0 * d4, K) * (1 - Poisson_coupling(d3, yref[d01][d02][0][d3 + 1][d4])); // occupation number outflux (lateral move)
                        fref[d01][d02][0][d3][d4] -= yref[d01][d02][0][d3][d4] * d4 * p.treedeath;                                                                                                                                                                // tree loss outflux
                        if (d4 > 0 && yref[d01][d02][0][d3 + 1][d4 - 1] > 0)
                            fref[d01][d02][0][d3][d4] += yref[d01][d02][0][d3][d4 - 1] * pow(yref[d01][d02][0][d3 + 1][d4 - 1], 0.5) * p.sapling2tree * Growth(1.0 * d4 - 1.0, K) * (1 - Poisson_coupling(d3, yref[d01][d02][0][d3 + 1][d4 - 1])); // sapling growth influx
                        fref[d01][d02][0][d3][d4] += yref[d01][d02][0][d3][d4 + 1] * (d4 + 1) * p.treedeath;                                                                                                                                       // tree loss influx

                        if (yref[d01][d02][0][d3 + 1][d4] > 0)
                            fref[d01][d02][0][d3 + 1][d4] = yref[d01][d02][2][d3][d4] * p.seedling2sapling - yref[d01][d02][0][d3 + 1][d4] * p.sapplingdeath + p.sapling2tree * pow(yref[d01][d02][0][d3 + 1][d4], 0.5) * Growth(1.0 * d4, K); // mean-field state
                        else
                            fref[d01][d02][0][d3 + 1][d4] = yref[d01][d02][2][d3][d4] * p.seedling2sapling - yref[d01][d02][0][d3 + 1][d4] * (p.sapplingdeath); // mean-field state
                    }

                    if (d3 == p.MF3 - 2 && d4 == p.MF4 - 3)
                    { // sapplings mean-field equation, final adult master equation state
                        fref[d01][d02][0][d3][d4] = 0.0;
                        if (d3 > 0)
                            fref[d01][d02][0][d3][d4] = yref[d01][d02][2][d3 - 1][d4] * p.seedling2sapling * yref[d01][d02][0][d3 - 1][d4];                  // occupation number influx
                        fref[d01][d02][0][d3][d4] -= yref[d01][d02][0][d3][d4] * d3 * p.sapplingdeath * Poisson_coupling(d3, yref[d01][d02][0][d3 + 1][d4]); // occupation number outflux (diagonal move)
                        if (d3 > 0)
                            fref[d01][d02][0][d3][d4] -= yref[d01][d02][0][d3][d4] * p.sapling2tree * pow(1.0 * d3, 0.5) * Growth(1.0 * d4, K) * Poisson_coupling(d3, yref[d01][d02][0][d3 + 1][d4]); // occupation number outflux (diagonal move)
                        if (yref[d01][d02][0][d3 + 1][d4] > 0)
                            fref[d01][d02][0][d3][d4] -= yref[d01][d02][0][d3][d4] * yref[d01][d02][0][d3 + 1][d4] * p.sapling2tree * pow(yref[d01][d02][0][d3 + 1][d4], -0.5) * Growth(1.0 * d4, K) * (1 - Poisson_coupling(d3, yref[d01][d02][0][d3 + 1][d4])); // occupation number outflux (lateral move)
                        fref[d01][d02][0][d3][d4] -= yref[d01][d02][0][d3][d4] * d4 * p.treedeath;                                                                                                                                                                // tree loss outflux
                        if (d4 > 0 && yref[d01][d02][0][d3 + 1][d4 - 1] > 0)
                            fref[d01][d02][0][d3][d4] += yref[d01][d02][0][d3][d4 - 1] * p.sapling2tree * pow(yref[d01][d02][0][d3 + 1][d4 - 1], 0.5) * Growth(1.0 * d4 - 1.0, K) * (1 - Poisson_coupling(d3, yref[d01][d02][0][d3 + 1][d4 - 1])); // sapling growth influx
                        fref[d01][d02][0][d3][d4] += yref[d01][d02][0][d3][d4 + 1] * (d4 + 1) * p.treedeath * Poisson_coupling(d4 + 1, yref[d01][d02][0][d3][d4 + 2]);                                                                             // tree loss influx

                        if (yref[d01][d02][0][d3 + 1][d4] > 0)
                            fref[d01][d02][0][d3 + 1][d4] = yref[d01][d02][2][d3][d4] * p.seedling2sapling - yref[d01][d02][0][d3 + 1][d4] * p.sapplingdeath + p.sapling2tree * pow(yref[d01][d02][0][d3 + 1][d4], 0.5) * Growth(1.0 * d4, K); // mean-field state
                        else
                            fref[d01][d02][0][d3 + 1][d4] = yref[d01][d02][2][d3][d4] * p.seedling2sapling - yref[d01][d02][0][d3 + 1][d4] * (p.sapplingdeath); // mean-field state
                    }

                    if (d3 < p.MF3 - 3 && d4 == p.MF4 - 2)
                    { // general sapplings master equation, adult tree mean-field equation
                        fref[d01][d02][0][d3][d4] = 0.0;
                        if (d4 > 0)
                            fref[d01][d02][0][d3][d4] = yref[d01][d02][0][d3 + 1][d4 - 1] * pow(1.0 + d3, 0.5) * p.sapling2tree * Growth(1.0 * d4 - 1.0, K); // occupation number influx
                        fref[d01][d02][0][d3][d4] -= yref[d01][d02][0][d3][d4] * d4 * p.treedeath * Poisson_coupling(d4, yref[d01][d02][0][d3][d4 + 1]);     // occupation number outflux
                        if (d3 > 0)
                            fref[d01][d02][0][d3][d4] += p.seedling2sapling * yref[d01][d02][2][d3 - 1][d4] * yref[d01][d02][0][d3 - 1][d4];
                        fref[d01][d02][0][d3][d4] -= p.seedling2sapling * yref[d01][d02][2][d3][d4] * yref[d01][d02][0][d3][d4]; // seedling to sapling mean-field forcing
                        if (d3 > 0)
                            fref[d01][d02][0][d3][d4] -= yref[d01][d02][0][d3][d4] * d3 * p.sapplingdeath + yref[d01][d02][0][d3][d4] * p.sapling2tree * pow(1.0 * d3, 0.5) * Growth(yref[d01][d02][0][d3][d4 + 1], K); // sapling loss outflux
                        fref[d01][d02][0][d3][d4] += yref[d01][d02][0][d3 + 1][d4] * p.sapling2tree * pow(1.0 + d3, 0.5) * Growth(yref[d01][d02][0][d3 + 1][d4 + 1], K);                                                // sapling loss influx
                        fref[d01][d02][0][d3][d4] += yref[d01][d02][0][d3 + 1][d4] * (d3 + 1) * p.sapplingdeath;                                                                                                        // sapling death influx

                        if (d3 > 0)
                            fref[d01][d02][0][d3][d4 + 1] = 1.0 * pow(1.0 * d3, 0.5) * p.sapling2tree * Growth(yref[d01][d02][0][d3][d4 + 1], K) - yref[d01][d02][0][d3][d4 + 1] * p.treedeath; // mean-field state
                        else
                            fref[d01][d02][0][d3][d4 + 1] = -yref[d01][d02][0][d3][d4 + 1] * p.treedeath; // mean-field state
                    }

                    if (d3 == p.MF3 - 3 && d4 == p.MF4 - 2)
                    { // final sapplings master equation state, adult tree mean-field equation
                        fref[d01][d02][0][d3][d4] = 0.0;
                        if (d4 > 0)
                            fref[d01][d02][0][d3][d4] = yref[d01][d02][0][d3 + 1][d4 - 1] * p.sapling2tree * pow(1.0 + d3, 0.5) * Growth(1.0 * d4 - 1.0, K) * Poisson_coupling(d3 + 1, yref[d01][d02][0][d3 + 2][d4 - 1]); // occupation number influx
                        fref[d01][d02][0][d3][d4] -= yref[d01][d02][0][d3][d4] * d4 * p.treedeath * Poisson_coupling(d4, yref[d01][d02][0][d3][d4 + 1]);                                                                   // occupation number outflux
                        if (d3 > 0)
                            fref[d01][d02][0][d3][d4] += p.seedling2sapling * yref[d01][d02][2][d3 - 1][d4] * yref[d01][d02][0][d3 - 1][d4];
                        fref[d01][d02][0][d3][d4] -= p.seedling2sapling * yref[d01][d02][2][d3][d4] * yref[d01][d02][0][d3][d4]; // seedling to sapling mean-field forcing
                        if (d3 > 0)
                            fref[d01][d02][0][d3][d4] -= yref[d01][d02][0][d3][d4] * d3 * p.sapplingdeath + yref[d01][d02][0][d3][d4] * p.sapling2tree * pow(1.0 * d3, 0.5) * Growth(yref[d01][d02][0][d3][d4 + 1], K); // sapling loss outflux
                        fref[d01][d02][0][d3][d4] += yref[d01][d02][0][d3 + 1][d4] * p.sapling2tree * pow(1.0 + d3, 0.5) * Growth(yref[d01][d02][0][d3 + 1][d4 + 1], K) * Poisson_coupling(d3 + 1, yref[d01][d02][0][d3 + 2][d4]);
                        ;                                                                                                                                                  // sapling loss influx
                        fref[d01][d02][0][d3][d4] += yref[d01][d02][0][d3 + 1][d4] * (d3 + 1) * p.sapplingdeath * Poisson_coupling(d3 + 1, yref[d01][d02][0][d3 + 2][d4]); // sapling death influx

                        if (d3 > 0)
                            fref[d01][d02][0][d3][d4 + 1] = 1.0 * pow(1.0 * d3, 0.5) * p.sapling2tree * Growth(yref[d01][d02][0][d3][d4 + 1], K) - yref[d01][d02][0][d3][d4 + 1] * p.treedeath; // mean-field state
                        else
                            fref[d01][d02][0][d3][d4 + 1] = -yref[d01][d02][0][d3][d4 + 1] * p.treedeath; // mean-field state
                    }

                    if (d3 == p.MF3 - 2 && d4 == p.MF4 - 2)
                    { // double mean-regime
                        fref[d01][d02][0][d3][d4] = 0.0;
                        if (d3 > 0)
                            fref[d01][d02][0][d3][d4] = yref[d01][d02][2][d3 - 1][d4] * p.seedling2sapling * yref[d01][d02][0][d3 - 1][d4]; // sappling MF occupation number influx
                        if (d3 > 0)
                            fref[d01][d02][0][d3][d4] -= (yref[d01][d02][0][d3][d4] * d3 * p.sapplingdeath + yref[d01][d02][0][d3][d4] * p.sapling2tree * pow(1.0 * d3, 0.5) * Growth(yref[d01][d02][0][d3][d4 + 1], K)) * Poisson_coupling(d3, yref[d01][d02][0][d3 + 1][d4]); // sappling MF occupation number outflux
                        if (d4 > 0 && yref[d01][d02][0][d3 + 1][d4 - 1] > 0)
                            fref[d01][d02][0][d3][d4] += yref[d01][d02][0][d3][d4 - 1] * p.sapling2tree * pow(yref[d01][d02][0][d3 + 1][d4 - 1], 0.5) * Growth(1.0 * d4 - 1.0, K) * (1 - Poisson_coupling(d3, yref[d01][d02][0][d3 + 1][d4 - 1])); // tree MF occupation number influx
                        fref[d01][d02][0][d3][d4] -= yref[d01][d02][0][d3][d4] * d4 * p.treedeath * Poisson_coupling(d4, yref[d01][d02][0][d3][d4 + 1]);                                                                                           // tree MF occupation number outflux

                        if (yref[d01][d02][0][d3 + 1][d4] > 0)
                            fref[d01][d02][0][d3 + 1][d4] = yref[d01][d02][2][d3][d4] * p.seedling2sapling - yref[d01][d02][0][d3 + 1][d4] * p.sapplingdeath - p.sapling2tree * pow(yref[d01][d02][0][d3 + 1][d4], 0.5) * Growth(yref[d01][d02][0][d3][d4 + 1], K); // sappling mean-field state
                        else
                            fref[d01][d02][0][d3 + 1][d4] = yref[d01][d02][2][d3][d4] * p.seedling2sapling - yref[d01][d02][0][d3 + 1][d4] * (p.sapplingdeath); // sappling mean-field state
                        if (yref[d01][d02][0][d3 + 1][d4] > 0)
                            fref[d01][d02][0][d3][d4 + 1] = pow(yref[d01][d02][0][d3 + 1][d4], 0.5) * p.sapling2tree * Growth(yref[d01][d02][0][d3][d4 + 1], K) - yref[d01][d02][0][d3][d4 + 1] * p.treedeath; // tree mean-field state
                        else
                            fref[d01][d02][0][d3][d4 + 1] = -yref[d01][d02][0][d3][d4 + 1] * p.treedeath; // tree mean-field state
                    }
                }
            }
        }
    }

    return GSL_SUCCESS;

} //********** end function dydt definition ********************************************************

#endif // DYN_DIFF_HPP_INCLUDED
