#ifndef DYN_DIFF_HPP_INCLUDED
#define DYN_DIFF_HPP_INCLUDED

#include <boost/multi_array.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <math.h>

using namespace std;

/**
 * @brief  ODE system for species dispersion in a coupled master-equation + mean-field system
 * @author  LHD
 * @since   2023-10-16
 */

struct Sparam
{
    const double tree2seed;
    const double seed2seedling;
    const double seed_death;
    const double seedling2sapling;
    const double seedling_death;
    const double sapling2tree;
    const double sapling_death;
    const double tree_death;
    const int n_l;
    const int n_me1;
    const int n_me2;
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

double carrying_capacity(double time, double x, double y)
{
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

    // Cast parameters
    Sparam &p = *static_cast<Sparam *>(param);

    // Create multi_array reference to y and f
    typedef boost::multi_array_ref<const double, 5> CSTmatref_type;
    typedef boost::multi_array_ref<double, 5> matref_type;
    typedef CSTmatref_type::index indexref;

    CSTmatref_type yref(y, boost::extents[p.n_l][p.n_l][3][p.n_me1][p.n_me2]);
    matref_type fref(f, boost::extents[p.n_l][p.n_l][3][p.n_me1][p.n_me2]);

    // Calculate mean-field diffusion
    std::vector<std::vector<double>> diffusion_matrix(p.n_l, std::vector<double>(p.n_l, 0.0));
    for (int d1 = 0; d1 < p.n_l; ++d1)
        for (int d2 = 0; d2 < p.n_l; ++d2)
        { // loop over sites
            // Draw climate-dependent diffusion parameter
            const double diff_rate = p.tree2seed;
            for (int d3 = 0; d3 < p.n_me1 - 1; ++d3)
            { // loop over sapling numbers
                // Diffusion from master equation of adult trees
                for (int d4 = 0; d4 < p.n_me2 - 2; ++d4)
                {
                    // Loop over diffusion kernel (currently von neumann neighborhood)
                    if (d1 > 0)
                        diffusion_matrix[d1][d2] += diff_rate * d4 * yref[d1 - 1][d2][0][d3][d4];
                    if (d1 < p.n_l - 1)
                        diffusion_matrix[d1][d2] += diff_rate * d4 * yref[d1 + 1][d2][0][d3][d4];
                    if (d2 > 0)
                        diffusion_matrix[d1][d2] += diff_rate * d4 * yref[d1][d2 - 1][0][d3][d4];
                    if (d2 < p.n_l - 1)
                        diffusion_matrix[d1][d2] += diff_rate * d4 * yref[d1][d2 + 1][0][d3][d4];
                }
                // Diffusion from mean-field limit of adult trees
                // Loop over diffusion kernel (currently von neumann neighborhood)
                if (d1 > 0)
                    diffusion_matrix[d1][d2] += diff_rate * yref[d1 - 1][d2][0][d3][p.n_me2 - 1] * yref[d1 - 1][d2][0][d3][p.n_me2 - 2];
                if (d1 < p.n_l - 1)
                    diffusion_matrix[d1][d2] += diff_rate * yref[d1 + 1][d2][0][d3][p.n_me2 - 1] * yref[d1 + 1][d2][0][d3][p.n_me2 - 2];
                if (d2 > 0)
                    diffusion_matrix[d1][d2] += diff_rate * yref[d1][d2 - 1][0][d3][p.n_me2 - 1] * yref[d1][d2 - 1][0][d3][p.n_me2 - 2];
                if (d2 < p.n_l - 1)
                    diffusion_matrix[d1][d2] += diff_rate * yref[d1][d2 + 1][0][d3][p.n_me2 - 1] * yref[d1][d2 + 1][0][d3][p.n_me2 - 2];
            }
        }

    // Add external seeding and compute seed + seedling mean-fields
    for (int d1 = 0; d1 < p.n_l; ++d1)
    {
        for (int d2 = 0; d2 < p.n_l; ++d2)
        { // loop over sites
            // Draw climate-dependent diffusion parameter
            const double diff_rate = p.tree2seed;
            // double K = Suitability(t,d1,d2)*100;
            for (int d3 = 0; d3 < p.n_me1 - 1; ++d3)
            {
                // Dynamics with master number of adult trees
                for (int d4 = 0; d4 < p.n_me2 - 2; ++d4)
                {
                    fref[d1][d2][1][d3][d4] = diffusion_matrix[d1][d2] + diff_rate * d4 * yref[d1][d2][0][d3][d4] - yref[d1][d2][1][d3][d4] * (p.seed_death + p.seed2seedling);
                    fref[d1][d2][2][d3][d4] = yref[d1][d2][1][d3][d4] * p.seed2seedling - yref[d1][d2][2][d3][d4] * (p.seedling_death + p.seedling2sapling);
                }
                // Plus dynamics with mean-field adult trees
                fref[d1][d2][1][d3][p.n_me2 - 2] = diffusion_matrix[d1][d2] + diff_rate * yref[d1][d2][0][d3][p.n_me2 - 1] * yref[d1][d2][0][d3][p.n_me2 - 2] - yref[d1][d2][1][d3][p.n_me2 - 2] * (p.seed_death + p.seed2seedling);
                fref[d1][d2][2][d3][p.n_me2 - 2] = yref[d1][d2][1][d3][p.n_me2 - 2] * p.seed2seedling - yref[d1][d2][2][d3][p.n_me2 - 2] * (p.seedling_death + p.seedling2sapling);
            }
        }
    }

    // Process all internal dynamics for saplings and trees
    for (int d1 = 0; d1 < p.n_l; ++d1)
    {
        for (int d2 = 0; d2 < p.n_l; ++d2)
        {

            double K = carrying_capacity(t, d2, d1);

            for (int d3 = 0; d3 < p.n_me1 - 1; ++d3)
            {
                for (int d4 = 0; d4 < p.n_me2 - 1; ++d4)
                {

                    if (d3 < p.n_me1 - 3 && d4 < p.n_me2 - 3)
                    { // general master equation
                        fref[d1][d2][0][d3][d4] = 0.0;
                        fref[d1][d2][0][d3][d4] = yref[d1][d2][0][d3 + 1][d4] * (d3 + 1) * p.sapling_death; // sapling death influx
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] += p.seedling2sapling * yref[d1][d2][2][d3 - 1][d4] * yref[d1][d2][0][d3 - 1][d4];
                        fref[d1][d2][0][d3][d4] -= p.seedling2sapling * yref[d1][d2][2][d3][d4] * yref[d1][d2][0][d3][d4]; // seedling to sapling mean-field forcing
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] -= yref[d1][d2][0][d3][d4] * d3 * p.sapling_death + yref[d1][d2][0][d3][d4] * p.sapling2tree * pow(1.0 * d3, 0.5) * Growth(1.0 * d4, K); // sapling loss outflux
                        if (d4 > 0)
                            fref[d1][d2][0][d3][d4] += yref[d1][d2][0][d3 + 1][d4 - 1] * p.sapling2tree * pow(1.0 + d3, 0.5) * Growth(1.0 * d4 - 1.0, K); // sapling growth influx
                        fref[d1][d2][0][d3][d4] -= yref[d1][d2][0][d3][d4] * d4 * p.tree_death;                                                            // tree loss outflux
                        fref[d1][d2][0][d3][d4] += yref[d1][d2][0][d3][d4 + 1] * (d4 + 1) * p.tree_death;                                                  // tree loss influx
                    }

                    if (d3 == p.n_me1 - 3 && d4 < p.n_me2 - 3)
                    { // final states of master equation for sapling
                        fref[d1][d2][0][d3][d4] = 0.0;
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] = p.seedling2sapling * (yref[d1][d2][2][d3 - 1][d4] * yref[d1][d2][0][d3 - 1][d4] - yref[d1][d2][2][d3][d4] * yref[d1][d2][0][d3][d4]); // seedling to sapling mean-field forcing
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] -= yref[d1][d2][0][d3][d4] * d3 * p.sapling_death + yref[d1][d2][0][d3][d4] * p.sapling2tree * pow(1.0 * d3, 0.5) * Growth(1.0 * d4, K); // sapling loss outflux
                        fref[d1][d2][0][d3][d4] += yref[d1][d2][0][d3 + 1][d4] * (d3 + 1) * p.sapling_death * Poisson_coupling(d3 + 1, yref[d1][d2][0][d3 + 2][d4]);                         // sapling death influx
                        if (d4 > 0)
                            fref[d1][d2][0][d3][d4] += yref[d1][d2][0][d3 + 1][d4 - 1] * pow(1.0 + d3, 0.5) * p.sapling2tree * Growth(1.0 * d4 - 1.0, K) * Poisson_coupling(d3 + 1, yref[d1][d2][0][d3 + 2][d4 - 1]); // sapling growth influx
                        fref[d1][d2][0][d3][d4] -= yref[d1][d2][0][d3][d4] * d4 * p.tree_death;                                                                                                                          // tree loss outflux
                        fref[d1][d2][0][d3][d4] += yref[d1][d2][0][d3][d4 + 1] * (d4 + 1) * p.tree_death;                                                                                                                // tree loss influx
                    }

                    if (d3 < p.n_me1 - 3 && d4 == p.n_me2 - 3)
                    { // final states of master equation for trees
                        fref[d1][d2][0][d3][d4] = 0.0;
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] = -yref[d1][d2][0][d3][d4] * d3 * p.sapling_death - yref[d1][d2][0][d3][d4] * p.sapling2tree * pow(1.0 * d3, 0.5) * Growth(1.0 * d4, K); // sapling loss outflux
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] += p.seedling2sapling * yref[d1][d2][2][d3 - 1][d4] * yref[d1][d2][0][d3 - 1][d4];
                        fref[d1][d2][0][d3][d4] -= p.seedling2sapling * yref[d1][d2][2][d3][d4] * yref[d1][d2][0][d3][d4]; // seedling to sapling mean-field forcing
                        fref[d1][d2][0][d3][d4] += yref[d1][d2][0][d3 + 1][d4] * (d3 + 1) * p.sapling_death;                 // sapling death influx
                        if (d4 > 0)
                            fref[d1][d2][0][d3][d4] += yref[d1][d2][0][d3 + 1][d4 - 1] * p.sapling2tree * pow(1.0 + d3, 0.5) * Growth(1.0 * d4 - 1.0, K);          // sapling growth influx
                        fref[d1][d2][0][d3][d4] -= yref[d1][d2][0][d3][d4] * d4 * p.tree_death;                                                                     // tree loss outflux
                        fref[d1][d2][0][d3][d4] += yref[d1][d2][0][d3][d4 + 1] * (d4 + 1) * p.tree_death * Poisson_coupling(d4 + 1, yref[d1][d2][0][d3][d4 + 2]); // tree loss influx
                    }

                    if (d3 == p.n_me1 - 3 && d4 == p.n_me2 - 3)
                    { // final states of master equations for sapling + trees
                        fref[d1][d2][0][d3][d4] = 0.0;
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] = -yref[d1][d2][0][d3][d4] * d3 * p.sapling_death - yref[d1][d2][0][d3][d4] * p.sapling2tree * pow(1.0 * d3, 0.5) * Growth(1.0 * d4, K); // sapling loss outflux
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] += p.seedling2sapling * yref[d1][d2][2][d3 - 1][d4] * yref[d1][d2][0][d3 - 1][d4];
                        fref[d1][d2][0][d3][d4] -= p.seedling2sapling * yref[d1][d2][2][d3][d4] * yref[d1][d2][0][d3][d4];                                           // seedling to sapling mean-field forcing
                        fref[d1][d2][0][d3][d4] += yref[d1][d2][0][d3 + 1][d4] * (d3 + 1) * p.sapling_death * Poisson_coupling(d3 + 1, yref[d1][d2][0][d3 + 2][d4]); // sapling death influx
                        if (d4 > 0)
                            fref[d1][d2][0][d3][d4] += yref[d1][d2][0][d3 + 1][d4 - 1] * p.sapling2tree * pow(1.0 + d3, 0.5) * Growth(1.0 * d4 - 1.0, K) * Poisson_coupling(d3 + 1, yref[d1][d2][0][d3 + 2][d4 - 1]); // sapling growth influx
                        fref[d1][d2][0][d3][d4] -= yref[d1][d2][0][d3][d4] * d4 * p.tree_death;                                                                                                                          // tree loss outflux
                        fref[d1][d2][0][d3][d4] += yref[d1][d2][0][d3][d4 + 1] * (d4 + 1) * p.tree_death * Poisson_coupling(d4 + 1, yref[d1][d2][0][d3][d4 + 2]);                                                      // tree loss influx
                    }

                    if (d3 == p.n_me1 - 2 && d4 < p.n_me2 - 3)
                    { // saplings mean-field equation, general adult master equation
                        fref[d1][d2][0][d3][d4] = 0.0;
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] = yref[d1][d2][2][d3 - 1][d4] * p.seedling2sapling * yref[d1][d2][0][d3 - 1][d4];                  // occupation number influx
                        fref[d1][d2][0][d3][d4] -= yref[d1][d2][0][d3][d4] * d3 * p.sapling_death * Poisson_coupling(d3, yref[d1][d2][0][d3 + 1][d4]); // occupation number outflux
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] -= yref[d1][d2][0][d3][d4] * p.sapling2tree * pow(1.0 * d3, 0.5) * Growth(1.0 * d4, K) * Poisson_coupling(d3, yref[d1][d2][0][d3 + 1][d4]); // occupation number outflux (diagonal move)
                        if (yref[d1][d2][0][d3 + 1][d4] > 0)
                            fref[d1][d2][0][d3][d4] -= yref[d1][d2][0][d3][d4] * yref[d1][d2][0][d3 + 1][d4] * p.sapling2tree * pow(yref[d1][d2][0][d3 + 1][d4], -0.5) * Growth(1.0 * d4, K) * (1 - Poisson_coupling(d3, yref[d1][d2][0][d3 + 1][d4])); // occupation number outflux (lateral move)
                        fref[d1][d2][0][d3][d4] -= yref[d1][d2][0][d3][d4] * d4 * p.tree_death;                                                                                                                                                                // tree loss outflux
                        if (d4 > 0 && yref[d1][d2][0][d3 + 1][d4 - 1] > 0)
                            fref[d1][d2][0][d3][d4] += yref[d1][d2][0][d3][d4 - 1] * pow(yref[d1][d2][0][d3 + 1][d4 - 1], 0.5) * p.sapling2tree * Growth(1.0 * d4 - 1.0, K) * (1 - Poisson_coupling(d3, yref[d1][d2][0][d3 + 1][d4 - 1])); // sapling growth influx
                        fref[d1][d2][0][d3][d4] += yref[d1][d2][0][d3][d4 + 1] * (d4 + 1) * p.tree_death;                                                                                                                                       // tree loss influx

                        if (yref[d1][d2][0][d3 + 1][d4] > 0)
                            fref[d1][d2][0][d3 + 1][d4] = yref[d1][d2][2][d3][d4] * p.seedling2sapling - yref[d1][d2][0][d3 + 1][d4] * p.sapling_death + p.sapling2tree * pow(yref[d1][d2][0][d3 + 1][d4], 0.5) * Growth(1.0 * d4, K); // mean-field state
                        else
                            fref[d1][d2][0][d3 + 1][d4] = yref[d1][d2][2][d3][d4] * p.seedling2sapling - yref[d1][d2][0][d3 + 1][d4] * (p.sapling_death); // mean-field state
                    }

                    if (d3 == p.n_me1 - 2 && d4 == p.n_me2 - 3)
                    { // saplings mean-field equation, final adult master equation state
                        fref[d1][d2][0][d3][d4] = 0.0;
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] = yref[d1][d2][2][d3 - 1][d4] * p.seedling2sapling * yref[d1][d2][0][d3 - 1][d4];                  // occupation number influx
                        fref[d1][d2][0][d3][d4] -= yref[d1][d2][0][d3][d4] * d3 * p.sapling_death * Poisson_coupling(d3, yref[d1][d2][0][d3 + 1][d4]); // occupation number outflux (diagonal move)
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] -= yref[d1][d2][0][d3][d4] * p.sapling2tree * pow(1.0 * d3, 0.5) * Growth(1.0 * d4, K) * Poisson_coupling(d3, yref[d1][d2][0][d3 + 1][d4]); // occupation number outflux (diagonal move)
                        if (yref[d1][d2][0][d3 + 1][d4] > 0)
                            fref[d1][d2][0][d3][d4] -= yref[d1][d2][0][d3][d4] * yref[d1][d2][0][d3 + 1][d4] * p.sapling2tree * pow(yref[d1][d2][0][d3 + 1][d4], -0.5) * Growth(1.0 * d4, K) * (1 - Poisson_coupling(d3, yref[d1][d2][0][d3 + 1][d4])); // occupation number outflux (lateral move)
                        fref[d1][d2][0][d3][d4] -= yref[d1][d2][0][d3][d4] * d4 * p.tree_death;                                                                                                                                                                // tree loss outflux
                        if (d4 > 0 && yref[d1][d2][0][d3 + 1][d4 - 1] > 0)
                            fref[d1][d2][0][d3][d4] += yref[d1][d2][0][d3][d4 - 1] * p.sapling2tree * pow(yref[d1][d2][0][d3 + 1][d4 - 1], 0.5) * Growth(1.0 * d4 - 1.0, K) * (1 - Poisson_coupling(d3, yref[d1][d2][0][d3 + 1][d4 - 1])); // sapling growth influx
                        fref[d1][d2][0][d3][d4] += yref[d1][d2][0][d3][d4 + 1] * (d4 + 1) * p.tree_death * Poisson_coupling(d4 + 1, yref[d1][d2][0][d3][d4 + 2]);                                                                             // tree loss influx

                        if (yref[d1][d2][0][d3 + 1][d4] > 0)
                            fref[d1][d2][0][d3 + 1][d4] = yref[d1][d2][2][d3][d4] * p.seedling2sapling - yref[d1][d2][0][d3 + 1][d4] * p.sapling_death + p.sapling2tree * pow(yref[d1][d2][0][d3 + 1][d4], 0.5) * Growth(1.0 * d4, K); // mean-field state
                        else
                            fref[d1][d2][0][d3 + 1][d4] = yref[d1][d2][2][d3][d4] * p.seedling2sapling - yref[d1][d2][0][d3 + 1][d4] * (p.sapling_death); // mean-field state
                    }

                    if (d3 < p.n_me1 - 3 && d4 == p.n_me2 - 2)
                    { // general saplings master equation, adult tree mean-field equation
                        fref[d1][d2][0][d3][d4] = 0.0;
                        if (d4 > 0)
                            fref[d1][d2][0][d3][d4] = yref[d1][d2][0][d3 + 1][d4 - 1] * pow(1.0 + d3, 0.5) * p.sapling2tree * Growth(1.0 * d4 - 1.0, K); // occupation number influx
                        fref[d1][d2][0][d3][d4] -= yref[d1][d2][0][d3][d4] * d4 * p.tree_death * Poisson_coupling(d4, yref[d1][d2][0][d3][d4 + 1]);     // occupation number outflux
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] += p.seedling2sapling * yref[d1][d2][2][d3 - 1][d4] * yref[d1][d2][0][d3 - 1][d4];
                        fref[d1][d2][0][d3][d4] -= p.seedling2sapling * yref[d1][d2][2][d3][d4] * yref[d1][d2][0][d3][d4]; // seedling to sapling mean-field forcing
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] -= yref[d1][d2][0][d3][d4] * d3 * p.sapling_death + yref[d1][d2][0][d3][d4] * p.sapling2tree * pow(1.0 * d3, 0.5) * Growth(yref[d1][d2][0][d3][d4 + 1], K); // sapling loss outflux
                        fref[d1][d2][0][d3][d4] += yref[d1][d2][0][d3 + 1][d4] * p.sapling2tree * pow(1.0 + d3, 0.5) * Growth(yref[d1][d2][0][d3 + 1][d4 + 1], K);                                                // sapling loss influx
                        fref[d1][d2][0][d3][d4] += yref[d1][d2][0][d3 + 1][d4] * (d3 + 1) * p.sapling_death;                                                                                                        // sapling death influx

                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4 + 1] = 1.0 * pow(1.0 * d3, 0.5) * p.sapling2tree * Growth(yref[d1][d2][0][d3][d4 + 1], K) - yref[d1][d2][0][d3][d4 + 1] * p.tree_death; // mean-field state
                        else
                            fref[d1][d2][0][d3][d4 + 1] = -yref[d1][d2][0][d3][d4 + 1] * p.tree_death; // mean-field state
                    }

                    if (d3 == p.n_me1 - 3 && d4 == p.n_me2 - 2)
                    { // final saplings master equation state, adult tree mean-field equation
                        fref[d1][d2][0][d3][d4] = 0.0;
                        if (d4 > 0)
                            fref[d1][d2][0][d3][d4] = yref[d1][d2][0][d3 + 1][d4 - 1] * p.sapling2tree * pow(1.0 + d3, 0.5) * Growth(1.0 * d4 - 1.0, K) * Poisson_coupling(d3 + 1, yref[d1][d2][0][d3 + 2][d4 - 1]); // occupation number influx
                        fref[d1][d2][0][d3][d4] -= yref[d1][d2][0][d3][d4] * d4 * p.tree_death * Poisson_coupling(d4, yref[d1][d2][0][d3][d4 + 1]);                                                                   // occupation number outflux
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] += p.seedling2sapling * yref[d1][d2][2][d3 - 1][d4] * yref[d1][d2][0][d3 - 1][d4];
                        fref[d1][d2][0][d3][d4] -= p.seedling2sapling * yref[d1][d2][2][d3][d4] * yref[d1][d2][0][d3][d4]; // seedling to sapling mean-field forcing
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] -= yref[d1][d2][0][d3][d4] * d3 * p.sapling_death + yref[d1][d2][0][d3][d4] * p.sapling2tree * pow(1.0 * d3, 0.5) * Growth(yref[d1][d2][0][d3][d4 + 1], K); // sapling loss outflux
                        fref[d1][d2][0][d3][d4] += yref[d1][d2][0][d3 + 1][d4] * p.sapling2tree * pow(1.0 + d3, 0.5) * Growth(yref[d1][d2][0][d3 + 1][d4 + 1], K) * Poisson_coupling(d3 + 1, yref[d1][d2][0][d3 + 2][d4]);
                        ;                                                                                                                                                  // sapling loss influx
                        fref[d1][d2][0][d3][d4] += yref[d1][d2][0][d3 + 1][d4] * (d3 + 1) * p.sapling_death * Poisson_coupling(d3 + 1, yref[d1][d2][0][d3 + 2][d4]); // sapling death influx

                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4 + 1] = 1.0 * pow(1.0 * d3, 0.5) * p.sapling2tree * Growth(yref[d1][d2][0][d3][d4 + 1], K) - yref[d1][d2][0][d3][d4 + 1] * p.tree_death; // mean-field state
                        else
                            fref[d1][d2][0][d3][d4 + 1] = -yref[d1][d2][0][d3][d4 + 1] * p.tree_death; // mean-field state
                    }

                    if (d3 == p.n_me1 - 2 && d4 == p.n_me2 - 2)
                    { // double mean-regime
                        fref[d1][d2][0][d3][d4] = 0.0;
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] = yref[d1][d2][2][d3 - 1][d4] * p.seedling2sapling * yref[d1][d2][0][d3 - 1][d4]; // sapling MF occupation number influx
                        if (d3 > 0)
                            fref[d1][d2][0][d3][d4] -= (yref[d1][d2][0][d3][d4] * d3 * p.sapling_death + yref[d1][d2][0][d3][d4] * p.sapling2tree * pow(1.0 * d3, 0.5) * Growth(yref[d1][d2][0][d3][d4 + 1], K)) * Poisson_coupling(d3, yref[d1][d2][0][d3 + 1][d4]); // sapling MF occupation number outflux
                        if (d4 > 0 && yref[d1][d2][0][d3 + 1][d4 - 1] > 0)
                            fref[d1][d2][0][d3][d4] += yref[d1][d2][0][d3][d4 - 1] * p.sapling2tree * pow(yref[d1][d2][0][d3 + 1][d4 - 1], 0.5) * Growth(1.0 * d4 - 1.0, K) * (1 - Poisson_coupling(d3, yref[d1][d2][0][d3 + 1][d4 - 1])); // tree MF occupation number influx
                        fref[d1][d2][0][d3][d4] -= yref[d1][d2][0][d3][d4] * d4 * p.tree_death * Poisson_coupling(d4, yref[d1][d2][0][d3][d4 + 1]);                                                                                           // tree MF occupation number outflux

                        if (yref[d1][d2][0][d3 + 1][d4] > 0)
                            fref[d1][d2][0][d3 + 1][d4] = yref[d1][d2][2][d3][d4] * p.seedling2sapling - yref[d1][d2][0][d3 + 1][d4] * p.sapling_death - p.sapling2tree * pow(yref[d1][d2][0][d3 + 1][d4], 0.5) * Growth(yref[d1][d2][0][d3][d4 + 1], K); // sapling mean-field state
                        else
                            fref[d1][d2][0][d3 + 1][d4] = yref[d1][d2][2][d3][d4] * p.seedling2sapling - yref[d1][d2][0][d3 + 1][d4] * (p.sapling_death); // sapling mean-field state
                        if (yref[d1][d2][0][d3 + 1][d4] > 0)
                            fref[d1][d2][0][d3][d4 + 1] = pow(yref[d1][d2][0][d3 + 1][d4], 0.5) * p.sapling2tree * Growth(yref[d1][d2][0][d3][d4 + 1], K) - yref[d1][d2][0][d3][d4 + 1] * p.tree_death; // tree mean-field state
                        else
                            fref[d1][d2][0][d3][d4 + 1] = -yref[d1][d2][0][d3][d4 + 1] * p.tree_death; // tree mean-field state
                    }
                }
            }
        }
    }

    return GSL_SUCCESS;

} //********** end function dydt definition ********************************************************

#endif // DYN_DIFF_HPP_INCLUDED
