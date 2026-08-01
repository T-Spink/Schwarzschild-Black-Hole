#include "integrator.h"

#include <cmath>
#include <cstring>
#include <algorithm>

// ======================================================================================= //
// data

static constexpr double RK45_B45[5][5] = {
    {2./9.   , 0.        , 0.      , 0.     , 0.      },
    {1./12.  , 1./4.     , 0.      , 0.     , 0.      },
    {69./128., -243./128., 135./64., 0.     , 0.      },
    {-17./12., 27./4.    , -27./5. , 16./15., 0.      },
    {65./432., -5./16.   , 13./16. , 4./27. , 5./144. }
};

static constexpr double RK45_c5[6] = {47./450., 0., 12./25., 32./225., 1./30., 6./25.};
static constexpr double RK45_c4[6] = {1./9.   , 0., 9./20. , 16./45. , 1./12., 0.    };


// ======================================================================================= //
// method definitions

void integrator::multiply_by_const(const double c, double* const u)
{
    for (int ii = 0; ii < n; ++ii)
        u[ii] *= c;
};

bool integrator::integrate(const double delta, double* const state)
{
    switch (int_method)
    {
        case IntegrationMethod::RK2: 
            rk2(delta, state);
            return true;

        case IntegrationMethod::RK4: 
            rk4(delta, state);
            return true;
            
        case IntegrationMethod::RK45: 
            return rk45(delta, state);

        case IntegrationMethod::RK1:
        default:
            rk1(delta, state);
            return true;

    }
}

void integrator::rk1(const double delta, double* const state)
{
    double k[n];

    f(state, k);

    for (int ii = 0; ii < n; ++ii)
        state[ii] += delta * k[ii];

}

void integrator::rk2(const double delta, double* const state)
{
    double k1[n], k2[n], u[n];

    // calculate k1
    f(state, &k1[0]);

    // calculate k2
    for (int ii = 0; ii < n; ++ii)
        u[ii] = state[ii] + delta * k1[ii];
    f(&u[0], &k2[0]);   

    // calculate approx.
    for (int ii = 0; ii < n; ++ii)
        state[ii] += (delta / 2.) * ( k1[ii] + k2[ii]);
    
}

void integrator::rk4(const double delta, double* const state)
{
    double u_[n], k1[n], k2[n], k3[n], k4[n];

    // calculate k1
    f(state, &k1[0]);

    // calculate k2
    for (int ii = 0; ii < n; ++ii)
        u_[ii] = state[ii] + (delta / 2.) * k1[ii];
    f(&u_[0], &k2[0]);

    // calculate k3
    for (int ii = 0; ii < n; ++ii)
        u_[ii] = state[ii] + (delta / 2.) * k2[ii];
    f(&u_[0], &k3[0]);

        // calculate k2
    for (int ii = 0; ii < n; ++ii)
        u_[ii] = state[ii] + delta * k3[ii];
    f(&u_[0], &k4[0]);

    // calculate approx.
    for (int ii = 0; ii < n; ++ii)
        state[ii] += (delta / 6.) * ( k1[ii] + 2. * k2[ii] + 2. * k3[ii] + k4[ii]);
    
}

bool integrator::rk45(const double delta0, double* const state)
{
    double delta_ = delta0;

    for (int cc = 0; cc < adaptive_max_iterations; ++cc)
    {        
        double u_[n], k[6][n]; // k is 6 rows of length n state vectors

        // calculate k0
        f(state, &k[0][0]);
        multiply_by_const(delta_, &k[0][0]);

        // calculate k1
        for (int ii = 0; ii < n; ++ii)
            u_[ii] = state[ii] + RK45_B45[0][0] * k[0][ii];
        f(&u_[0], &k[1][0]);
        multiply_by_const(delta_, &k[1][0]);

        // calculate k2
        for (int ii = 0; ii < n; ++ii)
            u_[ii] = state[ii] + RK45_B45[1][0] * k[0][ii] + RK45_B45[1][1] * k[1][ii];
        f(&u_[0], &k[2][0]);
        multiply_by_const(delta_, &k[2][0]);

        // calculate k3
        for (int ii = 0; ii < n; ++ii)
            u_[ii] = state[ii] + RK45_B45[2][0] * k[0][ii] + RK45_B45[2][1] * k[1][ii] + RK45_B45[2][2] * k[2][ii];
        f(&u_[0], &k[3][0]);
        multiply_by_const(delta_, &k[3][0]);

        // calculate k4
        for (int ii = 0; ii < n; ++ii)
            u_[ii] = state[ii] + RK45_B45[3][0] * k[0][ii] + RK45_B45[3][1] * k[1][ii] + RK45_B45[3][2] * k[2][ii] + RK45_B45[3][3] * k[3][ii];
        f(&u_[0], &k[4][0]);
        multiply_by_const(delta_, &k[4][0]);

        // calculate k5
        for (int ii = 0; ii < n; ++ii)
            u_[ii] = state[ii] + RK45_B45[4][0] * k[0][ii] + RK45_B45[4][1] * k[1][ii] + RK45_B45[4][2] * k[2][ii] + RK45_B45[4][3] * k[3][ii] + RK45_B45[4][4] * k[4][ii];
        f(&u_[0], &k[5][0]);
        multiply_by_const(delta_, &k[5][0]);

        // calculate local truncation error
        double local_truncation_error = 0.;
        for (int ii = 0; ii < n; ++ii)
        {
            double temp = 0.;
            for (int jj = 0; jj < 6; ++jj)
                temp += (RK45_c5[jj] - RK45_c4[jj]) * k[jj][ii];
            local_truncation_error += temp * temp;
        }
        local_truncation_error = std::sqrt(local_truncation_error);

        // evaluate error
        if (local_truncation_error > adaptive_local_truncation_error_tolerance){ // error is too large
            delta_ *= 0.9 * std::pow(adaptive_local_truncation_error_tolerance / local_truncation_error, 0.2); // decrease step size
            continue; // repeat step
        }

        // error is acceptable, update state
        for (int ii = 0; ii < n; ++ii)
        {
            for (int jj = 0; jj < 6; ++jj)
                state[ii] += RK45_c5[jj] * k[jj][ii];
        }
        return true;

    }
        
    return false;
}


// ======================================================================================= //