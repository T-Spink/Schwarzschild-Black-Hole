#include <cmath>
#include <cstring>
#include <algorithm>

static const float RK45_B45[5][5] = {
    {2.f/9.f   , 0.f         , 0.f       , 0.f      , 0.f       },
    {1.f/12.f  , 1.f/4.f     , 0.f       , 0.f      , 0.f       },
    {69.f/128.f, -243.f/128.f, 135.f/64.f, 0.f      , 0.f       },
    {-17.f/12.f, 27.f/4.f    , -27.f/5.f , 16.f/15.f, 0.f       },
    {65.f/432.f, -5.f/16.f   , 13.f/16.f , 4.f/27.f , 5.f/144.f }
};

static const float RK45_c5[6] = {47.f/450.f, 0.f, 12.f/25.f, 32.f/225.f, 1.f/30.f, 6.f/25.f};
static const float RK45_c4[6] = {1.f/9.f   , 0.f, 9.f/20.f , 16.f/45.f , 1.f/12.f, 0.f     };

template<typename T> 
bool integrator<T>::integrate(const double delta, double* const state)
{
    switch (int_method)
    {
        case RK2: 
            rk2(delta, state);
            return true;

        case RK4: 
            rk4(delta, state);
            return true;
            
        case RK45: 
            return rk45(delta, state);

        case RK1:
        default:
            rk1(delta, state);
            return true;

    }
}

template<typename T>
void integrator<T>::rk1(const double delta, double* const state)
{
    double k[n];

    (obj->*f)(state, &k[0]);

    for (int ii=0; ii<n; ii++)
        state[ii] += delta * k[ii];

}

template<typename T>
void integrator<T>::rk2(const double delta, double* const state)
{
    double k1[n], k2[n];

    // calculate k1
    (obj->*f)(state, &k1[0]);

    // calculate k2
    double u_[n];
    for (int ii=0; ii<n; ii++)
        u_[ii] = state[ii] + delta * k1[ii];
    (obj->*f)(&u_[0], &k2[0]);

    // calculate approx.
    for (int ii=0; ii<n; ii++)
        state[ii] += (delta/2) * ( k1[ii] + k2[ii]);
    
}

template<typename T>
void integrator<T>::rk4(const double delta, double* const state)
{
    double u_[n], k1[n], k2[n], k3[n], k4[n];

    // calculate k1
    (obj->*f)(state, &k1[0]);

    // calculate k2
    for (int ii=0; ii<n; ii++)
        u_[ii] = state[ii] + (delta/2.f) * k1[ii];
    (obj->*f)(&u_[0], &k2[0]);

    // calculate k3
    for (int ii=0; ii<n; ii++)
        u_[ii] = state[ii] + (delta/2.f) * k2[ii];
    (obj->*f)(&u_[0], &k3[0]);

        // calculate k2
    for (int ii=0; ii<n; ii++)
        u_[ii] = state[ii] + delta * k3[ii];
    (obj->*f)(&u_[0], &k4[0]);

    // calculate approx.
    for (int ii=0; ii<n; ii++)
        state[ii] += (delta/6) * ( k1[ii] + 2*k2[ii] + 2*k3[ii] + k4[ii]);
    
}

template<typename T>
bool integrator<T>::rk45(const double delta0, double* const state)
{
    double delta_ = delta0;

    int cc=0;
    for (;;)
    {        

        // check for max adaptations
        if (cc++>max_rk45its) return false;

        double u_[n], k[6][n]; // k is 6 rows of length n state vectors

        // calculate k0
        (obj->*f)(state, &k[0][0]);
        for (int ii=0; ii<n; ii++) k[0][ii] *= delta_;

        // calculate k1
        for (int ii=0; ii<n; ii++)
            u_[ii] = state[ii] + RK45_B45[0][0] * k[0][ii];
        (obj->*f)(&u_[0], &k[1][0]);
        for (int ii=0; ii<n; ii++) k[1][ii] *= delta_;

        // calculate k2
        for (int ii=0; ii<n; ii++)
            u_[ii] = state[ii] + RK45_B45[1][0] * k[0][ii] + RK45_B45[1][1] * k[1][ii];
        (obj->*f)(&u_[0], &k[2][0]);
        for (int ii=0; ii<n; ii++) k[2][ii] *= delta_;

        // calculate k3
        for (int ii=0; ii<n; ii++)
            u_[ii] = state[ii] + RK45_B45[2][0] * k[0][ii] + RK45_B45[2][1] * k[1][ii] + RK45_B45[2][2] * k[2][ii];
        (obj->*f)(&u_[0], &k[3][0]);
        for (int ii=0; ii<n; ii++) k[3][ii] *= delta_;

        // calculate k4
        for (int ii=0; ii<n; ii++)
            u_[ii] = state[ii] + RK45_B45[3][0] * k[0][ii] + RK45_B45[3][1] * k[1][ii] + RK45_B45[3][2] * k[2][ii] + RK45_B45[3][3] * k[3][ii];
        (obj->*f)(&u_[0], &k[4][0]);
        for (int ii=0; ii<n; ii++) k[4][ii] *= delta_;

        // calculate k5
        for (int ii=0; ii<n; ii++)
            u_[ii] = state[ii] + RK45_B45[4][0] * k[0][ii] + RK45_B45[4][1] * k[1][ii] + RK45_B45[4][2] * k[2][ii] + RK45_B45[4][3] * k[3][ii] + RK45_B45[4][4] * k[4][ii];
        (obj->*f)(&u_[0], &k[5][0]);
        for (int ii=0; ii<n; ii++) k[5][ii] *= delta_;

        // calculate local truncation error
        double TE=0.;
        for (int ii=0; ii<n; ii++){
            double TE_ = 0.;
            for (int jj=0; jj<6; jj++)
                TE_ += (RK45_c5[jj]-RK45_c4[jj])*k[jj][ii];
            TE += TE_*TE_;
        }
        TE = sqrt(TE);

        // evaluate error
        if (TE>rk45_tol){ // error is too large
            delta_ *= 0.9*pow(rk45_tol/TE, 0.2); // decrease step size
            continue; // repeat step
        }

        // error is acceptable, calculate new step 
        for (int ii=0; ii<n; ii++){
            for (int jj=0; jj<6; jj++)
                state[ii] += RK45_c5[jj] * k[jj][ii];
        }
        return true;

    }
        
}
