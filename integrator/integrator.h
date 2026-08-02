#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <iostream>
#include <cstring>

enum class IntegrationMethod{
    NotSet,
    RK1,
    RK2,
    RK4,
    RK45
};

typedef void (*rhs)(const double* const state, double* const f);

class integrator
{    

public:

    integrator(rhs f_, int n_, IntegrationMethod int_method_) : 
    f(f_), n(n_), int_method(int_method_){};

    integrator(rhs f_, int n_, IntegrationMethod int_method_, double adaptive_local_truncation_error_tolerance_, int adaptive_max_iterations_) :
     f(f_), n(n_), int_method(int_method_), 
     adaptive_max_iterations(adaptive_max_iterations_), 
     adaptive_local_truncation_error_tolerance(adaptive_local_truncation_error_tolerance_){};

    bool integrate(const double delta, double* const state);


private:

    const int n;
    const rhs f;
    const IntegrationMethod int_method = IntegrationMethod::NotSet;
    const double adaptive_local_truncation_error_tolerance = 1e-6;
    const int adaptive_max_iterations = 1e6;

    void multiply_by_const(const double c, double* const u);
    
    void rk1  (const double delta, double* const state);
    void rk2  (const double delta, double* const state);
    void rk4  (const double delta, double* const state);
    bool rk45 (const double delta, double* const state);
    

};


#endif