#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <iostream>
#include <cstring>

enum IntegrationMethod{
    RK1,
    RK2,
    RK4,
    RK45
};

template<typename T> class integrator
{    
    typedef void (T::*rhs)(const double* const state, double* const f);

public:

    integrator(T* obj_, rhs f_, int n_):obj(obj_),f(f_),n(n_){};
    void init_rk45_tol(const double tol){ rk45_tol = tol; };
    void init_integration(const IntegrationMethod int_method_){ int_method=int_method_; };

    bool integrate(const double delta, double* const state);

private:
    int n;
    T* obj;
    rhs f;
    IntegrationMethod int_method;
    double rk45_tol;
    const int max_rk45its = 1e5;

    void rk1      (const double delta, double* const state);
    void rk2      (const double delta, double* const state);
    void rk4      (const double delta, double* const state);
    bool rk45     (const double delta, double* const state);
    
};

#include "integrator.tpp"

#endif