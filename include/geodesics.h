#ifndef GEODESICS_H
#define GEODESICS_H

#include "integrator.h"
#include <cstring>

enum EndDomain{
        Error,
        BlackHole,
        Disk,
        Space
    };

class geodesics
{

public:

    geodesics(void);
    void init_state(const double (&u)[8]);
    void init_Rs(const double Rs_){ Rs = Rs_;};
    void init_disk(const double R0, const double Rf);
    void init_solver(const double delta, const IntegrationMethod int_method, const double tol);
    const EndDomain GetEnd(void){return end;};
    void reset(void);
    ~geodesics(void);

    void GetRGB(const double (&IC)[8], unsigned char (&rgb)[3]);

private:

    double state[8];  // t,r,theta,phi,dt,dR,dtheta,dphi

    EndDomain end;

    double g[4], ginv[4];
    double dg[4][4][4];
    
    double initial_state[8];
    double Rs, disk_R0, disk_Rf;
    const unsigned char rgb0[3] = {255, 230, 200}; // colour at disk0
    const unsigned char rgbf[3] = {80, 30, 10};  // colour at diskf
    double delta;

    integrator<geodesics>* solver;

    EndDomain FindEndDomain(const double (&IC)[8]);

    void rhs(const double* const u, double* const f);

    void calc_metric(const double (&state)[8]);
    void calc_dmetric(const double (&state)[8]);
    void calc_metric_inv(const double (&state)[8]);
};

#endif