#include "integrator.h"
#include "geodesics.h"
#include "globals.h"
#include <Eigen/Dense>
#include <numbers>

namespace globals
{

    const double pixel_sf = observer_R * std::tan(fov_degs * rads_per_degs / 2.) / (n_pixels / 2.);
    const Eigen::Vector<double, n_3_vector> observer_position = {observer_R, observer_theta, observer_phi}; 

    static const double st = std::sin(observer_theta);
    static const double ct = std::cos(observer_theta);
    static const double sp = std::sin(observer_phi);
    static const double cp = std::cos(observer_phi);

    // determine normal to and basis vectors of the image plane 
    const Eigen::Vector<double, n_3_vector> n = {            
        observer_R * st * cp, 
        observer_R * st * sp,
        observer_R * ct
    };
    const Eigen::Vector<double, n_3_vector> u = {-n[2] / n.norm(), 0, n[0] / n.norm()};
    const Eigen::Vector<double, n_3_vector> v = n.cross(u).normalized();

    static constexpr double f = 1. - Rs / observer_R;
    const Eigen::Matrix<double, n_4_vector, n_4_vector> tetrad{
            {1. / std::sqrt(f), 0.          , 0.             , 0.                    },
            {0.               , std::sqrt(f), 0.             , 0.                    },
            {0.               , 0.          , 1. / observer_R, 0.                    },
            {0.               , 0.          , 0.             , 1. / (observer_R * st)}
        };
    const Eigen::Matrix<double, n_3_vector, n_3_vector> cart_to_spherical{
            { (st * cp), (st * sp), (ct     ) },
            { (ct * cp), (ct * sp), (-st    ) },
            { (-sp    ), (cp     ), (0.     ) }
    };


    // ========================================= // 
    // variables

    unsigned char pixels_array[n_pixels * n_pixels * rgb_length] = {};
    unsigned char* const pixels = &pixels_array[0];
    std::atomic<int> render_count = 0;

    // numerical solver - can be used wherever these arguments are constant and includes no state
    integrator solver_(&f_geodesic, n_8_vector, int_method, rk45_tolerance, rk45_max_steps);
    integrator * const solver = &solver_;

    // ========================================= // 
}
