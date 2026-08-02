#include "geodesics.h"
#include "integrator.h"
#include "globals.h"
#include <iostream>
#include <cmath>
#include <thread> 
#include <atomic>
#include <Eigen/Dense>
#include <GLFW/glfw3.h>
#include <OpenGL/glu.h>
#include <chrono>
#include <cstring>
#include <mutex>
#include <vector>

std::mutex cout_mutex;

void worker(void)
{
    // create a ray
    geodesics ray(globals::solver, globals::delta);        
    ray.init_annulus(globals::annulus_inner_radius, globals::annulus_outer_radius);
    ray.init_Rs(globals::Rs);

    // assign a pixel until none left
    unsigned char pixel_rgb[globals::rgb_length];
    EndDomain end;
    int pixel;
    while ((pixel = globals::render_count.fetch_add(1)) < globals::N_pixels)
    {
        
        // unflatten the coordinates
        const int ll = pixel;
        const int ii = ll % globals::n_pixels;
        const int jj = (ll - ii) / globals::n_pixels;

        // translate coordinates to range [-n/2, n/2]
        const double ii_ = ii - 0.5 * (globals::n_pixels - 1);
        const double jj_ = jj - 0.5 * (globals::n_pixels - 1);

        // calculate direction of camera/ ray
        const Eigen::Vector<double, globals::n_3_vector> point = globals::pixel_sf * (jj_ * globals::u + ii_ * globals::v);
        const Eigen::Vector<double, globals::n_3_vector> ray_direction = (point - globals::n).normalized();
        const Eigen::Vector<double, globals::n_3_vector> d = (globals::cart_to_spherical * ray_direction).normalized();

        // apply transformation
        double p[globals::n_4_vector] = {};
        p[0] = 1.;
        for (int kk = 0; kk < globals::n_3_vector; ++kk)
            p[kk + 1] = globals::tetrad(kk + 1, kk + 1) * d[kk];

        // finialise state
        double IC[globals::n_8_vector] = {0., 
            globals::observer_position[0], globals::observer_position[1], globals::observer_position[2], 
            p[0], p[1], p[2], p[3]
        };

        // integrate geodesic equations to find end domain and rgb
        ray.traverse_geodesic(IC, end, pixel_rgb);

        // illustrate the edge of the black hole
        double r = std::sqrt(ii_ * ii_ + jj_ * jj_) * globals::pixel_sf / globals::Rs;
        if (globals::show_Rs && end != EndDomain::Disk && r < 1.05 && r > 0.95)
            memcpy(&pixel_rgb[0], &globals::red_rgb[0], globals::rgb_length * sizeof(unsigned char));

        // set the pixel colour
        memcpy(&globals::pixels[globals::rgb_length * ll], &pixel_rgb, globals::rgb_length * sizeof(unsigned char));

        // print update
        if ((ll + 1) % (globals::N_pixels / 100) == 0){
            std::lock_guard<std::mutex> lock(cout_mutex);
            std::cout << (((double)ll + 1.) / globals::N_pixels * 100.) << "%.." << std::endl;
        }       
        //printf("%d\n", pixel);

    }

}



int main()
{

    if (!glfwInit()){
        std::cerr << "Failed to initialise GLFW.\n";
        return 0;
    }

    GLFWwindow* window = glfwCreateWindow(globals::n_window,globals::n_window, "BlackHole", NULL, NULL);
    if (!window){
        std::cerr << "Failed to create window.\n";
        glfwTerminate();
        return 0;
    }
    glfwSetWindowAspectRatio(window, 1, 1);  // ensures the objects arent stretched
    glfwMakeContextCurrent(window);

    // threading
    const int n_threads = std::thread::hardware_concurrency();
    std::cout << "Using " << n_threads << " threads." << std::endl;
    std::vector<std::thread> threads;
    for (int ii = 0; ii < n_threads; ++ii)
        threads.emplace_back(worker);

    // run simulation to find all pixel colours
    std::cout << "Generating pixels..." << std::endl;
    for (auto& t: threads)
        t.join();
     std::cout << "Done!" << std::endl;

    // create texture
    GLuint textureID;
    glGenTextures(1, &textureID);
    glBindTexture(GL_TEXTURE_2D, textureID);

    // set texture parameters (important!)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    // upload pixel data to texture
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, globals::n_pixels, globals::n_pixels, 0, GL_RGB, GL_UNSIGNED_BYTE, globals::pixels);

    // create window and display pi
    while (!glfwWindowShouldClose(window)){
        
        // display
        glfwPollEvents();
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        glLoadIdentity();

        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, textureID);

        glBegin(GL_QUADS);
        glTexCoord2f(0.f, 0.f); glVertex2f(-1.f, -1.f);
        glTexCoord2f(1.f, 0.f); glVertex2f(1.f, -1.f);
        glTexCoord2f(1.f, 1.f); glVertex2f(1.f, 1.f);
        glTexCoord2f(0.f, 1.f); glVertex2f(-1.f, 1.f);
        glEnd();

        glDisable(GL_TEXTURE_2D);

        // check if enter key is pressed to close the window
        if (glfwGetKey(window, GLFW_KEY_ENTER) == GLFW_PRESS) {
            std::cout << "Enter key pressed. Closing the window..." << std::endl;
            glfwSetWindowShouldClose(window, GLFW_TRUE);
        }

        glfwSwapBuffers(window);
    }

    glfwDestroyWindow(window);
    glfwTerminate();

    return 1;
}
