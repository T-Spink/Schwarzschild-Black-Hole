#include "geodesics.h"
#include "constants.h"
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

void render_pixel(const globals& globs, geodesics& ray, const int ll)
{
     // unflatten the coordinates
    const int ii = ll%globs.n_pixels;
    const int jj = (ll-ii) / globs.n_pixels;

    // translate coordinates to range [-n/2, n/2]
    const double ii_ = ii - 0.5*(globs.n_pixels-1);
    const double jj_ = jj - 0.5*(globs.n_pixels-1);

    // calculate direction of camera/ ray
    const Eigen::Vector<double,3> point = globs.pixel_sf*(jj_*globs.u + ii_*globs.v);
    const Eigen::Vector<double,3> ray_direction = (point - globs.n).normalized();
    const Eigen::Vector<double,3> d = (globs.cart_to_spherical*ray_direction).normalized();

    // apply transformation
    double p[4] = {};
    p[0] = 1.;
    for (int kk=0; kk<3; kk++)
        p[kk+1] = globs.tetrad(kk+1,kk+1)*d[kk];

    // finialise state
    double IC[8] = {0.};
    memcpy(&IC[1], globs.position.data(), 3*sizeof(double));  // initialise pos
    memcpy(&IC[4], p                    , 4*sizeof(double));  // initialise momentum

    // integrate geodesic equations to find end domain and rgb 
    unsigned char rgb[3];
    ray.GetRGB(IC, rgb);

    // illustrate the edge of the black hole
    EndDomain end = ray.GetEnd();
    double r = std::sqrt(ii_*ii_ + jj_*jj_)*globs.pixel_sf/globs.Rs;
    if (globs.show_Rs && end!=Disk && r<1.05 && r>0.95){
        rgb[0] = 255;
        rgb[1] = 0.;
        rgb[2] = 0.;
    }

    // set the pixel colour
    memcpy(&globs.pixels[3*ll], &rgb, 3*sizeof(unsigned char));

}

void worker(globals& globs)
{
    // create a ray
    geodesics ray;        
    ray.init_disk(globs.disk0, globs.diskf);
    ray.init_Rs(globs.Rs);
    ray.init_solver(globs.delta, globs.int_method, globs.rk45_tol);

    // assign a pixel until none left
    int pixel=-1;
    while ((pixel=globs.render_count.fetch_add(1))<globs.N_pixels){

        // render the next pixel
        render_pixel(globs, ray, pixel);

        // reset the ray
        ray.reset();

        // print update
        if ((pixel+1)%(globs.N_pixels/100)==0){
            std::lock_guard<std::mutex> lock(cout_mutex);
            std::cout << (((double)pixel+1.)/globs.N_pixels*100.) << "%..." << std::endl;
        }       
    }
}


int main()
{
    // structure containing all global variables
    globals globs;

    if (!glfwInit()){
        std::cerr << "Failed to initialise GLFW.\n";
        return -1;
    }

    GLFWwindow* window = glfwCreateWindow(globs.n_window,globs.n_window, "BlackHole", NULL, NULL);
    if (!window){
        std::cerr << "Failed to create window.\n";
        glfwTerminate();
        return -1;
    }
    glfwSetWindowAspectRatio(window, 1, 1);  // ensures the objects arent stretched
    glfwMakeContextCurrent(window);

    // threading
    const int n_threads = std::thread::hardware_concurrency();
    std::cout << "Using " << n_threads << " threads." << std::endl;
    std::vector<std::thread> threads;
    for (int ii=0; ii<n_threads; ii++)
        threads.emplace_back(worker, std::ref(globs));

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
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, globs.n_pixels, globs.n_pixels, 0, GL_RGB, GL_UNSIGNED_BYTE, globs.pixels);

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

        // Check if Enter key is pressed to close the window
        if (glfwGetKey(window, GLFW_KEY_ENTER) == GLFW_PRESS) {
            std::cout << "Enter key pressed. Closing the window..." << std::endl;
            glfwSetWindowShouldClose(window, GLFW_TRUE);
        }

        glfwSwapBuffers(window);
    }

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}
