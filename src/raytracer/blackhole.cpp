#include "blackhole.h"

// The algorithm seen here is primarily taken from
// https://20k.github.io/c++/2024/05/31/schwarzschild.html

BlackHole::BlackHole() {}

static float getSchwarzchildRadius(){

}

glm::mat4 BlackHole::getKerrMetric(float t, float r, float theta, float phi){

}

glm::mat4 BlackHole::getSchwarzchildMetric(float t, float r, float theta, float phi){
    // This is independant of phi or t but those numbers are here to stay
    return glm::mat4{
        -(1-rs/r), 0,0,0,
        0, 1/(1-rs/r),0,0,
        0,0,r*r,0,
        0,0,0,r*r * sin(theta)*sin(theta)
    };
}

glm::mat4 BlackHole::getTetradBasis(float t, float r, float theta, float phi){
    glm::mat4 tetradSet = {
        1/(1 - sqrt(rs/r)), 0, 0, 0,
        0, sqrt(1 - (rs/r)), 0, 0,
        0, 0, rs/r, 0,
        0, 0, 0, 1/(r*sin(theta))
    };
    return tetradSet;
}

auto getGradient(auto&& function, glm::vec4 position, int dimension){
    float epsilon = 0.0000001;
    glm::vec4 delta = position;

    delta[dimension] += epsilon;

    return (function(delta) - function(position))/epsilon;
}

// There isn't a really good way of explaining this in the papers, as it is assumed as common
// knowledge for the type of person for this type of paper. This is the one part that has tripped
// me up every single time.

/**
 * @brief BlackHole::calculateChristoffel Calculates the christoffel symbol of the second kind for a specific position and metric
 * @param position
 * @param metric
 * @return
 */
std::vector<glm::mat4> BlackHole::calculateSchwarzchildChristoffel(glm::vec4 position){
    // get the metric
    //glm::mat4 metInv = glm::inverse(metric);

    std::vector<glm::mat4> gamma;
    gamma.assign(4,glm::mat4(0));

    std::vector<glm::mat4> metricDifference;
    metricDifference.assign(4,glm::mat4(0));

    // Equation is:
    // Gamma^mu_alpha,beta = 1/2 * sum sigma from 0->3 (g^mu,sigma * ())

    // Calculate the metric difference - Difference in overall force, therefore curvature of spacetime

    // for(int i = 0; i < 4; ++i){
    //     glm::mat4 gradient = getGradient(metric, position, i);
    //     // for(int j = 0; j < 4; ++j){
    //     //     for(int k = 0; k < 4; ++k){
    //     //         metricDifference[i][j][k] = gradient[j][k];
    //     //     }
    //     // }
    //     metricDifference[i] = gradient;
    // }

    // calculate the christoffel symbols per element
    // THIS: https://physics.stackexchange.com/questions/733433/christoffel-symbols-for-schwarzschild-metric

    // Due to math reasons about t not interacting with this, the only symbols I need to worry about are the following 9:
    // gamma_x_y_z = gamma_x_z_y
    float r = position[1];
    float theta = position[2];
    float phi = position[3];
    float gamma_r_t_t = (M*(r - 2*M))/(pow(r,3));
    float gamma_r_r_r = -(M)/(r*(r-2*M));
    float gamma_r_theta_theta = -(r - 2*M);
    float gamma_r_phi_phi = -(r - 2*M)*(sin(theta)*sin(theta));
    float gamma_t_r_t = M/(r*(r - 2*M));
    float gamma_theta_r_theta = 1/r;
    float gamma_theta_phi_phi = -sin(theta)*cos(theta);
    float gamma_phi_r_phi = 1/r;
    float gamma_phi_theta_phi = cos(theta)/sin(theta);

    // Gamma^t_XX
    gamma[0] = glm::mat4{
        0, gamma_t_r_t,0,0,
        gamma_t_r_t,0,0,0,
        0,0,0,0,
        0,0,0,0
    };

    // Gamma^r_XX
    gamma[1] = glm::mat4{
        gamma_r_t_t,0,0,0,
        0,gamma_r_r_r,0,0,
        0,0,gamma_r_theta_theta,0,
        0,0,0,gamma_r_phi_phi
    };

    // Gamma^theta_XX
    gamma[2] = glm::mat4{
        0,0,0,0,
        0,0,gamma_theta_r_theta,0,
        0,gamma_theta_r_theta,0,0,
        0,0,0,gamma_theta_phi_phi
    };

    // Gamma^phi_XX
    gamma[3] = glm::mat4{
        0,0,0,0,
        0,0,0,gamma_phi_r_phi,
        0,0,0,gamma_phi_theta_phi,
        0,gamma_phi_r_phi,gamma_phi_theta_phi,0
    };

    return gamma;
}

/**
 * @brief BlackHole::calculateChristoffel Calculates the christoffel symbol of the second kind for a specific position and metric
 * @param position
 * @param metric
 * @return
 */
std::vector<glm::mat4> BlackHole::calculateKerrChristoffel(glm::vec4 position){
    // get the metric
    //glm::mat4 metInv = glm::inverse(metric);

    std::vector<glm::mat4> gamma;
    gamma.assign(4,glm::mat4(0));

    std::vector<glm::mat4> metricDifference;
    metricDifference.assign(4,glm::mat4(0));

    // Equation is:
    // Gamma^mu_alpha,beta = 1/2 * sum sigma from 0->3 (g^mu,sigma * ())

    // Calculate the metric difference - Difference in overall force, therefore curvature of spacetime

    // for(int i = 0; i < 4; ++i){
    //     glm::mat4 gradient = getGradient(metric, position, i);
    //     // for(int j = 0; j < 4; ++j){
    //     //     for(int k = 0; k < 4; ++k){
    //     //         metricDifference[i][j][k] = gradient[j][k];
    //     //     }
    //     // }
    //     metricDifference[i] = gradient;
    // }

    // calculate the christoffel symbols per element
    // THIS: https://physics.stackexchange.com/questions/733433/christoffel-symbols-for-schwarzschild-metric

    // Due to math reasons about t not interacting with this, the only symbols I need to worry about are the following 9:
    // gamma_x_y_z = gamma_x_z_y
    float r = position[1];
    float theta = position[2];
    float phi = position[3];
    float gamma_r_t_t = (M*(r - 2*M))/(pow(r,3));
    float gamma_r_r_r = -(M)/(r*(r-2*M));
    float gamma_r_theta_theta = -(r - 2*M);
    float gamma_r_phi_phi = -(r - 2*M)*(sin(theta)*sin(theta));
    float gamma_t_r_t = M/(r*(r - 2*M));
    float gamma_theta_r_theta = 1/r;
    float gamma_theta_phi_phi = -sin(theta)*cos(theta);
    float gamma_phi_r_phi = 1/r;
    float gamma_phi_theta_phi = cos(theta)/sin(theta);

    // Gamma^t_XX
    gamma[0] = glm::mat4{
        0, gamma_t_r_t,0,0,
        gamma_t_r_t,0,0,0,
        0,0,0,0,
        0,0,0,0
    };

    // Gamma^r_XX
    gamma[1] = glm::mat4{
        gamma_r_t_t,0,0,0,
        0,gamma_r_r_r,0,0,
        0,0,gamma_r_theta_theta,0,
        0,0,0,gamma_r_phi_phi
    };

    // Gamma^theta_XX
    gamma[2] = glm::mat4{
        0,0,0,0,
        0,0,gamma_theta_r_theta,0,
        0,gamma_theta_r_theta,0,0,
        0,0,0,gamma_theta_phi_phi
    };

    // Gamma^phi_XX
    gamma[3] = glm::mat4{
        0,0,0,0,
        0,0,0,gamma_phi_r_phi,
        0,0,0,gamma_phi_theta_phi,
        0,gamma_phi_r_phi,gamma_phi_theta_phi,0
    };

    return gamma;
}

/**
     * @brief getBoyerLindquistCoordinate Converts a cartesian coordinate to Boyer Lindquist coordinates.
     * Less expensive than Kerr Schild coordinates but cannot be used within the outer event horizon.
     * @param cartesianCoordinate An x,y,z,n homogenous cartesian coordinate
     * @return A Boyer Lindquist coordinate
     */
glm::vec4 BlackHole::getBoyerLindquistCoordinate(glm::vec4& cartesianCoordinate){

}

/**
     * @brief getKerrSchildCoordinate Converts a cartesian coordinate into Kerr Schild coordinates.
     * Expensive. Use only within the outer event horizon.
     * @param cartesianCoordinate An x,y,z,n homogenous cartesian coordinate
     * @return A Kerr Schild coordinate
     */
glm::vec4 BlackHole::getKerrSchildCoordinate(glm::vec4& cartesianCoordinate){

}

glm::vec4 BlackHole::stepLightlikeGeodesic(glm::vec4 position, glm::vec3 direction){
    // We assume that steps 1 and 2 of the pipeline have been completed upon construction.
    // We assume that step 3 has been completed in raytracer.cpp
    // Step 4, construct curved spacetime

    glm::vec3 minkowskiDirection = {-direction[2], direction[1], direction[0]};

    //glm::vec4 velocity = tetradSet * glm::vec4(-1,direction);
}

    /**
    Current implementation issue. Currently trying to find a way to interfere as little
    as possible in the raymarcher. The current plan is to make a function that steps once across
    a lightlike geodesic, but that imposes significant overhead from the constant conversion
    to and from world space / spherical / locally flat minkowski metric.
    **/
void BlackHole::renderBlackHole(float t, float r, float theta, float phi){
    // Step 1: Calculate metric tensor g_mu,v

    glm::mat4 metricTensor = getSchwarzchildMetric(t,r,theta,phi);

    // Step 2: calculate tetrad basis vectors
    glm::mat4 tetradSet = {
        1/(1 - sqrt(rs/r)), 0, 0, 0,
        0, sqrt(1 - (rs/r)), 0, 0,
        0, 0, rs/r, 0,
        0, 0, 0, 1/(r*sin(theta))
    };

    // Step 3: Construct initial ray direction

    // 3a: pick direction d^k, with time component (v0) set to -1, spatial is set to d

    // For now, we're just going to shoot a dummy vector in the direction -1,1,1,0

    glm::vec4 dk = {-1,1,1,0};

    // Step 4: Construct curved spacetime -> v_curved = e^u_i * v^i_flat

    glm::vec4 position = {10,0,0,0};
    glm::vec3 direction = {-dk[2], dk[1], dk[0]};

    glm::vec4 velocity = tetradSet * glm::vec4(-1,direction);

    // Step 5: Plug into metric tensor. Step a certain amount


    // Step 6: Run a set number of steps,

    // Step 6a: If inside EH, we're done. Pixel is black

    // Step 6b: Otherwise, we have a ray direction and position, render as normal;


}


// Getters/Setters

float BlackHole::getSpin(){
    return this->J;
}

float BlackHole::getMass(){
    return this->M;
}
