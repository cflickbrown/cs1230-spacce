#include "blackhole.h"

// The algorithm seen here is primarily taken from
// https://20k.github.io/c++/2024/05/31/schwarzschild.html

BlackHole::BlackHole() {}

static float getSchwarzchildRadius(){

}

glm::mat4 BlackHole::getKerrMetric(float r, float theta, float phi, float t){

}

glm::mat4 BlackHole::getSchwarzchildMetric(float r, float theta, float phi, float t){
    // This is independant of phi or t but those numbers are here to stay
    return glm::mat4{
        -(1-rs/r), 0,0,0,
        0, 1/(1-rs/r),0,0,
        0,0,r*r,0,
        0,0,0,r*r * sin(theta)*sin(theta)
    };
}

glm::mat4 BlackHole::getTetradBasis(float r, float theta, float phi, float t){
    glm::mat4 tetradSet = {
        1/(1 - sqrt(rs/r)), 0, 0, 0,
        0, sqrt(1 - (rs/r)), 0, 0,
        0, 0, rs/r, 0,
        0, 0, 0, 1/(r*sin(theta))
    };
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
std::vector<glm::mat4> BlackHole::calculateChristoffel(glm::vec4 position, glm::mat4 metric){
    // get the metric
    glm::mat4 metInv = glm::inverse(metric);

    std::vector<glm::mat4> gamma;
    gamma.assign(4,glm::mat4(0));

    std::vector<glm::mat4> metricDifference;
    metricDifference.assign(4,glm::mat4(0));

    // Calculate the metric difference

    // calculate the christoffel symbols per element

    // Equation is:
    // Gamma^mu_alpha,beta = 1/2 * sum sigma from 0->3 (g^mu,sigma * ())

    // return
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
