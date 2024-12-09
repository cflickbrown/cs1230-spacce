#ifndef BLACKHOLE_H
#define BLACKHOLE_H

#include <glm/glm.hpp>
#include <vector>

class BlackHole
{
public:
    BlackHole();

    // One and done static methods

    static float getSchwarzchildRadius();

    glm::mat4 getTetradBasis(float r, float theta, float phi, float t);

    // Bit of a hack here but if we can approximate a 4x4x4 tensor with an array of 4x4 mats...
    std::vector<glm::mat4> calculateChristoffel(glm::vec4 position, auto&& metric);

    glm::mat4 getSchwarzchildMetric(float r, float theta, float phi, float t);

    glm::mat4 getKerrMetric(float r, float theta, float phi, float t);

    void renderBlackHole(float t, float r, float theta, float phi);

    /**
     * @brief getBoyerLindquistCoordinate Converts a cartesian coordinate to Boyer Lindquist coordinates.
     * Less expensive than Kerr Schild coordinates but cannot be used within the outer event horizon.
     * @param cartesianCoordinate An x,y,z,n homogenous cartesian coordinate
     * @return A Boyer Lindquist coordinate
     */
    glm::vec4 getBoyerLindquistCoordinate(glm::vec4& cartesianCoordinate);

    /**
     * @brief getKerrSchildCoordinate Converts a cartesian coordinate into Kerr Schild coordinates.
     * Expensive. Use only within the outer event horizon.
     * @param cartesianCoordinate An x,y,z,n homogenous cartesian coordinate
     * @return A Kerr Schild coordinate
     */
    glm::vec4 getKerrSchildCoordinate(glm::vec4& cartesianCoordinate);

    /**
     * @brief stepLightlikeGeodesic
     * @param position
     * @param direction
     * @return A homogenized direction vector for usage in standard spacetime-agnostic world-space coordinates.
     * This contains the exact direction and magnitude of the step.
     */
    glm::vec4 stepLightlikeGeodesic(glm::vec4 position, glm::vec3 direction);


    // Getters/Setters

    float getSpin();
    float getMass();

    /**
     * @brief setSpin Sets the spin parameter of the black hole.
     * @param spin A spin parameter∈ (-1,1). Parameters outside these values will result in naked singularities.
     */
    void setSpin(float spin);

    /**
     * @brief setMass Sets the mass of the black hole.
     * @param mass A value in terms of M, where M equals 200,000 solar masses.
     */
    void setMass(float mass);

    glm::vec4 position;

    float rs;   // Schwarschild radius

    float M;
    float J;
    float Q; // Most likely unused here, as we're not going for a Kerr-Newman black hole

};

#endif // BLACKHOLE_H
