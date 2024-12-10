#include "perlinnoise.h"

PerlinNoise::PerlinNoise() {

    // Task 8: turn off wireframe shading
    // m_wireshade = true; // STENCIL CODE
    m_wireshade = false; // TA SOLUTION

    // Define resolution of terrain generation
    m_resolution = 200;

    // Generate random vector lookup table
    m_lookupSize = 2000;
    m_randVecLookup.reserve(m_lookupSize);

    // Initialize random number generator
    std::srand(20);

    // Populate random vector lookup table
    for (int i = 0; i < m_lookupSize; i++)
    {
        m_randVecLookup.push_back(glm::vec3(std::rand() * 2.0 / RAND_MAX - 1.0,
                                            std::rand() * 2.0 / RAND_MAX - 1.0,
                                            std::rand() * 2.0 / RAND_MAX - 1.0));
    }
}

// Destructor
PerlinNoise::~PerlinNoise()
{
    m_randVecLookup.clear();
}

// Samples the (infinite) random vector grid at (row, col)
glm::vec3 PerlinNoise::sampleRandomVector(int x, int y, int z)
{
    std::hash<int> intHash;
    int index = intHash(x * 41 + y * 43 + z * 43) % m_lookupSize;
    return m_randVecLookup.at(index);
}

// Helper for computePerlin() and, possibly, getColor()
float interpolate(float A, float B, float alpha) {
    // Task 4: implement your easing/interpolation function below

    float easedAlpha = (3.0 * alpha * alpha) - (2.0 * alpha * alpha * alpha);
    return A + (easedAlpha * (B - A));

    // Return 0 as placeholder
}

// Computes the intensity of Perlin noise at some point
float PerlinNoise::computePerlin3d(float x, float y, float z) {
    // Task 1: get grid indices (as ints)

    int zInt = floor(z);
    int yInt = floor(y);
    int xInt = floor(x);

    // Task 2: compute offset vectors
    glm::vec3 gOffset1 = glm::vec3(x - xInt, y - yInt, z - zInt);
    glm::vec3 gOffset2 = glm::vec3(x - (xInt + 1), y - yInt, z - zInt);
    glm::vec3 gOffset3 = glm::vec3(x - xInt, y - (yInt + 1), z - zInt);
    glm::vec3 gOffset4 = glm::vec3(x - xInt, y - yInt, z - (zInt  + 1));
    glm::vec3 gOffset5 = glm::vec3(x - (xInt + 1), y - (yInt + 1), z - zInt);
    glm::vec3 gOffset6 = glm::vec3(x - (xInt + 1), y - yInt, z - (zInt + 1));
    glm::vec3 gOffset7 = glm::vec3(x - xInt, y - (yInt + 1), z - (zInt + 1));
    glm::vec3 gOffset8 = glm::vec3(x - (xInt + 1), y - (yInt + 1), z - (zInt + 1));

    // Task 3: compute the dot product between the grid point direction vectors and its offset vectors
    float A = dot(sampleRandomVector(xInt, yInt, zInt), gOffset1); // dot product between top-left direction and its offset
    float B = dot(sampleRandomVector(xInt + 1, yInt, zInt), gOffset2); // dot product between top-right direction and its offset
    float C = dot(sampleRandomVector(xInt, yInt + 1, zInt), gOffset3); // dot product between bottom-right direction and its offset
    float D = dot(sampleRandomVector(xInt + 1, yInt + 1, zInt), gOffset4); // dot product between bottom-left direction and its offset

    float E = dot(sampleRandomVector(xInt, yInt, zInt + 1), gOffset5); // dot product between top-left direction and its offset
    float F = dot(sampleRandomVector(xInt + 1, yInt, zInt + 1), gOffset6); // dot product between top-right direction and its offset
    float G = dot(sampleRandomVector(xInt, yInt + 1, zInt + 1), gOffset7); // dot product between bottom-right direction and its offset
    float H = dot(sampleRandomVector(xInt + 1, yInt + 1, zInt + 1), gOffset8); // dot product between bottom-left direction and its offset

    // Task 5: Debug this line to properly use your interpolation function to produce the correct value

    return interpolate(
            interpolate(interpolate(G, H, x - xInt), interpolate(E, F, x - xInt), y - yInt),
            interpolate(interpolate(C, D, x - xInt), interpolate(A, B, x - xInt), y - yInt),
            z - zInt);
    //change: C, D flipped, 0.5 -> dx or dy

    // Return 0 as a placeholder
    // return 0;
}
