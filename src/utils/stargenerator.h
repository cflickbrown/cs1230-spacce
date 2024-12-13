#pragma once

#include <glm/glm.hpp>
#include <random>

#include "utils/sceneparser.h"
#include "utils/perlinnoise.h"

namespace StarGenerator {

// Called once in main to generate stars and add them to the renderData.
void generateStars(RenderData& renderData);

// Generates positions for all stars to be added in the scene using a combination of stratfied sampling and
// Perlin Noise-based rejection sampling to create a Perlin Noise density position generation.
void generatePosition(std::vector<glm::vec3>& positions, std::mt19937& generator, PerlinNoise noise);

// Uses a Gaussian sampling to generate a size for a star.
float generateSize(std::mt19937& generator, std::normal_distribution<float>& size);

// Uses uniformly randomness to choose one of five base colors.
glm::vec3 generateColor(std::mt19937& generator, std::uniform_int_distribution<>& color);

// Uses Perlin Noise to generate an intensity and shininess for the star based on the position give.
// Classifies stars based on their Perlin Noise value as dim, medium, or bright, then samples intensity/shininess
// values using Gaussian sampling with means/std-devs based on their classification.
glm::vec2 generateIntensity(std::mt19937& generator, PerlinNoise noise, glm::vec3 position);

// Returns the intensity based on the values given.
glm::vec2 getIntensity(std::mt19937& generator, float intensityMean, float intensitySD, float shineMean, float shineSD);

// Creates the actual star and puts it into renderData along with the volumetric "halo".
void generatePrimitive(std::mt19937& generator, PerlinNoise noise, RenderData& renderData, std::vector<glm::vec3>& positions);

}
