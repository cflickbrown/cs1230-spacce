#pragma once

#include <glm/glm.hpp>
#include <random>

#include "utils/sceneparser.h"
#include "utils/perlinnoise.h"

namespace StarGenerator {

void generateStars(RenderData& renderData);

void generatePosition(std::vector<glm::vec3>& positions, std::mt19937& generator, PerlinNoise noise);

float generateSize(std::mt19937& generator, std::normal_distribution<float>& size);

glm::vec3 generateColor(std::mt19937& generator, std::uniform_int_distribution<>& color);

glm::vec2 generateIntensity(std::mt19937& generator, PerlinNoise noise, glm::vec3 position);

glm::vec2 getIntensity(std::mt19937& generator, float intensityMean, float intensitySD, float shineMean, float shineSD);

float lighten(const float color);

void generatePrimitive(std::mt19937& generator, PerlinNoise noise, RenderData& renderData, std::vector<glm::vec3>& positions);

}
