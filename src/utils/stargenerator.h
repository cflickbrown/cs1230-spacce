#pragma once

#include <glm/glm.hpp>
#include <random>

#include "utils/sceneparser.h"

namespace StarGenerator {

void generateStars(RenderData& renderData);

glm::vec3 generatePosition(std::vector<glm::vec3>& positions, std::mt19937& generator);

float generateSize(std::mt19937& generator, std::normal_distribution<float>& size);

glm::vec3 generateColor(std::mt19937& generator, std::uniform_int_distribution<>& color);

glm::vec2 generateIntensity(std::mt19937& generator,
                            std::normal_distribution<float>& intensity,
                            std::normal_distribution<float>& shine);

float lighten(const float color);

void generatePrimitive(std::mt19937& generator, RenderData& renderData, std::vector<glm::vec3>& positions);

}
