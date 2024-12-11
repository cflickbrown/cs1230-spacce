#pragma once

#include <glm/glm.hpp>

#include "utils/sceneparser.h"

namespace StarGenerator {

void generateStars(RenderData& renderData);

glm::vec3 generatePosition(std::vector<glm::vec3>& positions);

void generatePrimitive(RenderData& renderData, std::vector<glm::vec3>& positions);

}
