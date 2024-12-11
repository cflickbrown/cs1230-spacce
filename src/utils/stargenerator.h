#pragma once

#include <glm/glm.hpp>

#include "utils/sceneparser.h"

namespace StarGenerator {

void generateStars(RenderData& renderData);

glm::vec3 generatePosition();

void generatePrimitive(RenderData& renderData, glm::vec3& pos);

}
