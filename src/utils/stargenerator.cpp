#include "utils/stargenerator.h"

#include <random>
#include <glm/gtx/transform.hpp>

namespace StarGenerator {

void generateStars(RenderData& renderData) {

    std::vector<glm::vec3> positions;
    generatePosition(positions);
    generatePrimitive(renderData, positions);
}

glm::vec3 generatePosition(std::vector<glm::vec3>& positions) {

    std::random_device rd;
    std::mt19937 generator(rd());

    // std::mt19937 generator();
    // generator.seed(1);

    std::uniform_real_distribution<> distance(0.0f, 1.0f);

    int gridDivisions = 5;

    float rangeMax = 5.0f;
    float stepSize = (2 * rangeMax) / gridDivisions;

    for (int i = 0; i < gridDivisions; i++) {
        for (int j = 0; j < gridDivisions; j++) {
            for (int k = 0; k < gridDivisions; k++) {

                float xMin = -rangeMax + i * stepSize;
                float xMax = xMin + stepSize;

                float yMin = -rangeMax + j * stepSize;
                float yMax = yMin + stepSize;

                float zMin = -rangeMax + k * stepSize;
                float zMax = zMin + stepSize;

                float x = xMin + distance(generator) * (xMax - xMin);
                float y = yMin + distance(generator) * (yMax - yMin);
                float z = zMin + distance(generator) * (zMax - zMin);

                positions.push_back(glm::vec3(x, y, z));
            }
        }
    }
}

void generatePrimitive(RenderData& renderData, std::vector<glm::vec3>& positions) {

    for (int i = 0; i < positions.size(); i++) {
        glm::mat4 ctm = glm::mat4(1.0f);

        ctm = translate(ctm, positions[i]);

        ctm = scale(ctm, glm::vec3(0.5, 0.5, 0.5));

        ScenePrimitive* star = new ScenePrimitive();
        SceneMaterial &material = star->material;
        material.clear();
        star->type = PrimitiveType::PRIMITIVE_SPHERE;

        for (int i = 0; i < 3; i++) {
            material.cAmbient[i] = 0.2;
        }

        for (int i = 0; i < 3; i++) {
            material.cDiffuse[i] = 0.8;
        }

        for (int i = 0; i < 3; i++) {
            material.cSpecular[i] = 0.5;
        }

        material.shininess = 50.0f;

        renderData.shapes.push_back(RenderShapeData(*star, ctm));
    }
}

}
