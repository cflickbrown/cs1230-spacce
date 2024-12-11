#include "utils/stargenerator.h"

#include <random>
#include <glm/gtx/transform.hpp>

namespace StarGenerator {

// okay so you need a function to generate position first
// then a function that generates the matrices and that will construct ScenePrimitives and put them in the data

void generateStars(RenderData& renderData) {

    // calls the position generation
    // then generates matrices and pushes back

    // glm::vec3 pos = generatePosition();
    glm::vec3 pos(0, 0, 0);
    generatePrimitive(renderData, pos);
}

glm::vec3 generatePosition() {

    std::mt19937 generator;
    generator.seed(1);

    std::uniform_real_distribution<> distance(0.0f, 1.0f);

    double x = distance(generator);
    double y = distance(generator);
    double z = distance(generator);

    return glm::vec3(x, y, z);
}

void generatePrimitive(RenderData& renderData, glm::vec3& pos) {

    // generate the matrix and make the primitive and add it to the data

    glm::mat4 ctm = glm::mat4(1.0f);

    ctm = translate(ctm, pos);

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
