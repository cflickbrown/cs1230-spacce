#include "utils/stargenerator.h"

#include <glm/gtx/transform.hpp>

namespace StarGenerator {

void generateStars(RenderData& renderData) {

    std::random_device rd;
    std::mt19937 generator(rd());

    // std::mt19937 generator;
    // generator.seed(1230);

    PerlinNoise noise = PerlinNoise();

    std::vector<glm::vec3> positions;
    generatePosition(positions, generator, noise);
    generatePrimitive(generator, noise, renderData, positions);
}

void generatePosition(std::vector<glm::vec3>& positions, std::mt19937& generator, PerlinNoise noise) {

    std::uniform_real_distribution<> distance(0.0f, 1.0f);

    int gridDivisions = 9;

    float rangeMax = 3.0f;
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

                float probability = noise.computePerlin3d(x, y, z);

                if (distance(generator) < probability) {
                    positions.push_back(glm::vec3(x, y, z));
                }
            }
        }
    }
}

float generateSize(std::mt19937& generator, std::normal_distribution<float>& size) {
    return size(generator);
}

glm::vec3 generateColor(std::mt19937& generator, std::uniform_int_distribution<>& color) {
    int choice = color(generator);

    switch(choice) {
    case 0:
        return glm::vec3(1.0f, 1.0f, 1.0f); // white
    case 1:
        return glm::vec3(0.5f, 0.5f, 1.0f); // blue
    case 2:
        return glm::vec3(0.6f, 0.3f, 0.8f); // purple
    case 3:
        return glm::vec3(1.0f, 0.5f, 0.5f); // red
    case 4:
        return glm::vec3(1.0f, 0.3f, 0.8f); // purple-y red
    default:
        return glm::vec3(1.0f, 1.0f, 1.0f);
    }
}

glm::vec2 generateIntensity(std::mt19937& generator, PerlinNoise noise, glm::vec3 position) {

    float noiseVal = noise.computePerlin3d(position.x, position.y, position.z);

    if (noiseVal < 0.1f) {
        return getIntensity(generator, 0.3, 0.1, 15, 5);
    }
    else if (noiseVal < 0.3f) {
        return getIntensity(generator, 0.5, 0.2, 50, 10);

    }
    else if (noiseVal < 1.0f) {
        return getIntensity(generator, 0.8, 0.15, 80, 15);
    }
}

glm::vec2 getIntensity(std::mt19937& generator, float intensityMean, float intensitySD, float shineMean, float shineSD) {
    std::normal_distribution<float> intensity(intensityMean, intensitySD);
    std::normal_distribution<float> shine(shineMean, shineSD);

    float intensityVal = intensity(generator);
    float shininess = shine(generator);

    return glm::vec2(intensityVal, shininess);
}

void generatePrimitive(std::mt19937& generator, PerlinNoise noise, RenderData& renderData, std::vector<glm::vec3>& positions) {
    std::normal_distribution<float> size(0.04, 0.02);

    std::uniform_int_distribution<> colorChoice(0, 4);

    for (int i = 0; i < positions.size(); i++) {
        glm::mat4 ctm = glm::mat4(1.0f);

        ctm = rotate(ctm, glm::radians(-45.0f), glm::vec3(1.0f, 0.0f, 0.0f));
        ctm = rotate(ctm, glm::radians(45.0f), glm::vec3(0.0f, 0.0f, 1.0f));

        ctm = translate(ctm, positions[i]);

        float starSize = generateSize(generator, size);
        glm::mat4 starCTM = scale(ctm, glm::vec3(starSize));
        glm::mat4 glowCTM = scale(ctm, glm::vec3(starSize * 1.3));

        ScenePrimitive star;
        SceneMaterial &starMat = star.material;
        starMat.clear();
        star.type = PrimitiveType::PRIMITIVE_SPHERE;

        ScenePrimitive starGlow;
        SceneMaterial &glowMat = starGlow.material;
        glowMat.clear();
        starGlow.type = PrimitiveType::PRIMITIVE_CUBE;

        glm::vec3 color = generateColor(generator, colorChoice);
        glm::vec2 intensities = generateIntensity(generator, noise, positions[i]);

        for (int i = 0; i < 3; i++) {
            starMat.cAmbient[i] = color[i] * intensities[0];
            glowMat.cAmbient[i] = color[i] * intensities[0];
        }

        for (int i = 0; i < 3; i++) {
            starMat.cDiffuse[i] = color[i] * intensities[0];
            glowMat.cDiffuse[i] = color[i] * intensities[0];
        }

        for (int i = 0; i < 3; i++) {
            starMat.cSpecular[i] = 1.0f;
        }

        starMat.shininess = intensities[1];

        glowMat.solid = false;
        glowMat.density = 0.2;

        renderData.shapes.push_back(RenderShapeData(star, starCTM));
        renderData.shapes.push_back(RenderShapeData(starGlow, glowCTM));
    }
}

float lighten(float color) {
    float blend = 0.1f;
    return (color * (1.0f - blend) + 1.0f * blend);
}

}
