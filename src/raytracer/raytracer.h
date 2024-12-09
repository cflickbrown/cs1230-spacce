#pragma once

#include <glm/glm.hpp>
#include "utils/rgba.h"
#include "utils/scenedata.h"
#include "utils/sceneparser.h"

// A forward declaration for the RaytraceScene class

class RayTraceScene;

// A class representing a ray-tracer

class RayTracer
{
public:
    struct Config {
        bool enableShadow        = false;
        bool enableReflection    = false;
        bool enableRefraction    = false;
        bool enableTextureMap    = false;
        bool enableTextureFilter = false;
        bool enableParallelism   = false;
        bool enableSuperSample   = false;
        bool enableAcceleration  = false;
        bool enableDepthOfField  = false;
        int maxRecursiveDepth    = 4;
        bool onlyRenderNormals   = false;
    };

public:
    RayTracer(Config config);

    // Renders the scene synchronously.
    // The ray-tracer will render the scene and fill imageData in-place.
    // @param imageData The pointer to the imageData to be filled.
    // @param scene The scene to be rendered.
    void render(RGBA *imageData, const RayTraceScene &scene);

    glm::vec4 recurRayTrace(glm::vec4 rayOrigin, glm::vec4 rayDirection, const RayTraceScene &scene, std::vector<RenderShapeData> relevantShapes, int droppedShapeIdx, int recLevel);

    glm::vec4 PhongIllumination(glm::vec3  position,
                           glm::vec3  normal,
                           glm::vec3  directionToCamera,
                           SceneMaterial  &material,
                           glm::vec2 uvCoords,
                           std::vector<SceneLightData> &lights,
                           SceneGlobalData globalData,
                           std::vector<RenderShapeData> shapes,
                           int droppedShapeIdx,
                           const RayTraceScene &scene,
                           int recLevel);

    glm::vec4 recurRayMarch(glm::vec4 rayOrigin, glm::vec4 rayDirection, const RayTraceScene &scene, float stepSize, int numOfSteps, float transparency, glm::vec4 pixelResult);

    glm::vec4 traceMarchOrBackground(glm::vec4 rayOrigin, glm::vec4 rayDirection, const RayTraceScene &scene, std::vector<RenderShapeData> relevantShapes, int droppedShapeIdx, int recLevel, float transparency);


private:
    const Config m_config;
};

