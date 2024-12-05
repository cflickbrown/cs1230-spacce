#include "sceneparser.h"
#include "scenefilereader.h"
#include <glm/gtx/transform.hpp>

#include <chrono>
#include <iostream>


void traverseSceneGraph(SceneNode currentNode, glm::mat4 ctm, RenderData &renderData) {
    glm::mat4 newCtm = ctm;

    std::vector<SceneTransformation*> allTransforms = currentNode.transformations;
    std::vector<SceneTransformation> translTransforms;
    std::vector<SceneTransformation> rotTransforms;
    std::vector<SceneTransformation> scaleTransforms;

    for (SceneTransformation* transf : allTransforms) {
        switch (transf->type) {
        case TransformationType::TRANSFORMATION_TRANSLATE:
            newCtm = translate(newCtm, transf->translate);
            break;
        case TransformationType::TRANSFORMATION_ROTATE:
            newCtm = rotate(newCtm, transf->angle, transf->rotate);
            break;
        case TransformationType::TRANSFORMATION_SCALE:
            newCtm = scale(newCtm, transf->scale);
            break;
        default:
            newCtm = newCtm * transf->matrix;
            break;
        }
    }

    for (ScenePrimitive* prim : currentNode.primitives) {
        renderData.shapes.push_back(RenderShapeData(*prim, newCtm));
    }


    for (SceneLight* light : currentNode.lights) {
        SceneLightData sceneLD;

        sceneLD.id = light->id;
        sceneLD.type = light->type;
        sceneLD.color = light->color;
        sceneLD.function = light->function;

        switch(light->type) {
        case LightType::LIGHT_POINT:
            sceneLD.pos = glm::vec4(0.f, 0.f, 0.f, 1.f) * glm::transpose(newCtm);
            break;
        case LightType::LIGHT_DIRECTIONAL:
            sceneLD.dir = light->dir * newCtm;
            break;
        case LightType::LIGHT_SPOT:
            sceneLD.pos = glm::vec4(0.f, 0.f, 0.f, 1.f) * glm::transpose(newCtm);
            sceneLD.dir = light->dir * newCtm;
            sceneLD.penumbra = light->penumbra;
            sceneLD.angle = light->angle;
            break;
        default:
            break;
        }

        renderData.lights.push_back(sceneLD);
    }

    for (SceneNode* child : currentNode.children) {
        traverseSceneGraph(*child, newCtm, renderData);
    }
}

bool SceneParser::parse(std::string filepath, RenderData &renderData) {
    ScenefileReader fileReader = ScenefileReader(filepath);
    bool success = fileReader.readJSON();
    if (!success) {
        return false;
    }

    // TODO: Use your Lab 5 code here

    renderData.globalData = fileReader.getGlobalData();
    renderData.cameraData = fileReader.getCameraData();

    SceneNode* rootNode = fileReader.getRootNode();
    SceneNode actRoot = *rootNode;

    renderData.shapes.clear();

    traverseSceneGraph(actRoot, glm::mat4(1.f), renderData);

    return true;
}
