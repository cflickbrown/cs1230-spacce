#include <stdexcept>
#include "raytracescene.h"
#include "utils/sceneparser.h"

RayTraceScene::RayTraceScene(int width, int height, const RenderData &metaData)
    : w(width), h(height), globalData(metaData.globalData),
    camera(Camera()), lights(metaData.lights) {

    SceneCameraData cameraInfo = metaData.cameraData;

    Camera newCamera = Camera(cameraInfo.pos, cameraInfo.look, cameraInfo.up, RayTraceScene::height(), RayTraceScene::width(),
                              cameraInfo.heightAngle, cameraInfo.focalLength, cameraInfo.aperture);

    Camera& newCameraAddr = newCamera;

    RayTraceScene::camera = newCameraAddr;

    std::vector<RenderShapeData> tempShapes(metaData.shapes.size());

    for (int i = 0; i < tempShapes.size(); i++) {
        tempShapes[i] = metaData.shapes[i];
        if(metaData.shapes[i].primitive.material.textureMap.isUsed) {
            tempShapes[i].primitive.material.texture = loadImageFromFile(metaData.shapes[i].primitive.material.textureMap.filename);
        }
    }

    RayTraceScene::shapes = tempShapes;

}

const int& RayTraceScene::width() const {
    return w;
}

const int& RayTraceScene::height() const {
    return h;
}

const SceneGlobalData& RayTraceScene::getGlobalData() const {
    return globalData;
}

const Camera& RayTraceScene::getCamera() const {
    return camera;
}

const std::vector<RenderShapeData>& RayTraceScene::getShapes() const {
    return shapes;
}

const std::vector<SceneLightData>& RayTraceScene::getLights() const {
    return lights;
}
