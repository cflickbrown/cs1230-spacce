#include <stdexcept>
#include "camera.h"


Camera::Camera() {
    //This space intentionally left blank.
}

Camera::Camera(glm::vec3 pos, glm::vec3 look, glm::vec3 up, int viewHeight, int viewWidth, float heightAngle, float focalLength, float aperture) {
    glm::vec4 viewVec = glm::vec4(-pos, 1.0f);
    glm::mat4 mTransl = glm::mat4(1.0f);
    mTransl[3] = viewVec;

    glm::vec3 vecW = look;
    vecW = -normalize(vecW);
    glm::vec3 vecV = up - dot(up, vecW)*vecW;
    vecV = normalize(vecV);
    glm::vec3 vecU = cross(vecV, vecW);

    glm::vec4 mRotCol1 = glm::vec4(vecU[0], vecV[0], vecW[0], 0);
    glm::vec4 mRotCol2 = glm::vec4(vecU[1], vecV[1], vecW[1], 0);
    glm::vec4 mRotCol3 = glm::vec4(vecU[2], vecV[2], vecW[2], 0);
    glm::vec4 mRotCol4 = glm::vec4(0, 0, 0, 1);

    glm::mat4 mRotat = glm::mat4(mRotCol1,mRotCol2,mRotCol3,mRotCol4);

    glm::mat4 viewMatrix = mRotat * mTransl;

    Camera::viewMatrix = viewMatrix;

    Camera::aspectRatio = (1.0 * viewWidth)/(1.0 * viewHeight);
    Camera::heightAngle = heightAngle;
    Camera::focalLength = focalLength;
    Camera::aperture = aperture;
}

glm::mat4 Camera::getViewMatrix() const {
    return viewMatrix;
}

float Camera::getAspectRatio() const {
    return aspectRatio;

}

float Camera::getHeightAngle() const {
    return heightAngle;
}

float Camera::getWidthAngle() const {
    return 2.0 * atan(getAspectRatio() * tan(getHeightAngle() / 2));
}

float Camera::getFocalLength() const {
    return focalLength;
}

float Camera::getAperture() const {
    return aperture;
}
