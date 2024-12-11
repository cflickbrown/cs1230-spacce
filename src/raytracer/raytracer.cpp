#include "raytracer.h"
#include "raytracescene.h"
#include "utils/shapes.h"

#include "utils/imagereader.h"


RayTracer::RayTracer(Config config) :
    m_config(config)
{}

float TRANSPARENCY_THRESH = 0.001;
int STEP_THRESH = 300;
glm::vec4 BACKGROUND_COLOR = glm::vec4(1,1,1,1);

IntersectionData findIntersectDataForShape(RenderShapeData shape, glm::vec4 origin, glm::vec4 direction) {
    glm::mat4 objSpaceConv = glm::inverse(shape.ctm);
    glm::vec4 rayOriginInObj = objSpaceConv * origin;
    glm::vec4 rayDirInObj = objSpaceConv * direction;

    float bestTForObj;
    ScenePrimitive currPrimitive;

    //set shape to be some specific primitive
    switch(shape.primitive.type) {
    case PrimitiveType::PRIMITIVE_CUBE:
        currPrimitive = Cube(shape.primitive.type, shape.primitive.material, shape.primitive.meshfile);
        break;
    case PrimitiveType::PRIMITIVE_SPHERE:
        currPrimitive = Sphere(shape.primitive.type, shape.primitive.material, shape.primitive.meshfile);
        break;
    case PrimitiveType::PRIMITIVE_CYLINDER:
        currPrimitive = Cylinder(shape.primitive.type, shape.primitive.material, shape.primitive.meshfile);
        break;
    case PrimitiveType::PRIMITIVE_CONE:
        currPrimitive = Cone(shape.primitive.type, shape.primitive.material, shape.primitive.meshfile);
        break;
    default:
        break;
    }

    return currPrimitive.getIntersectData(rayOriginInObj, rayDirInObj);
}

/**
 * @brief findIntersectDataForShapes gets all nearest intersections for a ray through a list of shapes
 * @param shapes a list of shapes to be intersected
 * @param origin the origin point of the ray in world space
 * @param direction the (normal) direction of the ray in world space
 * @return a list of all intersects with the ray that are closest to the origin
 */
std::vector<IntersectionData> findIntersectDataForShapes(std::vector<RenderShapeData> shapes, glm::vec4 origin, glm::vec4 direction) {
    std::vector<IntersectionData> smallestIntsForObjects(shapes.size());

    for (int i = 0; i < shapes.size(); i++) {
        RenderShapeData shape = shapes[i];

        //get the potential nearest intersect for this shape
        smallestIntsForObjects[i] = findIntersectDataForShape(shape, origin, direction);
    }

    return smallestIntsForObjects;
}


/**
 * @brief toRGBA Clamps a vec4 down to a valid RGBA value
 * @param value some vec4
 * @return an RGBA clamped to [0, 255]
 */
RGBA toRGBA(const glm::vec4 &value) {
    std::uint8_t redColor = 255 * fmin(fmax(value.r, 0), 1);
    std::uint8_t greenColor = 255 * fmin(fmax(value.g, 0), 1);
    std::uint8_t blueColor = 255 * fmin(fmax(value.b, 0), 1);

    return RGBA{redColor, greenColor, blueColor};
}


/**
 * @brief RayTracer::render Renders an image on some canvas of RGBA, given some scene
 * @param imageData the canvas of RGBA to render upon
 * @param scene The scene information, containing primitives, lights, and a camera
 */
void RayTracer::render(RGBA *imageData, const RayTraceScene &scene) {

    //extract data
    Camera useCamera = scene.getCamera();
    SceneGlobalData useGlobalData = scene.getGlobalData();

    //get viewplane
    float kDepth = 1.0f;
    float vpWidth = 2.0 * kDepth * tan(useCamera.getWidthAngle() / 2.0);
    float vpHeight = 2.0 * kDepth * tan(useCamera.getHeightAngle() / 2.0);

    //get the camera position in world space
    glm::vec4 cameraInWorld = glm::inverse(useCamera.getViewMatrix()) * glm::vec4(0,0,0,1);

    for(int ri = 0; ri < scene.height(); ri++) {
        std::cout << ri << std::endl;
        for(int ci = 0; ci < scene.width(); ci++) {
            //get the ray for the pixel (todo: allow for upsampling)
            float rayXCam = vpWidth * (((ci + 0.5) / scene.width()) - 0.5);
            float rayYCam = vpHeight * (((scene.height() - 1.0 - ri + 0.5) / scene.height()) - 0.5);
            float rayZCam = -1.0 * kDepth;

            //find where the ray intersects the viewplane in camera space
            glm::vec4 vpPointInCam = glm::vec4(rayXCam, rayYCam, rayZCam, 1.0);

            //get the direction of the ray in world space
            glm::vec4 rayDirInWorld = glm::normalize((glm::inverse(useCamera.getViewMatrix()) * vpPointInCam) - cameraInWorld);

            // imageData[ri * scene.width() + ci] = toRGBA(recurRayTrace(cameraInWorld, rayDirInWorld, scene, scene.getShapes(), scene.getShapes().size(), 0));

            // imageData[ri * scene.width() + ci] = toRGBA(recurRayMarch(cameraInWorld, rayDirInWorld, scene, 0.05, 0, 1.f, glm::vec4(0,0,0,0)));
            imageData[ri * scene.width() + ci] = toRGBA(traceMarchOrBackground(cameraInWorld, rayDirInWorld, scene, scene.getShapes(), scene.getShapes().size(), 0, 1.f));

        }
    }
}

/**
 * @brief spotlightFalloff calculates the falloff for a spotlight given some angle information
 * @param angleOffset the angle between the center of the light and the object being lit
 * @param lightAngle the cumulative, total angle of the spotlight
 * @param lightPenumbra the angle between the inner (solid) and outer (falloff) sections of the spotlight
 * @return the amount of light at the offset, in the range [0, 1]
 */
float spotlightFalloff(float angleOffset, float lightAngle, float lightPenumbra) {
    float innerAngle = lightAngle - lightPenumbra;
    float termOne = -2.0 * pow((angleOffset - innerAngle) / (lightAngle - innerAngle), 3);
    float termTwo = 3 * pow((angleOffset - innerAngle) / (lightAngle - innerAngle), 2);
    return termOne + termTwo;
}

/**
 * @brief RayTracer::PhongIllumination Applies illumination effects at an intersection point on an object primitive
 * @param position The position of the intersect in world space
 * @param normal The normal at the intersection point in world space
 * @param directionToCamera The direction from the intersection point to the camera in world space
 * @param material The material of the intersected primitive
 * @param uvCoord the uv coordinates on the texture map for the point
 * @param lights The vector of all lights in the scene
 * @param globalData The globalData of the scene, containing the lighting factors
 * @param relevantShapes all shapes that should be checked for potential shading
 * @param droppedShapeIdx the index of the current intersected shape (removed from relevantShapes) in the list of all shapes
 * @param scene the total scene data
 * @param reclevel the recursion level of the raytrace
 * @return A RGBA pixel to be displayed after lighting has been applied
 */
glm::vec4 RayTracer::PhongIllumination(glm::vec3  position,
           glm::vec3  normal,
           glm::vec3  directionToOrigin,
           SceneMaterial  &material,
           glm::vec2 uvCoord,
           std::vector<SceneLightData> &lights,
           SceneGlobalData globalData,
           std::vector<RenderShapeData> relevantShapes,
           int droppedShapeIdx,
           const RayTraceScene &scene,
           int recLevel) {
    // Normalizing directions
    normal            = glm::normalize(normal);
    directionToOrigin = glm::normalize(directionToOrigin);

    // Output illumination (we can ignore opacity)
    glm::vec4 illumination(0, 0, 0, 1);

    //add the ambient term
    illumination.r = globalData.ka * material.cAmbient.r;
    illumination.g = globalData.ka * material.cAmbient.g;
    illumination.b = globalData.ka * material.cAmbient.b;

    for (const SceneLightData &light : lights) {
        float fAtt = 1;
        glm::vec3 pointToLight = light.dir;
        float distToLight = FLT_MAX;
        float lightIntensity = 1.0;

        //set distance, intensity, attuation, position of light according to type
        switch(light.type) {
        case LightType::LIGHT_DIRECTIONAL:
            pointToLight = glm::normalize(pointToLight);
            pointToLight *= -1;
            break;
        case LightType::LIGHT_POINT:
            distToLight = glm::length((glm::vec3)light.pos - position);
            pointToLight = glm::normalize((glm::vec3)light.pos - position);
            fAtt = fmin(1, 1/(light.function[0] + (distToLight * light.function[1]) + (distToLight * distToLight * light.function[2])));
            break;
        case LightType::LIGHT_SPOT:
            distToLight = glm::length((glm::vec3)light.pos - position);
            pointToLight = glm::normalize((glm::vec3)light.pos - position);
            glm::vec3 lightToPoint = glm::normalize(position - (glm::vec3)light.pos);

            //calculate the amount of light on the point
            float angleOffset = acos(dot(glm::normalize(lightToPoint), glm::normalize((glm::vec3)light.dir)));

            if(angleOffset > light.angle) {
                lightIntensity = 0;
            } else if (angleOffset < light.angle - light.penumbra) {
                lightIntensity = 1;
            } else {
                lightIntensity = 1 - spotlightFalloff(angleOffset, light.angle, light.penumbra);
            }

            fAtt = fmin(1, 1/(light.function[0] + distToLight * light.function[1] + distToLight * distToLight * light.function[2]));
            break;
        }

        //find all possible intersections between the lit point and the light. if one exists, the light does not reach the point.
        std::vector<IntersectionData> allPossibleIntersects = findIntersectDataForShapes(relevantShapes, glm::vec4(position, 1.f), glm::vec4(pointToLight, 0.f));

        bool shaded = false;

        for(IntersectionData intersect : allPossibleIntersects) {
            //catch for if the light is closer to the origin than the ray intersection, and if the intersect is behind the point
            float distToIntersect = glm::distance(position + intersect.intersectT * pointToLight, position);
            if(intersect.hasIntersect && intersect.intersectT > 0.0 && distToIntersect < distToLight && intersect.solid) {
                shaded = true;
            }
        }

        if(!shaded || !m_config.enableShadow) {
            //calculate the diffuse values
            float diffuseRed = globalData.kd * material.cDiffuse.r;
            float diffuseGreen = globalData.kd * material.cDiffuse.g;
            float diffuseBlue = globalData.kd * material.cDiffuse.b;

            //calculate the texture
            if(material.textureMap.isUsed && m_config.enableTextureMap) {
                //get the intialized texture and blend with the diffuse value
                Image* intersectTexture = material.texture;
                int colIdx = (int)floor(fmin(fmax(0, uvCoord[0]), 1) * material.textureMap.repeatU * intersectTexture->width) % intersectTexture->width;
                int rowIdx = (int)floor((1 - fmin(1, uvCoord[1])) * material.textureMap.repeatV * intersectTexture->height) % intersectTexture->height;

                if(colIdx == material.textureMap.repeatU * intersectTexture->width) {
                    colIdx -= 1;
                }
                if(rowIdx == material.textureMap.repeatV * intersectTexture->height) {
                    rowIdx -= 1;
                }

                RGBA texturePixel = material.texture->data[(rowIdx * intersectTexture->width) + colIdx];

                diffuseRed = (material.blend * (texturePixel.r/255.0)) + ((1.0 - material.blend) * diffuseRed);
                diffuseGreen = (material.blend * (texturePixel.g/255.0)) + ((1.0 - material.blend) * diffuseGreen);
                diffuseBlue = (material.blend * (texturePixel.b/255.0)) + ((1.0 - material.blend) * diffuseBlue);
            }

            //finalize diffuse lighting on the calculated surface
            diffuseRed *= fAtt * light.color.r * lightIntensity * dot(normal, pointToLight);
            diffuseGreen *= fAtt * light.color.g * lightIntensity * dot(normal, pointToLight);
            diffuseBlue *= fAtt * light.color.b * lightIntensity * dot(normal, pointToLight);

            //apply diffuse lighting
            if(dot(normal, pointToLight) > 0) {
                illumination.r += fmax(diffuseRed, 0);
                illumination.g += fmax(diffuseGreen, 0);
                illumination.b += fmax(diffuseBlue, 0);
            }

            //calculate specular values
            glm::vec3 reflectionLightDirection = (2*dot(glm::normalize(pointToLight), normal) * normal) - glm::normalize(pointToLight);
            float specLightDistribution = pow(dot(reflectionLightDirection, directionToOrigin), material.shininess);

            float specRed = fAtt * light.color.r * lightIntensity * globalData.ks * material.cSpecular.r * specLightDistribution;
            float specGreen = fAtt * light.color.g * lightIntensity * globalData.ks * material.cSpecular.g * specLightDistribution;
            float specBlue = fAtt * light.color.b * lightIntensity * globalData.ks * material.cSpecular.b * specLightDistribution;

            if(dot(normal, pointToLight) > 0) {
                illumination.r += fmax(specRed, 0);
                illumination.g += fmax(specGreen, 0);
                illumination.b += fmax(specBlue, 0);
            }
        }
    }

    // //if there's a reflection, recur the raytrace and find the reflected value, adding it to the illumination
    // if(m_config.enableReflection && !glm::all(glm::equal(material.cReflective, glm::vec4(0,0,0,0)))) {
    //     glm::vec3 directionFromOrigin = directionToOrigin;
    //     directionFromOrigin *= -1;
    //     glm::vec3 reflectRayDir = glm::normalize(directionFromOrigin)  - 2.f*glm::dot(directionFromOrigin, normal)*normal;

    //     glm::vec4 refIllum = recurRayTrace(glm::vec4(position, 1.f), glm::normalize(glm::vec4(reflectRayDir, 0.f)), scene, relevantShapes, droppedShapeIdx, recLevel + 1);

    //     illumination.r += globalData.ks * refIllum.r * material.cReflective.r;
    //     illumination.g += globalData.ks * refIllum.g * material.cReflective.g;
    //     illumination.b += globalData.ks * refIllum.b * material.cReflective.b;
    // }

    return illumination;
}



glm::vec4 RayTracer::traceMarchOrBackground(glm::vec4 rayOrigin, glm::vec4 rayDirection, const RayTraceScene &scene, std::vector<RenderShapeData> relevantShapes, int droppedShapeIdx, int recLevel, float transparency) {

    std::vector<IntersectionData> smallestIntsForObjects(relevantShapes.size());

    //iterate through all shapes
    smallestIntsForObjects = findIntersectDataForShapes(relevantShapes, rayOrigin, rayDirection);

    //select the nearest (smallest) intersect across all shapes for the ray
    int nearestIntShapeIdx = -1;
    IntersectionData nearestIntersect = IntersectionData(false, FLT_MAX, glm::vec4(0));
    for(int i = 0; i < smallestIntsForObjects.size(); i++) {
        if(smallestIntsForObjects[i].hasIntersect && smallestIntsForObjects[i].intersectT < nearestIntersect.intersectT && smallestIntsForObjects[i].intersectT > 0) {
            nearestIntersect = smallestIntsForObjects[i];
            nearestIntShapeIdx = i;
        }
    }

    if(nearestIntShapeIdx < 0) {
        return BACKGROUND_COLOR;
    } else if (relevantShapes[nearestIntShapeIdx].primitive.material.solid) {
        return recurRayTrace(rayOrigin, rayDirection, scene, relevantShapes, droppedShapeIdx, recLevel);
    } else {
        return recurRayMarch(rayOrigin, rayDirection, scene, 0.05, 0, 1.f, glm::vec4(0,0,0,0), 0);
    }

}



/**
 * @brief RayTracer::recurRayTrace a recursive raytacer, giving the illumination for some ray in a scene with a relevant set of shapes
 * @param rayOrigin the origin of the ray to be traced
 * @param rayDirection the (normalized) direction of the ray to be traced
 * @param scene the total scene the ray exists in
 * @param relevantShapes the set of shapes that may be intersected for a given recursion
 * @param droppedShapeIdx the index of a removed shape from the list of shapes. if no shapes are removed, then this is the size of the list of all shapes
 * @param recLevel the recursion level of the raytrace
 * @return a vector representing the illumination values seen by the ray at some origin
 */
glm::vec4 RayTracer::recurRayTrace(glm::vec4 rayOrigin, glm::vec4 rayDirection, const RayTraceScene &scene, std::vector<RenderShapeData> relevantShapes, int droppedShapeIdx, int recLevel) {
    //4 is a typical heuristic, but could be something else here
    if(recLevel < m_config.maxRecursiveDepth) {

        std::vector<IntersectionData> smallestIntsForObjects(relevantShapes.size());

        //iterate through all shapes
        smallestIntsForObjects = findIntersectDataForShapes(relevantShapes, rayOrigin, rayDirection);

        //select the nearest (smallest) intersect across all shapes for the ray
        int nearestIntShapeIdx = 0;
        IntersectionData nearestIntersect = IntersectionData(false, FLT_MAX, glm::vec4(0));
        for(int i = 0; i < smallestIntsForObjects.size(); i++) {
            if(smallestIntsForObjects[i].hasIntersect && smallestIntsForObjects[i].intersectT < nearestIntersect.intersectT && smallestIntsForObjects[i].intersectT > 0) {
                nearestIntersect = smallestIntsForObjects[i];
                nearestIntShapeIdx = i;
            }
        }

        if(!nearestIntersect.hasIntersect) {
            return glm::vec4(0.f,0.f,0.f,0.f);
        } else {
            RenderShapeData intersectedShape = relevantShapes[nearestIntShapeIdx];

            //get normal at intersect point, the intersect point, and the direction from the intersect point to the origin in world space
            glm::vec3 normalInWorld = glm::inverse(glm::transpose(intersectedShape.ctm)) * nearestIntersect.normAtIntForObj;
            glm::vec3 intPointInWorld = rayOrigin + nearestIntersect.intersectT * rayDirection;
            glm::vec3 intToOriginInWorld = glm::normalize(glm::vec3(rayOrigin) - intPointInWorld);

            std::vector<SceneLightData> useLights = scene.getLights();

            //get the set of relevant shapes, which are all shapes except the one currently intersected by the ray
            std::vector<RenderShapeData> newRelevantShapes = scene.getShapes();

            if(nearestIntShapeIdx >= droppedShapeIdx) {
                nearestIntShapeIdx += 1;
            }

            newRelevantShapes.erase(newRelevantShapes.begin() + nearestIntShapeIdx);

            if(m_config.onlyRenderNormals) {
                int normaledRed = fmin(255, 255 * (normalInWorld.r + 1) / 2);
                int normaledGreen = fmin(255, 255 * (normalInWorld.g + 1) / 2);
                int normaledBlue = fmin(255, 255 * (normalInWorld.b + 1) / 2);
                return glm::vec4(normaledRed/255.0, normaledGreen/255.0, normaledBlue/255.0, 0.f);
            }

            glm::vec4 illuminatedPixel = RayTracer::PhongIllumination(intPointInWorld,
                                                                      normalInWorld,
                                                                      intToOriginInWorld,
                                                                      intersectedShape.primitive.material,
                                                                      nearestIntersect.uvCoords,
                                                                      useLights,
                                                                      scene.getGlobalData(),
                                                                      newRelevantShapes,
                                                                      nearestIntShapeIdx,
                                                                      scene,
                                                                      recLevel);
            return illuminatedPixel;
        }
    } else {
        return glm::vec4(0.f,0.f,0.f,0.f);
    }
}


/**
 * @brief findIntersectDataForShapes gets all nearest intersections for a ray through a list of shapes
 * @param shapes a list of shapes to be intersected
 * @param origin the origin point of the ray in world space
 * @param direction the (normal) direction of the ray in world space
 * @return a list of all intersects with the ray that are closest to the origin
 */
std::vector<float> findDensityDataForShapes(std::vector<RenderShapeData> shapes, glm::vec4 point) {
    std::vector<float> objectDensities(shapes.size());

    for (int i = 0; i < shapes.size(); i++) {
        RenderShapeData shape = shapes[i];

        //get origin, ray direcion in object space
        glm::mat4 objSpaceConv = glm::inverse(shape.ctm);
        glm::vec4 pointInObj = objSpaceConv * point;

        float bestTForObj;
        ScenePrimitive currPrimitive;

        //set shape to be some specific primitive
        switch(shape.primitive.type) {
        case PrimitiveType::PRIMITIVE_CUBE:
            currPrimitive = Cube(shape.primitive.type, shape.primitive.material, shape.primitive.meshfile);
            break;
        case PrimitiveType::PRIMITIVE_SPHERE:
            currPrimitive = Sphere(shape.primitive.type, shape.primitive.material, shape.primitive.meshfile);
            break;
        case PrimitiveType::PRIMITIVE_CYLINDER:
            currPrimitive = Cylinder(shape.primitive.type, shape.primitive.material, shape.primitive.meshfile);
            break;
        case PrimitiveType::PRIMITIVE_CONE:
            currPrimitive = Cone(shape.primitive.type, shape.primitive.material, shape.primitive.meshfile);
            break;
        default:
            break;
        }

        //get the potential nearest intersect for this shape
        objectDensities[i] = currPrimitive.getDensityData(pointInObj);
    }

    return objectDensities;
}


/**
 * @brief RayTracer::recurRayMarch a recursive raymarcher, giving the value [0, 1] of some pixel at some point
 * @param rayOrigin the origin of the ray to be marched
 * @param rayDirection the (normalized) direction of the ray to be marched
 * @param scene the total scene the ray exists in
 * @param stepSize the length of the ray being marched along the direction
 * @param numOfSteps the number of marched steps completed
 * @param transparency the amount of transparency at some the marched point
 * @param pixelResult the running result of all sampled points
 * @return a vector representing the illumination values seen by the ray at some origin
 */

glm::vec4 RayTracer::recurRayMarch(glm::vec4 rayOrigin, glm::vec4 rayDirection, const RayTraceScene &scene, float stepSize, int numOfSteps, float transparency, glm::vec4 pixelResult, float accContactDist) {
    if(transparency < TRANSPARENCY_THRESH) {
        if(accContactDist > 0) {
            return glm::vec4(glm::vec3(pixelResult / accContactDist), 1);
        } else {
            return glm::vec4(glm::vec3(pixelResult), 1);
        }
    } else if (numOfSteps > STEP_THRESH) {
        if(pixelResult.r > 1 || pixelResult.g > 1 || pixelResult.b > 1) {
            return glm::vec4(1,0,0,1);
        }
        if(accContactDist > 0) {
            return glm::vec4(glm::vec3((BACKGROUND_COLOR * transparency) + (pixelResult / accContactDist)), 1.f);
        } else {
            return glm::vec4(glm::vec3((BACKGROUND_COLOR * transparency) + (pixelResult)), 1.f);
        }
    } else {
        std::vector<float> objectDensities(scene.getShapes().size());

        glm::vec4 rayEndpoint = rayOrigin + (stepSize * rayDirection);

        //In reality this should be the middle of the ray step, not the end, but also who cares for now

        glm::vec4 raySelectPoint = rayOrigin + (((float)rand() / ((float)RAND_MAX)) * stepSize * rayDirection);

        objectDensities = findDensityDataForShapes(scene.getShapes(), raySelectPoint);

        //select the highest density across all shapes
        int highestDensityIdx = 0;
        float highestDensity = 0;
        for(int i = 0; i < objectDensities.size(); i++) {
            if(objectDensities[i] > highestDensity) {
                highestDensity = objectDensities[i];
                highestDensityIdx = i;
            }
        }

        //do lighting stuff for that object

        float currentTransparency = transparency;
        if (highestDensity > 100) {
            //object is basically solid, treat it as such

            std::vector<RenderShapeData> relevantShapes = scene.getShapes();
            glm::vec4 cameraInWorld = glm::inverse(scene.getCamera().getViewMatrix()) * glm::vec4(0,0,0,1);

            IntersectionData intForObj = findIntersectDataForShape(scene.getShapes()[highestDensityIdx], cameraInWorld, rayDirection);

            RenderShapeData intersectedShape = relevantShapes[highestDensityIdx];

            //get normal at intersect point, the intersect point, and the direction from the intersect point to the origin in world space
            glm::vec3 normalInWorld = glm::inverse(glm::transpose(intersectedShape.ctm)) * intForObj.normAtIntForObj;
            glm::vec3 intPointInWorld = cameraInWorld + intForObj.intersectT * rayDirection;
            glm::vec3 intToOriginInWorld = glm::normalize(glm::vec3(cameraInWorld) - intPointInWorld);

            std::vector<SceneLightData> useLights = scene.getLights();

            //get the set of relevant shapes, which are all shapes except the one currently intersected by the ray
            std::vector<RenderShapeData> newRelevantShapes = scene.getShapes();

            newRelevantShapes.erase(newRelevantShapes.begin() + highestDensityIdx);

            if(m_config.onlyRenderNormals) {
                int normaledRed = fmin(255, 255 * (normalInWorld.r + 1) / 2);
                int normaledGreen = fmin(255, 255 * (normalInWorld.g + 1) / 2);
                int normaledBlue = fmin(255, 255 * (normalInWorld.b + 1) / 2);
                return glm::vec4(normaledRed/255.0, normaledGreen/255.0, normaledBlue/255.0, 0.f);
            }

            glm::vec4 illuminatedPixel = RayTracer::PhongIllumination(intPointInWorld,
                                                                      normalInWorld,
                                                                      intToOriginInWorld,
                                                                      intersectedShape.primitive.material,
                                                                      intForObj.uvCoords,
                                                                      useLights,
                                                                      scene.getGlobalData(),
                                                                      newRelevantShapes,
                                                                      highestDensityIdx,
                                                                      scene,
                                                                      0);

            //todo: blend!
            if(accContactDist > 0) {
                pixelResult = (1 - transparency) * (pixelResult / accContactDist);
            }
            pixelResult += transparency * illuminatedPixel;
            // return pixelResult;
            return recurRayMarch(rayEndpoint, rayDirection, scene, stepSize, numOfSteps + 1, 0, pixelResult, 0);

        } else if(highestDensity > 0) {
            //do lighting...

            float attentuationAtPoint = exp(-stepSize * highestDensity);
            currentTransparency *= attentuationAtPoint;

            glm::vec4 accLightEffect = glm::vec4(0,0,0,1);

            for(SceneLightData light : scene.getLights()) {
                glm::vec4 pointToLight = light.dir;
                float fAtt = 1.0;

                if(light.type == LightType::LIGHT_POINT || light.type == LightType::LIGHT_SPOT) {
                    pointToLight = glm::vec4(glm::normalize(glm::vec3(light.pos - raySelectPoint)), 0);
                    float distToLight = glm::length(glm::vec3(light.pos - raySelectPoint));
                    fAtt = fmin(1, 1/(light.function[0] + distToLight * light.function[1] + distToLight * distToLight * light.function[2]));
                } else {
                    pointToLight = glm::normalize(pointToLight);
                    pointToLight *= -1;
                }


                IntersectionData intForObj = findIntersectDataForShape(scene.getShapes()[highestDensityIdx], raySelectPoint, pointToLight);
                accLightEffect += (float)exp(-intForObj.intersectT * highestDensity * stepSize) * light.color * scene.getGlobalData().kd * stepSize;
            }

            //this is REALLY basic for now...
            //
            // pixelResult += currentTransparency * accLightEffect * scene.getShapes().at(highestDensityIdx).primitive.material.cDiffuse * stepSize;
            pixelResult += accLightEffect * transparency * scene.getShapes().at(highestDensityIdx).primitive.material.cDiffuse;
            // pixelResult += scene.getGlobalData().ka * scene.getShapes().at(highestDensityIdx).primitive.material.cAmbient * stepSize;
            return recurRayMarch(rayEndpoint, rayDirection, scene, stepSize, numOfSteps + 1, currentTransparency, pixelResult, accContactDist + stepSize);
        }



        return recurRayMarch(rayEndpoint, rayDirection, scene, stepSize, numOfSteps + 1, currentTransparency, pixelResult, accContactDist);
    }
}

//TODO:
//- better lighting effects
    //shadows for
//- speed-ups (skip pixels with no intersects)

