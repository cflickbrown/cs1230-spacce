#include "raytracer.h"
#include "raytracescene.h"
#include "utils/shapes.h"

#include "utils/imagereader.h"


RayTracer::RayTracer(Config config) :
    m_config(config)
{}

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

enum LightHitType{EVENT_HORIZON, ESCAPED, UNFINISHED};

struct relativisticOutcome{
    LightHitType outcome;
    std::vector<glm::vec4> directions;
    std::vector<glm::vec4> positions;
};

/**
 * @brief integrateGeodesic Integrates a given geodesic around a single black hole.
 * @param position Geodesic position in curved spacetime
 * @param direction Geodesic direction in curved spacetime
 * @param blackHole
 * @return A struct with an outcome and a vector or positions and directions. The raymarcher will integrate along every Nth point.
 */
relativisticOutcome integrateGeodesic(glm::vec4 position, glm::vec4 direction, BlackHole& blackHole, float stepSize){
    relativisticOutcome ret = {UNFINISHED, std::vector<glm::vec4>(), std::vector<glm::vec4>()};

    glm::vec4 posRelative = position - blackHole.position;
    //float time = position[0],r = glm::length(posRelative),theta = atan(posRelative.y/posRelative.x),phi = acos(posRelative.z/r);
    glm::mat4 schwarzchildMetric; //= blackHole.getSchwarzchildMetric(time,r,theta,phi);

    glm::vec4 currentPosition = posRelative, currentDirection = direction;

    for(int i = 0; i < 100000; ++i){
        //blackHole.getSchwarzchildMetric(currentPosition[0], currentPosition[1], currentPosition[2], currentPosition[3]);

        std::vector<glm::mat4> christoff2 = blackHole.calculateSchwarzchildChristoffel(currentPosition);

        // Now that we have all neccesary parts of the geodesic equation, we can finally integrate along the geodesic... for one single step of one single ray.
        // This is a "bit" expensive
        glm::vec4 acceleration = glm::vec4(0);
        for(int mu = 0; mu < 4; ++mu){
            float sum =0.f;
            for(int al = 0; al < 4; ++al){
                for(int be = 0; be < 4; ++be){
                    sum += -christoff2[mu][al][be] * direction[al] * direction[be];
                }
            }
            acceleration[mu] = sum;
        }

        //std::cout << "MIDDLE: " << direction[0] << " " << direction[1] << " " << direction[2] << " " << direction[3] << std::endl;

        //acceleration = glm::normalize(acceleration);
        currentDirection += acceleration * stepSize;
        currentPosition += currentDirection * stepSize;

        currentPosition[0] -= stepSize;

        ret.positions.push_back(currentPosition);
        ret.directions.push_back(currentDirection);

        if(abs(currentPosition[1]) < blackHole.rs){  // Inside EH
            std::cout << position[1] << " " << blackHole.rs << std::endl;
            ret.outcome = EVENT_HORIZON;
            return ret;
        }

        if(abs(currentPosition[1]) > 30.f){  // Escaped
            ret.outcome = ESCAPED;
            return ret;
        }
    }

    return ret;
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

    m_blackHole = BlackHole();
    m_blackHole.J = 0;
    m_blackHole.M = 0.1f;
    m_blackHole.Q = 0;
    m_blackHole.position = glm::vec4(0,0,0,1);
    m_blackHole.rs = m_blackHole.M*2;

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

            // BEGIN BLACK HOLE HELL
            float
                t = 0.f,
                r = glm::length(cameraInWorld),
                theta = cameraInWorld.x == 0 && cameraInWorld.y == 0 ? 2*3.1415926535 : atan(cameraInWorld.y/cameraInWorld.x),
                phi = cameraInWorld.z == 0 ? 2*3.1415926535 : acos(cameraInWorld.z/r);
            //theta = 3.1415926535f/2.f;
            //phi = -3.1415926535f/2.f;
            glm::mat4 schwarzchildMetric = this->m_blackHole.getSchwarzchildMetric(t,r,theta,phi);
            glm::mat4 tetradBasis = m_blackHole.getTetradBasis(t,r,theta,phi);

            glm::vec4 position = {t,r,theta,phi};
            // Weird NaN stuff sometimes happens here
            //float deltaPhi = acos((posRelative.z+rayDirection.z)/r);
            //if((posRelative.z+rayDirection.z)/r == 0){
            //    deltaPhi = 0;
            //}

            //glm::vec4 positionDelta = {0,glm::length(posRelative+rayDirection),atan((posRelative.y+rayDirection.y)/(posRelative.x+rayDirection.x)), deltaPhi};
            //glm::vec4 direction = positionDelta-position;
            //direction = glm::normalize(direction);
            //direction[0] = -1;
            // Convert the flat space coordinate to out curved spacetime by multiplying it by multiplying it into a tetrad Basis
            glm::vec4 direction = {-1,-vpPointInCam[2], vpPointInCam[1], vpPointInCam[0]};

            //std::cout << "direction: " << direction[0] << " " << direction[1] << " " << direction[2] << " " << direction[3] << std::endl;
            //std::cout << "position: " << position[0] << " " << position[1] << " " << position[2] << " " << position[3] << std::endl;
            //std::cout << "positionDelta: " << positionDelta[0] << " " << positionDelta[1] << " " << positionDelta[2] << " " << positionDelta[3] << std::endl;
            //std::cout << "posRelative: " << posRelative[0] << " " << posRelative[1] << " " << posRelative[2] << " " << posRelative[3] << std::endl;
            //std::cout << "rayDirection: " << rayDirection[0] << " " << rayDirection[1] << " " << rayDirection[2] << " " << rayDirection[3] << std::endl;

            //std::cout << direction[0] << " " << direction[1] << " " << direction[2] << " " << direction[3] << std::endl;
            direction = tetradBasis * direction;


            relativisticOutcome geodesicPath = integrateGeodesic(position, direction, m_blackHole, 0.005f);

            glm::vec4 geodesicEndDir = geodesicPath.directions.back();
            glm::vec4 geodesicEndPos = geodesicPath.positions.back();

            switch(geodesicPath.outcome){
            case EVENT_HORIZON:
                //std::cout << "HORIZON" << "\n";
                imageData[ri * scene.width() + ci] = RGBA(255,255,255);
                break;
            case ESCAPED:
                //std::cout << "ESCAPED" << "\n";
                imageData[ri * scene.width() + ci] = toRGBA(glm::normalize(glm::vec4(geodesicEndDir[2], geodesicEndDir[3],0,0)));
                break;
            case UNFINISHED:
                std::cout << "UNFINISHED" << "\n";
                imageData[ri * scene.width() + ci] = RGBA(255,0,0);

                break;
            }

            //std::cout << "geodesicEndPos: " << geodesicEndPos[0] << " " << geodesicEndPos[1] << " " << geodesicEndPos[2] << " " << geodesicEndPos[3] << "\n";
            //std::cout << "geodesicEndDir: " << geodesicEndDir[0] << " " << geodesicEndDir[1] << " " << geodesicEndDir[2] << " " << geodesicEndDir[3] << "\n";

            // END BLACK HOLE HELL

            // imageData[ri * scene.width() + ci] = toRGBA(recurRayTrace(cameraInWorld, rayDirInWorld, scene, scene.getShapes(), scene.getShapes().size(), 0));

            // real
            //imageData[ri * scene.width() + ci] = toRGBA(recurRayMarch(cameraInWorld, rayDirInWorld, scene, 0.02, 0, 1.f, glm::vec4(0,0,0,0),0));
        }
        std::cout << std::endl;
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
            if(intersect.hasIntersect && intersect.intersectT > 0.0 && distToIntersect < distToLight) {
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

    //if there's a reflection, recur the raytrace and find the reflected value, adding it to the illumination
    if(m_config.enableReflection && !glm::all(glm::equal(material.cReflective, glm::vec4(0,0,0,0)))) {
        glm::vec3 directionFromOrigin = directionToOrigin;
        directionFromOrigin *= -1;
        glm::vec3 reflectRayDir = glm::normalize(directionFromOrigin)  - 2.f*glm::dot(directionFromOrigin, normal)*normal;

        glm::vec4 refIllum = recurRayTrace(glm::vec4(position, 1.f), glm::normalize(glm::vec4(reflectRayDir, 0.f)), scene, relevantShapes, droppedShapeIdx, recLevel + 1);

        illumination.r += globalData.ks * refIllum.r * material.cReflective.r;
        illumination.g += globalData.ks * refIllum.g * material.cReflective.g;
        illumination.b += globalData.ks * refIllum.b * material.cReflective.b;
    }

    return illumination;
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

float TRANSPARENCY_THRESH = 0.001;
int STEP_THRESH = 1000;
glm::vec4 BACKGROUND_COLOR = glm::vec4(0.5,0.5,0.5,1);

glm::vec4 RayTracer::recurRayMarch(glm::vec4 rayOrigin, glm::vec4 rayDirection, const RayTraceScene &scene, float stepSize, int numOfSteps, float transparency, glm::vec4 pixelResult, float time) {
    if(transparency < TRANSPARENCY_THRESH) {
        return glm::vec4(glm::vec3(pixelResult), 1);
    } else if (numOfSteps > STEP_THRESH || glm::length(rayOrigin) > 30) {
        if(pixelResult.r > 1 || pixelResult.g > 1 || pixelResult.b > 1) {
            return glm::vec4(1,0,0,1);
        }
        if(glm::length(rayOrigin) > 30){
            std::cout << "Ray Escaped" << std::endl;
        }
        //return glm::vec4(1,0,0,1); // TEMP: See if rays are getting anywhere
        return glm::vec4(glm::vec3((/*BACKGROUND_COLOR*/ glm::normalize(rayDirection) * transparency) + pixelResult), 1.f);
    } else {
        std::vector<float> objectDensities(scene.getShapes().size());

        glm::vec4 rayEndpoint = rayOrigin + (stepSize * rayDirection);

        //In reality this should be the middle of the ray step, not the end, but also who cares for now

        glm::vec4 raySelectPoint = rayOrigin + (((float)rand() / ((float)RAND_MAX)) * stepSize * rayDirection);

        objectDensities = findDensityDataForShapes(scene.getShapes(), raySelectPoint);

        //select the highest density across all shapes
        int highestDensityIdx = 0;
        float highestDensity = 0.f;
        for(int i = 0; i < objectDensities.size(); i++) {
            if(objectDensities[i] > highestDensity) {
                highestDensity = objectDensities[i];
                highestDensityIdx = i;
            }
        }

        //do lighting stuff for that object

        float currentTransparency = transparency;

        if(highestDensity > 0) {

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
        }

        // Run black hole geodesic integration here

        // if(numOfSteps == 100){
        //     std::cout << "ooga booga" << std::endl;
        // }

        // Convert to spherical coordinates
        glm::vec4 posRelative = rayEndpoint - m_blackHole.position;
        //std::cout << "posRelative at: " << numOfSteps << " : " << posRelative[0] << " " << posRelative[1] << " " << posRelative[2] << " " << posRelative[3] << std::endl;
        float t = time,r = glm::length(posRelative),theta = atan(posRelative.y/posRelative.x),phi = acos(posRelative.z/r);

        glm::mat4 schwarzchildMetric = this->m_blackHole.getSchwarzchildMetric(t,r,theta,phi);
        glm::mat4 tetradBasis = m_blackHole.getTetradBasis(t,r,theta,phi);

        // for(int i = 0; i < 4; ++i){
        //     for(int j = 0; j < 4; ++j){
        //         std::cout << tetradBasis[i][j] << " ";
        //     }
        //     std::cout << std::endl;
        // }

        // // convert ray such that z points to the black hole

        glm::vec4 position = {t,r,theta,phi};
        // Weird NaN stuff sometimes happens here
        float deltaPhi = acos((posRelative.z+rayDirection.z)/r);
        if((posRelative.z+rayDirection.z)/r == 0){
            deltaPhi = 0;
        }

        glm::vec4 positionDelta = {0,glm::length(posRelative+rayDirection),atan((posRelative.y+rayDirection.y)/(posRelative.x+rayDirection.x)), deltaPhi};
        glm::vec4 direction = positionDelta-position;
        //direction = glm::normalize(direction);
        direction[0] = -1;
        //glm::vec4 direction = {-1,-rayDirection[2], rayDirection[1], rayDirection[0]};

        //std::cout << "BEGINNING: " << direction[0] << " " << direction[1] << " " << direction[2] << " " << direction[3] << std::endl;
        //std::cout << "position: " << position[0] << " " << position[1] << " " << position[2] << " " << position[3] << std::endl;
        //std::cout << "positionDelta: " << positionDelta[0] << " " << positionDelta[1] << " " << positionDelta[2] << " " << positionDelta[3] << std::endl;
        //std::cout << "posRelative: " << posRelative[0] << " " << posRelative[1] << " " << posRelative[2] << " " << posRelative[3] << std::endl;
        //std::cout << "rayDirection: " << rayDirection[0] << " " << rayDirection[1] << " " << rayDirection[2] << " " << rayDirection[3] << std::endl;

        //std::cout << direction[0] << " " << direction[1] << " " << direction[2] << " " << direction[3] << std::endl;
        direction = tetradBasis * direction;



        // Calculate new ray direction, normalize velocity to the speed of light

        std::vector<glm::mat4> christoff2 = m_blackHole.calculateSchwarzchildChristoffel(position);

        // Now that we have all neccesary parts of the geodesic equation, we can finally integrate along the geodesic... for one single step of one single ray.
        // This is a "bit" expensive
        glm::vec4 acceleration = glm::vec4(0);
        for(int mu = 0; mu < 4; ++mu){
            float sum =0.f;
            for(int al = 0; al < 4; ++al){
                for(int be = 0; be < 4; ++be){
                    sum += -christoff2[mu][al][be] * direction[al] * direction[be];
                }
            }
            acceleration[mu] = sum;
        }

        //std::cout << "MIDDLE: " << direction[0] << " " << direction[1] << " " << direction[2] << " " << direction[3] << std::endl;

        //acceleration = glm::normalize(acceleration);
        direction+=acceleration * stepSize;
        position+=direction*stepSize;

        position[0] = time;

        if(position[1] < m_blackHole.rs && position[1] > 0){
            std::cout << position[1] << " " << m_blackHole.rs << std::endl;
            return {0,0,0,1};
        }

        // If not, convert back to cartesian coordinates and step the ray.

        //std::cout << sin(position[2]) << " " << cos(position[2]) << " " << sin(position[3]) << " " << cos(position[3]) << std::endl;
        //std::cout << position[1] << std::endl;
        glm::vec4 positionCartesian = {position[1] * sin(position[3]) * cos(position[2]),
                                         position[1] * sin(position[3]) * sin(position[2]),
                                        position[1] * cos(position[3]), 1};

        glm::vec4 directionCartesian = {direction[1] * sin(direction[3]) * cos(direction[2]),
                                      direction[1] * sin(direction[3]) * sin(direction[2]),
                                      direction[1] * cos(direction[3]), 0};

        //glm::vec4 directionCartesian = positionCartesian - rayOrigin;
        directionCartesian = glm::normalize(directionCartesian);

        //directionCartesian = glm::normalize(directionCartesian)*stepSize;

        //std::cout << glm::length(directionCartesian) << std::endl;
        //std::cout << position[0] << " " << position[1] << " " << position[2] << " " << position[3] << std::endl;
        //std::cout << "END: " << direction[0] << " " << direction[1] << " " << direction[2] << " " << direction[3] << std::endl;
        //std::cout << "END (CART): " << directionCartesian[0] << " " << directionCartesian[1] << " " << directionCartesian[2] << " " << directionCartesian[3] << std::endl;

        return recurRayMarch(positionCartesian, directionCartesian, scene, stepSize, numOfSteps + 1, currentTransparency, pixelResult, time+0.05);
    }
}

//TODO:
//- better lighting effects
//- speed-ups (skip pixels with no intersects)

