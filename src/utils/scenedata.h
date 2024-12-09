#pragma once

#include <vector>
#include <string>
#include <functional>

#include <glm/glm.hpp>

#include "rgba.h"

#include "imagereader.h"

// Enum of the types of virtual lights that might be in the scene
enum class LightType {
    LIGHT_POINT,
    LIGHT_DIRECTIONAL,
    LIGHT_SPOT,
};

// Enum of the types of primitives that might be in the scene
enum class PrimitiveType {
    PRIMITIVE_CUBE,
    PRIMITIVE_CONE,
    PRIMITIVE_CYLINDER,
    PRIMITIVE_SPHERE,
    PRIMITIVE_MESH,
    PRIMITIVE_DENSITY_CUBE,
    PRIMITIVE_DENSITY_CONE,
    PRIMITIVE_DENSITY_CYLINDER,
    PRIMITIVE_DENSITY_SPHERE
};

// Enum of the types of transformations that can be applied
enum class TransformationType {
    TRANSFORMATION_TRANSLATE,
    TRANSFORMATION_SCALE,
    TRANSFORMATION_ROTATE,
    TRANSFORMATION_MATRIX
};

// Type which can be used to store an RGBA color in floats [0,1]
using SceneColor = glm::vec4;

// Struct which contains the global color coefficients of a scene.
// These are multiplied with the object-specific materials in the lighting equation.
struct SceneGlobalData {
    float ka; // Ambient term
    float kd; // Diffuse term
    float ks; // Specular term
    float kt; // Transparency; used for extra credit (refraction)
};

// Struct which contains raw parsed data fro a single light
struct SceneLight {
    int id;
    LightType type;

    SceneColor color;
    glm::vec3 function; // Attenuation function
    glm::vec4 dir;      // Not applicable to point lights

    float penumbra; // Only applicable to spot lights, in RADIANS
    float angle;    // Only applicable to spot lights, in RADIANS

    float width, height; // No longer supported (area lights)
};

// Struct which contains data for a single light with CTM applied
struct SceneLightData {
    int id;
    LightType type;

    SceneColor color;
    glm::vec3 function; // Attenuation function

    glm::vec4 pos; // Position with CTM applied (Not applicable to directional lights)
    glm::vec4 dir; // Direction with CTM applied (Not applicable to point lights)

    float penumbra; // Only applicable to spot lights, in RADIANS
    float angle;    // Only applicable to spot lights, in RADIANS

    float width, height; // No longer supported (area lights)
};

// Struct which contains data for the camera of a scene
struct SceneCameraData {
    glm::vec4 pos;
    glm::vec4 look;
    glm::vec4 up;

    float heightAngle; // The height angle of the camera in RADIANS

    float aperture;    // Only applicable for depth of field
    float focalLength; // Only applicable for depth of field
};

// Struct which contains data for texture mapping files
struct SceneFileMap {
    SceneFileMap() : isUsed(false) {}

    bool isUsed;
    std::string filename;

    float repeatU;
    float repeatV;

    void clear()
    {
        isUsed = false;
        repeatU = 0.0f;
        repeatV = 0.0f;
        filename = std::string();
    }
};

// Struct which contains data for a material (e.g. one which might be assigned to an object)
struct SceneMaterial {
    SceneColor cAmbient;  // Ambient term
    SceneColor cDiffuse;  // Diffuse term
    SceneColor cSpecular; // Specular term
    float shininess;      // Specular exponent

    SceneColor cReflective; // Used to weight contribution of reflected ray lighting (via multiplication)

    SceneColor cTransparent; // Transparency;        used for extra credit (refraction)
    float ior;               // Index of refraction; used for extra credit (refraction)

    SceneFileMap textureMap; // Used for texture mapping
    float blend;             // Used for texture mapping
    Image* texture;

    SceneColor cEmissive; // Not used
    SceneFileMap bumpMap; // Not used

    float density = 10000;

    void clear()
    {
        cAmbient = glm::vec4(0);
        cDiffuse = glm::vec4(0);
        cSpecular = glm::vec4(0);
        shininess = 0;

        cReflective = glm::vec4(0);

        cTransparent = glm::vec4(0);
        ior = 0;

        textureMap.clear();
        blend = 0;

        cEmissive = glm::vec4(0);
        bumpMap.clear();
    }
};

//Struct which contains data for some potential intersection point on a primitive
struct IntersectionData {
    bool hasIntersect = false;
    float intersectT;
    glm::vec4 normAtIntForObj;
    glm::vec2 uvCoords;
    bool solid = true;
};

//Consuming some initial point P and vector ray direction D, return some IntersectionData object for the ScenePrimitive
using ImplicitEquationSignature = auto(glm::vec4 p, glm::vec4 d)->IntersectionData;

//Consuming some point P, return some transmission factor (sigma_a) for the ScenePrimitive at that point. Larger values = high density objects, 0 = transparent
using ImplicitDensitySignature = auto(glm::vec4 p)->float;

// Struct which contains data for a single primitive in a scene
struct ScenePrimitive {
    PrimitiveType type;
    SceneMaterial material;
    std::string meshfile; // Used for triangle meshes
    std::function<ImplicitEquationSignature> getIntersectData;
    std::function<ImplicitDensitySignature> getDensityData;
};


// Struct which contains data for a transformation.
struct SceneTransformation {
    TransformationType type;

    glm::vec3 translate; // Only applicable when translating. Defines t_x, t_y, and t_z, the amounts to translate by, along each axis.
    glm::vec3 scale;     // Only applicable when scaling.     Defines s_x, s_y, and s_z, the amounts to scale by, along each axis.
    glm::vec3 rotate;    // Only applicable when rotating.    Defines the axis of rotation; should be a unit vector.
    float angle;         // Only applicable when rotating.    Defines the angle to rotate by in RADIANS, following the right-hand rule.
    glm::mat4 matrix;    // Only applicable when transforming by a custom matrix. This is that custom matrix.
};

// Struct which represents a node in the scene graph/tree, to be parsed by the student's `SceneParser`.
struct SceneNode {
    std::vector<SceneTransformation*> transformations; // Note the order of transformations described in lab 5
    std::vector<ScenePrimitive*> primitives;
    std::vector<SceneLight*> lights;
    std::vector<SceneNode*> children;
};
