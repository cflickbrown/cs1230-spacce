#pragma once

#include <vector>
#include <functional>
#include <algorithm>
#include <cmath>

#include <glm/glm.hpp>

#include "scenedata.h"
#include "perlinnoise.h"

//absurdly large alpha to simulate solid-ness
float SOLID_ALPHA = 101;

PerlinNoise P_NOISE = PerlinNoise();

ScenePrimitive SolidCube(PrimitiveType type,
                         SceneMaterial material,
                         std::string meshfile) {
    return ScenePrimitive {
        type,
        material,
        meshfile,
        [=](glm::vec4 p, glm::vec4 d) {
            return IntersectionData();
        },
        [=](glm::vec4 p) {

            //if point outside the object, then object does not add anything
            if(p.x > 0.5 || p.x < -0.5 || p.y > 0.5 || p.y < -0.5 || p.z > 0.5 || p.z < -0.5) {
                return 0.f;
            }

            return SOLID_ALPHA;
        }
    };
}

ScenePrimitive SolidSphere(PrimitiveType type,
                         SceneMaterial material,
                         std::string meshfile) {
    return ScenePrimitive {
        type,
        material,
        meshfile,
        [=](glm::vec4 p, glm::vec4 d) {
            return IntersectionData();
        },
        [=](glm::vec4 p) {

            //if point outside the object, then object does not add anything
            if(sqrt((p.x * p.x) + (p.y * p.y) + (p.z * p.z)) > 0.5) {
                return 0.f;
            }

            return SOLID_ALPHA;
        }
    };
}

ScenePrimitive SolidCylinder(PrimitiveType type,
                           SceneMaterial material,
                           std::string meshfile) {
    return ScenePrimitive {
        type,
        material,
        meshfile,
        [=](glm::vec4 p, glm::vec4 d) {
            return IntersectionData();
        },
        [=](glm::vec4 p) {

            //if point outside the object, then object does not add anything
            if(sqrt((p.x * p.x) + (p.z * p.z)) > 0.5 || p.y < -0.5 || p.y > 0.5) {
                return 0.f;
            }

            return SOLID_ALPHA;
        }
    };
}

ScenePrimitive SolidCone(PrimitiveType type,
                             SceneMaterial material,
                             std::string meshfile) {
    return ScenePrimitive {
        type,
        material,
        meshfile,
        [=](glm::vec4 p, glm::vec4 d) {
            return IntersectionData();
        },
        [=](glm::vec4 p) {

            //if point outside the object, then object does not add anything
            if(sqrt((p.x * p.x) + (p.z * p.z)) > (0.5 - ((p.y + 0.5)/2)) || p.y < -0.5 || p.y > 0.5) {
                return 0.f;
            }

            return SOLID_ALPHA;
        }
    };
}


// Implict objects that return ScenePrimitives with intersect functions

ScenePrimitive Cube(PrimitiveType type,
          SceneMaterial material,
          std::string meshfile) {
    return ScenePrimitive{
        type,
        material,
        meshfile,
        [=](glm::vec4 p, glm::vec4 d) {

            float topXFaceT = (0.5 - p.x) / d.x;
            float topYFaceT = (0.5 - p.y) / d.y;
            float topZFaceT = (0.5 - p.z) / d.z;
            float bottomXFaceT = (-0.5 - p.x) / d.x;
            float bottomYFaceT = (-0.5 - p.y) / d.y;
            float bottomZFaceT = (-0.5 - p.z) / d.z;

            //get face intersects & uv coordinates
            IntersectionData topXFace = IntersectionData(false, topXFaceT, glm::vec4(1,0,0,0), glm::vec2(-1.0 * ((p.z + topXFaceT * d.z) - 0.5), (p.y + topXFaceT * d.y) + 0.5));
            IntersectionData bottomXFace = IntersectionData(false, bottomXFaceT, glm::vec4(-1,0,0,0), glm::vec2((p.z + bottomXFaceT * d.z) + 0.5, (p.y + bottomXFaceT * d.y) + 0.5));
            IntersectionData topYFace = IntersectionData(false, topYFaceT, glm::vec4(0,1,0,0), glm::vec2((p.x + topYFaceT * d.x) + 0.5, -1.0 * ((p.z + topYFaceT * d.z) + 0.5)));
            IntersectionData bottomYFace = IntersectionData(false, bottomYFaceT, glm::vec4(0,-1,0,0), glm::vec2((p.x + bottomYFaceT * d.x) + 0.5, (p.z + bottomYFaceT * d.z) + 0.5));
            IntersectionData topZFace = IntersectionData(false, topZFaceT, glm::vec4(0,0,1,0), glm::vec2((p.x + topZFaceT * d.x) + 0.5, (p.y + topZFaceT * d.y) + 0.5));
            IntersectionData bottomZFace = IntersectionData(false, bottomZFaceT, glm::vec4(0,0,-1,0), glm::vec2(-1.0*((p.x + bottomZFaceT * d.x) - 0.5), (p.y + bottomZFaceT * d.y) + 0.5));

            std::vector<IntersectionData> intersects = {topXFace, bottomXFace, topYFace, bottomYFace, topZFace, bottomZFace};

            std::vector<IntersectionData> validIntersects;

            float compPoint = 0.500005; //accounting for some float error

            //make sure intersects are on the cube
            for (IntersectionData intersect : intersects) {
                if (p.x + intersect.intersectT * d.x <= compPoint && p.x + intersect.intersectT * d.x >= -compPoint &&
                    p.y + intersect.intersectT * d.y <= compPoint && p.y + intersect.intersectT * d.y >= -compPoint &&
                    p.z + intersect.intersectT * d.z <= compPoint && p.z + intersect.intersectT * d.z >= -compPoint) {
                    validIntersects.push_back(IntersectionData(true, intersect.intersectT, intersect.normAtIntForObj, intersect.uvCoords, material.solid));
                } else {
                    validIntersects.push_back(IntersectionData(false, FLT_MAX, glm::vec4(0)));
                }
            }

            //select nearest intersect
            IntersectionData smallestIntersect = IntersectionData(false, FLT_MAX, glm::vec4(0));

            for(IntersectionData validInt : validIntersects) {
                if(smallestIntersect.intersectT > validInt.intersectT) {
                    smallestIntersect = validInt;
                }
            }

            return smallestIntersect;
        },
        [=](glm::vec4 p) {

            //if point outside the object, then object does not add anything
            if(p.x > 0.5 || p.x < -0.5 || p.y > 0.5 || p.y < -0.5 || p.z > 0.5 || p.z < -0.5) {
                return 0.f;
            }

            return material.density;
        }
    };
}


ScenePrimitive Sphere(PrimitiveType type,
          SceneMaterial material,
          std::string meshfile) {
    return ScenePrimitive{
        type,
        material,
        meshfile,
        [=](glm::vec4 p, glm::vec4 d) {

            float valA = d.x * d.x + d.y * d.y + d.z * d.z;
            float valB = 2.0 * (p.x * d.x + p.y * d.y + p.z * d.z);
            float valC = p.x * p.x + p.y * p.y + p.z * p.z - 0.25;
            //0.25 = sqrt(radius = 0.5);

            float discrim = valB * valB - (4.0 * valA * valC);

            if(discrim > 0) {
                float root1 = (-1.0 * valB + sqrt(discrim)) / (2.0 * valA);
                float root2 = (-1.0 * valB - sqrt(discrim)) / (2.0 * valA);
                float useRoot = fmin(root1, root2);

                glm::vec4 normAtPoint = p + useRoot * d;
                normAtPoint *= 2;

                float phiV = asin((p.y + useRoot * d.y) / 0.5);
                float vCoord = (phiV / M_PI) + 0.5;

                float thetaU = atan2(p.z + useRoot * d.z, p.x + useRoot * d.x);
                float uCoord = thetaU < 0 ? -1*thetaU/(2*M_PI) : 1 - (thetaU/(2*M_PI));

                if(vCoord == 0 || vCoord == 1) {
                    uCoord = 0.5;
                }

                return IntersectionData(true, useRoot, glm::normalize(normAtPoint), glm::vec2(uCoord, vCoord), material.solid);
            } else if (discrim == 0.0) {
                float useRoot = (-1.0 * valB + sqrt(discrim)) / (2.0 * valA);

                glm::vec4 normAtPoint = p + useRoot * d;
                normAtPoint *= 2;

                float phiV = asin((p.y + useRoot * d.y) / 0.5);
                float vCoord = (phiV / M_PI) + 0.5;

                float thetaU = atan2(p.z + useRoot * d.z, p.x + useRoot * d.x);
                float uCoord = thetaU < 0 ? -1*thetaU/(2*M_PI) : 1 - (thetaU/(2*M_PI));

                if(vCoord == 0 || vCoord == 1) {
                    uCoord = 0.5;
                }

                return IntersectionData(true, useRoot, glm::normalize(normAtPoint), glm::vec2(uCoord, vCoord), material.solid);
            } else {
                return IntersectionData(false, FLT_MAX, glm::vec4(0));
            }
        },
        [=](glm::vec4 p) {

            //if point outside the object, then object does not add anything
            if(sqrt((p.x * p.x) + (p.y * p.y) + (p.z * p.z)) > 0.5) {
                return 0.f;
            }

            // return P_NOISE.computePerlin3d(p.x, p.y, p.z);
            return material.density;
        }
    };
}

ScenePrimitive Cylinder(PrimitiveType type,
            SceneMaterial material,
            std::string meshfile) {
    return ScenePrimitive{
        type,
        material,
        meshfile,
        [=](glm::vec4 p, glm::vec4 d) {

            std::vector<IntersectionData> intersects;

            float valASide = d.x * d.x + d.z * d.z;
            float valBSide = 2.0 * (p.x * d.x + p.z * d.z);
            float valCSide = p.x * p.x + p.z * p.z - 0.25;
            //0.25 = sqrt(radius = 0.5);

            float discrim = valBSide * valBSide - 4.0 * valASide * valCSide;

            float compPoint = 0.50001; //accounting for some float error

            //get wall intersects
            if(discrim > 0.0) {
                float root1 = (-1.0 * valBSide + sqrt(discrim)) / (2.0 * valASide);
                float root2 = (-1.0 * valBSide - sqrt(discrim)) / (2.0 * valASide);

                if(p.y + root1 * d.y <= compPoint && p.y + root1 * d.y >= -compPoint) {
                    glm::vec4 normVec = p + root1 * d;
                    normVec *= 2;
                    normVec.y = 0;

                    float thetaU = atan2(p.z + root1 * d.z, p.x + root1 * d.x);
                    float uCoord = thetaU < 0 ? -1*thetaU/(2*M_PI) : 1 - (thetaU/(2*M_PI));

                    intersects.push_back(IntersectionData(true, root1, glm::normalize(normVec), glm::vec2(uCoord, 0.5 + p.y + root1 * d.y), material.solid));
                }

                if(p.y + root2 * d.y <= compPoint && p.y + root2 * d.y >= -compPoint) {
                    glm::vec4 normVec = p + root2 * d;
                    normVec *= 2;
                    normVec.y = 0;

                    float thetaU = atan2(p.z + root2 * d.z, p.x + root2 * d.x);
                    float uCoord = thetaU < 0 ? -1*thetaU/(2*M_PI) : 1 - (thetaU/(2*M_PI));

                    intersects.push_back(IntersectionData(true, root2, glm::normalize(normVec), glm::vec2(uCoord, 0.5 + p.y + root2 * d.y), material.solid));
                }
            } else if (discrim == 0.0) {
                float root1 = (-1.0 * valBSide + sqrt(discrim)) / (2.0 * valASide);

                if(p.y + root1 * d.y <= compPoint && p.y + root1 * d.y >= -compPoint) {
                    glm::vec4 normVec = p + root1 * d;
                    normVec *= 2;
                    normVec.y = 0;

                    float thetaU = atan2(p.z + root1 * d.z, p.x + root1 * d.x);
                    float uCoord = thetaU < 0 ? -1*thetaU/(2*M_PI) : 1 - (thetaU/(2*M_PI));

                    intersects.push_back(IntersectionData(true, root1, glm::normalize(normVec), glm::vec2(uCoord, 0.5 + p.y + root1 * d.y), material.solid));
                }
            } else {
                intersects.push_back(IntersectionData(false, FLT_MAX, glm::vec4(0)));
            }

            //get cap intersects
            float topCapT = (0.5 - p.y) / d.y;
            float bottomCapT = (-0.5 - p.y) / d.y;

            if(pow(p.x + topCapT * d.x, 2) + pow(p.z + topCapT * d.z, 2) <= 0.25) {
                intersects.push_back(IntersectionData(true, topCapT, glm::vec4(0, 1, 0, 0), glm::vec2((p.x + topCapT * d.x) + 0.5, -1.0 * ((p.z + topCapT * d.z) - 0.5)), material.solid));
            }

            if(pow(p.x + bottomCapT * d.x, 2) + pow(p.z + bottomCapT * d.z, 2) <= 0.25) {
                intersects.push_back(IntersectionData(true, bottomCapT, glm::vec4(0, -1, 0, 0), glm::vec2((p.x + bottomCapT * d.x) - 0.5, (p.z + bottomCapT * d.z) - 0.5), material.solid));
            }

            IntersectionData smallestIntersect = IntersectionData(false, FLT_MAX, glm::vec4(0));

            //get the nearest intersect
            for(IntersectionData intersect : intersects) {
                if(intersect.hasIntersect && intersect.intersectT < smallestIntersect.intersectT) {
                    smallestIntersect = intersect;
                }
            }

            return smallestIntersect;
        },
        [=](glm::vec4 p) {

            //if point outside the object, then object does not add anything
            if(sqrt((p.x * p.x) + (p.z * p.z)) > 0.5 || p.y < -0.5 || p.y > 0.5) {
                return 0.f;
            }

            // return (float)rand()/((float)RAND_MAX / 1.f);
            return material.density;
        }
    };
}

ScenePrimitive Cone(PrimitiveType type,
              SceneMaterial material,
              std::string meshfile) {
    return ScenePrimitive{
        type,
        material,
        meshfile,
        [=](glm::vec4 p, glm::vec4 d) {

            std::vector<IntersectionData> intersects;

            float valA = d.x * d.x + d.z * d.z - 0.25 * (d.y * d.y);
            float valB = 2.0 * p.x * d.x + 2.0 * p.z * d.z - 0.5 * (p.y * d.y) + 0.25 * d.y;
            float valC = p.x * p.x + p.z * p.z - (0.25 * (p.y * p.y)) + (0.25 * p.y) - 0.0625;

            float discrim = valB * valB - 4.0 * valA * valC;

            float compPoint = 0.5; //accounting for some float error

            //get cone top intersects
            if(discrim > 0.0) {
                float root1 = (-1.0 * valB + sqrt(discrim)) / (2.0 * valA);
                float root2 = (-1.0 * valB - sqrt(discrim)) / (2.0 * valA);

                if((p.y + root1 * d.y) <= compPoint && (p.y + root1 * d.y) >= -compPoint) {
                    glm::vec4 normVec = p + root1 * d;
                    normVec *= 2;
                    normVec.y *= compPoint; //comes out to just y in partial derivative
                    normVec.y = (compPoint - normVec.y) / 2.0;

                    float thetaU = atan2(p.z + root1 * d.z, p.x + root1 * d.x);
                    float uCoord = thetaU < 0 ? -1*thetaU/(2*M_PI) : 1 - (thetaU/(2*M_PI));

                    intersects.push_back(IntersectionData(true, root1, glm::normalize(normVec), glm::vec2(uCoord, 0.5 + p.y + root1 * d.y), material.solid));
                }

                if((p.y + root2 * d.y <= compPoint) && (p.y + root2 * d.y) >= -compPoint) {
                    glm::vec4 normVec = p + root2 * d;
                    normVec *= 2;
                    normVec.y *= 0.5;
                    normVec.y = (0.5 - normVec.y) / 2.0;

                    float thetaU = atan2(p.z + root2 * d.z, p.x + root2 * d.x);
                    float uCoord = thetaU < 0 ? -1*thetaU/(2*M_PI) : 1 - (thetaU/(2*M_PI));

                    intersects.push_back(IntersectionData(true, root2, glm::normalize(normVec), glm::vec2(uCoord, 0.5 + p.y + root2 * d.y), material.solid));
                }
            } else if (discrim == 0.0) {
                float root1 = (-1.0 * valB + sqrt(discrim)) / (2.0 * valA);

                if((p.y + root1 * d.y) <= compPoint && (p.y + root1 * d.y) >= -compPoint) {
                    glm::vec4 normVec = p + root1 * d;
                    normVec *= 2;
                    normVec.y *= 0.5; //comes out to just y in partial derivative
                    normVec.y = (0.5 - normVec.y) / 2.0;

                    float thetaU = atan2(p.z + root1 * d.z, p.x + root1 * d.x);
                    float uCoord = thetaU < 0 ? -1*thetaU/(2*M_PI) : 1 - (thetaU/(2*M_PI));

                    intersects.push_back(IntersectionData(true, root1, glm::normalize(normVec), glm::vec2(uCoord, 0.5 + p.y + root1 * d.y), material.solid));
                }
            } else {
                intersects.push_back(IntersectionData(false, FLT_MAX, glm::vec4(0)));
            }

            //get cone base intersect
            float baseT = (-0.5 - p.y) / d.y;

            if(pow(p.x + baseT * d.x, 2) + pow(p.z + baseT * d.z, 2) <= 0.25) {
                intersects.push_back(IntersectionData(true, baseT, glm::vec4(0, -1, 0, 0), glm::vec2((p.x + baseT * d.x) + 0.5, (p.z + baseT * d.z) + 0.5), material.solid));
            }

            IntersectionData smallestIntersect = IntersectionData(false, FLT_MAX, glm::vec4(0));

            //select the nearest intersect
            for(IntersectionData intersect : intersects) {
                if(intersect.hasIntersect == true && intersect.intersectT < smallestIntersect.intersectT) {
                    smallestIntersect = intersect;
                }
            }

            return smallestIntersect;
        },
        [=](glm::vec4 p) {

            //if point outside the object, then object does not add anything
            if(sqrt((p.x * p.x) + (p.z * p.z)) > (0.5 - ((p.y + 0.5)/2)) || p.y < -0.5 || p.y > 0.5) {
                return 0.f;
            }

            return material.density;
        }
    };
}



ScenePrimitive DensityCube(PrimitiveType type,
                           SceneMaterial material,
                           std::string meshfile) {
    return ScenePrimitive {
        type,
        material,
        meshfile,
        [=](glm::vec4 p, glm::vec4 d) {
            return IntersectionData();
        },
        [=](glm::vec4 p) {
            //if point outside the object, then object does not add anything
            if(p.x > 0.5 || p.x < -0.5 || p.y > 0.5 || p.y < -0.5 || p.z > 0.5 || p.z < -0.5) {
                return 0.f;
            }

            return material.density;
        }
    };
}


ScenePrimitive DensitySphere(PrimitiveType type,
                             SceneMaterial material,
                             std::string meshfile) {
    return ScenePrimitive {
        type,
        material,
        meshfile,
        [=](glm::vec4 p, glm::vec4 d) {
            return Sphere(type, material, meshfile).getIntersectData(p, d);;
        },
        [=](glm::vec4 p) {

            //if point outside the object, then object does not add anything
            float dist = sqrt((p.x * p.x) + (p.y * p.y) + (p.z * p.z));

            if(dist >= 0.5) {
                return 0.f;
            }

            return P_NOISE.computePerlin3d(p.x, p.y, p.z);
        }
    };
}

ScenePrimitive DensityCylinder(PrimitiveType type,
                               SceneMaterial material,
                               std::string meshfile) {
    return ScenePrimitive {
        type,
        material,
        meshfile,
        [=](glm::vec4 p, glm::vec4 d) {
            return IntersectionData();
        },
        [=](glm::vec4 p) {
            //if point outside the object, then object does not add anything
            if(sqrt((p.x * p.x) + (p.z * p.z)) > 0.5 || p.y < -0.5 || p.y > 0.5) {
                return 0.f;
            }

            return material.density;
        }
    };
}

ScenePrimitive DensityCone(PrimitiveType type,
                           SceneMaterial material,
                           std::string meshfile) {
    return ScenePrimitive {
        type,
        material,
        meshfile,
        [=](glm::vec4 p, glm::vec4 d) {
            return IntersectionData();
        },
        [=](glm::vec4 p) {

            //if point outside the object, then object does not add anything
            if(sqrt((p.x * p.x) + (p.z * p.z)) > (0.5 - ((p.y + 0.5)/2)) || p.y < -0.5 || p.y > 0.5) {
                return 0.f;
            }

            return material.density;
        }
    };
}



