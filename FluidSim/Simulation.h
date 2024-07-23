#pragma once
#include "raylib.h"
#include "Array.h"
#include "SpatialHashGrid.h"

struct SimData {
    Vector2 gravity;
    float targetDensity;
    float fixedTimeStep;
    int activeCount;
    float particleRadius;
};

struct SpawnArea {
    Vector4 bounds;
    Vector2 initVel;

    Vector2 getStartPos() {
        return { bounds.x, bounds.y };
    }

    Vector2 getEndPos() {
        return { bounds.z, bounds.w };
    }

    float getWidth() {
        return bounds.z - bounds.x;
    }

    float getHeight() {
        return bounds.w - bounds.y;
    }
};

struct DespawnArea {
    Vector4 bounds;

    Vector2 getStartPos() {
        return { bounds.x, bounds.y };
    }

    Vector2 getEndPos() {
        return { bounds.z, bounds.w };
    }

    float getWidth() {
        return bounds.z - bounds.x;
    }

    float getHeight() {
        return bounds.w - bounds.y;
    }
};

struct AirflowSpace {
    Vector2 pos;
    float width;
    float height;
    Vector2 force;
};

struct Ball {
    Vector2 pos;
    float radius;
};

class Simulation {
private:
    Vector4 bounds;
    float scale;
    Vector2 resolution;
    Texture texture;

    Vector2 gravity;
    float frictionCoefficient;
    float stickyDist;
    float stickyCoefficient;

    float viscLinear;
    float viscQuad;

    bool showSmoothingRadius;
    float smoothingRadius;
    float sqrRadius;
    float particleRadius;
    float targetDensity;
    float pressureMultiplier;
    float nearPressureMultiplier;
    float timeDilation;

    float mouseInteractRadius;
    float mouseInteractForce;

    int spawnAmount;
    SpawnArea spawnArea;
    DespawnArea despawnArea;
    Array<AirflowSpace> airflows;
    Array<Ball> balls;

    float defaultMass;
    Array<int> objectPool;
    int activeCount;

    float fixedTimeStep;
    float timePassed;
    float timeStepUpperBound;
    float timeStepLowerBound;

    // particle data
    Array<float> densities;
    Array<float> nearDensities;
    Array<Vector2> previousPositions;
    Array<Vector2> velocities;
    Array<Vector2> positions;
    Array<float> masses;

    float maxVelocity;

    SpatialHashGrid spatialHash;
    int2 cellOffsets[9];
    SpatialHashGrid spatialRendering;
    Array<Vector2> unscaledPositions;

    unsigned int updateParticleProgram;
    unsigned int gravProjectionProgram;
    unsigned int clearTextureBufferProgram;
    unsigned int updateTextureBufferProgram;
    Shader renderSimShader;
    int resUniformLoc;

    SimData simData;
    unsigned int simDataSSBO;
    unsigned int poolSSBO;
    unsigned int positionSSBO;
    unsigned int projectedPositionSSBO;
    unsigned int densitySSBO;

    unsigned int defaultTextureSSBO;
    unsigned int textureSSBO;

    unsigned int hashListSSBO;
    unsigned int lookupSSBO;

public:
    Simulation(Vector4 _bounds);
    ~Simulation();

    Vector4 getBounds() { return bounds; }
    float getWidth() { return bounds.z - bounds.x; }
    float getHeight() { return bounds.w - bounds.y; }
    float getScaledWidth() { return (bounds.z - bounds.x)/scale; }
    float getScaledHeight() { return (bounds.w - bounds.y)/scale; }
    bool outOfBounds(Vector2 pos) { return pos.x < 0 || pos.x >= bounds.z - bounds.x || pos.y < 0 || pos.y >= bounds.w - bounds.y; }
    Vector2 convertToSimPos(Vector2 screenPos) { return { (screenPos.x - bounds.x) / scale, (screenPos.y - bounds.y) / scale }; }
    Vector2 convertToScreenPos(Vector2 simPos) { return { bounds.x + (simPos.x * scale), bounds.y + (simPos.y * scale) }; }

    static float densityKernel(float radius, float dist);
    static float nearDensityKernel(float radius, float dist);
    std::pair<float, float> calculateDensity(int particleID);
    float convertDensityToPressure(float density);
    Vector2 calculatePressureForce(int particleID);
    void applyPressureDisplacements(int particleID, float deltaTime);
    Vector2 calculateViscosityImpulse(int particleID);

    void spawnParticle(float mass, Vector2 pos, Vector2 vel);
    void despawnParticle(int particleID);

    float calculateTimeStepSize();
    void update(float deltaTime);
    void stepForward();
    void draw();
};
