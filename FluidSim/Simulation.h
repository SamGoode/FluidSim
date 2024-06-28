#pragma once
#include "raylib.h"
#include "Particle.h"
#include "Array.h"
#include "SpatialHashGrid.h"

struct SimData {
    Vector2 gravity;
    float targetDensity;
    float fixedTimeStep;
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
    float collisionDampening;

    bool showSmoothingRadius;
    float smoothingRadius;
    float sqrRadius;
    float particleRadius;
    float targetDensity;
    float pressureMultiplier;
    float timeDilation;

    float mouseInteractRadius;
    float mouseInteractForce;

    int spawnAmount;
    SpawnArea spawnArea;
    DespawnArea despawnArea;
    Array<AirflowSpace> airflows;
    Array<Ball> balls;

    float defaultMass;
    Array<Particle> particles;
    Array<int> objectPool;
    int activeCount;

    float fixedTimeStep;
    float timePassed;

    Array<float> densities;
    Array<Vector2> previousPositions;

    SpatialHashGrid spatialHash;
    int2 cellOffsets[9];

    unsigned int updateParticleProgram;
    unsigned int gravProjectionProgram;
    unsigned int clearTextureBufferProgram;
    unsigned int updateTextureBufferProgram;
    Shader renderSimShader;
    int resUniformLoc;

    SimData simData;
    unsigned int simDataSSBO;
    unsigned int particleSSBO;
    unsigned int projectedPositionSSBO;
    unsigned int densitySSBO;
    unsigned int textureSSBO;

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

    static float smoothing(float radius, float dist);
    static float smoothingGradient(float radius, float dist);
    float calculateDensity(int particleID);
    float convertDensityToPressure(float density);
    float calculateSharedPressure(float densityA, float densityB);
    Vector2 calculateGradientVec(int particleID);

    void spawnParticle(float mass, Vector2 pos, Vector2 vel);
    void despawnParticle(int particleID);

    void update(float deltaTime);
    void draw();
};
