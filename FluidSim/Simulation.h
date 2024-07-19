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

struct Spring {
    bool isActive;
    float restLength;
};

class SpringBuffer {
private:
    int springCount;
    Spring* springs;

public:
    SpringBuffer() {}

    SpringBuffer(int _springCount) {
        springCount = _springCount;
        springs = new Spring[springCount * springCount];
    }

    SpringBuffer(const SpringBuffer& copy) {
        springCount = copy.springCount;
        int capacity = springCount * springCount;
        springs = new Spring[capacity];
        for (int i = 0; i < capacity; i++) {
            springs[i] = copy.springs[i];
        }
    }

    ~SpringBuffer() {
        delete[] springs;
    }

    SpringBuffer& operator=(const SpringBuffer& copy) {
        delete[] springs;

        springCount = copy.springCount;
        int capacity = springCount * springCount;
        springs = new Spring[capacity];
        for (int i = 0; i < capacity; i++) {
            springs[i] = copy.springs[i];
        }

        return *this;
    }

    // capacity = (springCount^2 - springCount) / 2

    // 01 02 03 04
    // 12 13 14
    // 23 24
    // 34

    // 00 01 02 03 04
    // 10 11 12 13 14
    // 20 21 22 23 24
    // 30 31 32 33 34
    // 40 41 42 43 44

    // 01
    // 02 12
    // 03 13 23
    // 04 14 24 34

    // 00 10 20 30 40
    // 01 11 21 31 41
    // 02 12 22 32 42
    // 03 13 23 33 43
    // 04 14 24 34 44

    int getCount() {
        return springCount;
    }

    int getCapacity() {
        return springCount * springCount;
    }

    int2 getIndices(int index) {
        return { index % springCount, index / springCount };
    }

    Spring& getSpring(int index) {
        if (index >= springCount * springCount) {
            throw "index out of range";
        }
        return springs[index];
    }

    Spring& getSpring(int particleIndexA, int particleIndexB) {
        if (particleIndexA >= particleIndexB) {
            throw "invalid pairing";
        }
        if (particleIndexA >= springCount || particleIndexB >= springCount) {
            throw "particleIndex out of range";
        }

        int springIndex = particleIndexA + particleIndexB * springCount;
        return springs[springIndex];
    }
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

    //Array<Spring> springs;
    SpringBuffer springs;

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

    // particle data
    Array<float> densities;
    Array<float> nearDensities;
    Array<Vector2> previousPositions;
    Array<Vector2> velocities;
    Array<Vector2> positions;
    Array<float> masses;

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
    static float densityGradient(float radius, float dist);
    static float nearDensityKernel(float radius, float dist);
    static float nearDensityGradient(float radius, float dist);
    float calculateDensity(int particleID);
    float calculateNearDensity(int particleID);
    float convertDensityToPressure(float density);
    float calculateSharedPressure(float densityA, float densityB);
    Vector2 calculatePressureForce(int particleID);

    void spawnParticle(float mass, Vector2 pos, Vector2 vel);
    void despawnParticle(int particleID);

    void update(float deltaTime);
    void stepForward(float timeStep);
    void draw();
};
