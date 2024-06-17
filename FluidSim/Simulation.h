#pragma once
#include "raylib.h"
#include "Particle.h"
#include "Array.h"
#include "SpatialHashGrid.h"

class Simulation {
private:
    Vector4 bounds;
    float scale;

    Vector2 gravity;
    float collisionDampening;

    bool showSmoothingRadius;
    float smoothingRadius;
    float targetDensity;
    float pressureMultiplier;
    float timeMultiplier;
    float mouseInteractRadius;
    float mouseInteractForce;

    float defaultMass;
    float defaultRadius;
    Array<Particle> particles;
    Array<Vector2> projectedPositions;
    Array<float> densities;

    SpatialHashGrid spatialHash;

    unsigned int particleUpdateProgram;
    unsigned int inSsbo;
    unsigned int outSsbo;

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

    void update(float deltaTime);
    void draw();
};
