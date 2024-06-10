#pragma once
#include "raylib.h"
#include "Particle.h"
#include "Array.h"

class Simulation {
private:
    Vector4 bounds;

    Vector2 gravity;
    float collisionDampening;

    bool showSmoothingRadius;
    float smoothingRadius;
    float targetDensity;
    float pressureMultiplier;
    float timeMultiplier;

    float defaultMass;
    float defaultRadius;
    Array<Particle> particles;
    Array<float> densities;

public:
    Simulation(Vector4 _bounds);

    Vector4 getBounds() { return bounds; }
    float getWidth() { return bounds.z - bounds.x; }
    float getHeight() { return bounds.w - bounds.y; }
    bool outOfBounds(Vector2 pos) { return pos.x < 0 || pos.x >= bounds.z - bounds.x || pos.y < 0 || pos.y >= bounds.w - bounds.y; }

    static float smoothing(float radius, float dist);
    static float smoothingGradient(float radius, float dist);
    float calculateDensity(Vector2 pos);
    float convertDensityToPressure(float density);
    float calculateSharedPressure(float densityA, float densityB);
    Vector2 calculateGradientVec(Vector2 pos);

    void update(float deltaTime);
    void draw();
};
