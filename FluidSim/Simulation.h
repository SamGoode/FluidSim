#pragma once
#include "raylib.h"
#include "Particle.h"
#include "Array.h"

class Simulation {
private:
    Vector4 bounds;
    Vector2 gravity;
    float collisionDampening;
    float smoothingRadius;
    Array<Particle> particles;

public:
    Simulation(Vector4 _bounds);

    Vector4 getBounds() { return bounds; }
    float getWidth() { return bounds.z - bounds.x; }
    float getHeight() { return bounds.w - bounds.y; }

    static float smoothing(float radius, float dist);
    float calculateDensity(Vector2 pos);
    Vector2 calculateGradientVec(Vector2 pos);

    void update(float deltaTime);
    void draw();
};
