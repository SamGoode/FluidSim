#pragma once
#include "raylib.h"
#include "Particle.h"
#include "Array.h"

class Simulation {
private:
    Vector4 bounds;
    Vector2 gravity;
    Array<Particle> particles;

public:
    Simulation(Vector4 _bounds);

    Vector4 getBounds();

    void draw();
};
