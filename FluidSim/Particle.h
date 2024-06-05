#pragma once
#include "raylib.h"

class Simulation;

class Particle {
private:
    Simulation* sim;
    float mass;
    float radius;
    Vector2 pos;
    Vector2 vel;

public:
    Particle() {}

    Particle(Simulation* _sim, Vector2 _pos);

    Particle(const Particle& copy);

    Particle& operator=(const Particle& copy);

    void draw();
};