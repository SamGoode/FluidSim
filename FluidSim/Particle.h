#pragma once
#include "raymath.h"

struct Particle {
    int isActive;
    float mass;
    Vector2 pos;
    Vector2 vel;
};