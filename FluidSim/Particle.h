#pragma once
#include "raylib.h"
#include "raymath.h"

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

    float getMass() { return mass; }
    float getRadius() { return radius; }
    Vector2 getPos() { return pos; }
    Vector2 getProjectedPos(float deltaTime) { return Vector2Add(pos, Vector2Scale(vel, deltaTime)); }
    Vector2 getVel() { return vel; }
    void setPos(Vector2 newPos) { pos = newPos; }
    void setX(float x) { pos.x = x; }
    void setY(float y) { pos.y = y; }
    void setVel(Vector2 newVel) { vel = newVel; }
    void setVelX(float x) { vel.x = x; }
    void setVelY(float y) { vel.y = y; }
    
    void applyForce(Vector2 force);

    void update(float deltaTime);
    void draw();
};