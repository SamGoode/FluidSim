#include "Particle.h"
#include "Simulation.h"
#include "raymath.h"

Particle::Particle(Simulation* _sim, Vector2 _pos) {
    sim = _sim;
    mass = 1;
    radius = 5;
    pos = _pos;
    vel = { 0, 0 };
}

Particle::Particle(const Particle& copy) {
    sim = copy.sim;
    mass = copy.mass;
    radius = copy.radius;
    pos = copy.pos;
    vel = copy.vel;
}

Particle& Particle::operator=(const Particle& copy) {
    sim = copy.sim;
    mass = copy.mass;
    radius = copy.radius;
    pos = copy.pos;
    vel = copy.vel;

    return *this;
}

void Particle::applyForce(Vector2 force) {
    vel = Vector2Add(vel, Vector2Divide(force, mass));
}

void Particle::update(float deltaTime) {
    pos = Vector2Add(pos, Vector2Scale(vel, deltaTime));
}

void Particle::draw() {
    Vector4 simBounds = sim->getBounds();
    DrawCircle(simBounds.x + pos.x, simBounds.y + pos.y, radius, BLUE);
}