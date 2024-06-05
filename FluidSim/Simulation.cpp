#include "Simulation.h"

Simulation::Simulation(Vector4 _bounds) {
    bounds = _bounds;
    particles = Array<Particle>(16, 16);

    for (int i = 0; i < particles.getCount(); i++) {
        float x = (rand() * (bounds.z - bounds.x)) / RAND_MAX;
        float y = (rand() * (bounds.w - bounds.y)) / RAND_MAX;
        particles[i] = Particle(this, { x, y });
    }
}

Vector4 Simulation::getBounds() {
    return bounds;
}

void Simulation::draw() {
    for (int i = 0; i < particles.getCount(); i++) {
        particles[i].draw();
    }
    DrawRectangleLines(bounds.x, bounds.y, bounds.z - bounds.x, bounds.w - bounds.y, BLACK);
}
