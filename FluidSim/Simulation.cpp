#include "Simulation.h"
#include "raymath.h"

Simulation::Simulation(Vector4 _bounds) {
    bounds = _bounds;

    gravity = { 0, 0 };
    collisionDampening = 0.2;

    smoothingRadius = 30;
    pressureMultiplier = 1000;

    //particles = Array<Particle>(512, 512);
    //for (int i = 0; i < particles.getCount(); i++) {
    //    // spawn particle at random location
    //    float x = (rand() * (bounds.z - bounds.x)) / RAND_MAX;
    //    float y = (rand() * (bounds.w - bounds.y)) / RAND_MAX;
    //    particles[i] = Particle(this, { x, y });

    //    // initiate particle with random velocity
    //    float velX = ((rand() * 40) / RAND_MAX) - 20;
    //    float velY = ((rand() * 40) / RAND_MAX) - 20;
    //    particles[i].setVel({ velX, velY });
    //}

    particles = Array<Particle>(400, 512);
    // initiate particles in a square formation
    for (int i = 0; i < particles.getCount(); i++) {
        particles[i] = Particle(this, { 100 + (float)(i % 20) * 10, 100 + (float)(i / 20) * 10 });
    }
}

// smoothing function
float Simulation::smoothing(float radius, float dist) {
    float volume = PI * 0.5 * radius * radius * radius * radius;
    float value = std::max(0.f, radius - dist);
    return (value * value * value)/volume;
}

// calculates density level at a particle's position
float Simulation::calculateDensity(Vector2 pos) {
    //Vector2 pos = particles[particleIndex].getPos();
    //if (outOfBounds(pos)) {
    //    return 0;
    //}

    float density = 0;

    for (int i = 0; i < particles.getCount(); i++) {
        float dist = Vector2Distance(pos, particles[i].getPos());
        density += smoothing(smoothingRadius, dist) * particles[i].getMass();
    }

    return density;
}

// calculates gradient vector within the density field at a particle's position
Vector2 Simulation::calculateGradientVec(Vector2 pos) {
    //Vector2 pos = particles[particleIndex]
    
    //float originalDensity = calculateDensity(pos);

    float dx = 0.001;
    float dy = 0.001;
    float dxDensity = calculateDensity({ pos.x + dx/2, pos.y }) - calculateDensity({ pos.x - dx/2, pos.y });
    float dyDensity = calculateDensity({ pos.x, pos.y + dy/2 }) - calculateDensity({ pos.x, pos.y - dy/2 });

    return { dxDensity / dx, dyDensity / dy };
}

//force = change in mass * velocity / change in time
void Simulation::update(float deltaTime) {
    for (int i = 0; i < particles.getCount(); i++) {
        // gravity application
        particles[i].applyForce(Vector2Scale(gravity, particles[i].getMass() * deltaTime));

        Vector2 gradientVec = calculateGradientVec(particles[i].getPos());
        // invert the gradientVec because gradient vector aims in direction of maximising function
        Vector2 forceVec = Vector2Scale(gradientVec, -1);
        particles[i].applyForce(Vector2Scale(forceVec, pressureMultiplier));

        // boundary collisions
        if (particles[i].getPos().y + particles[i].getRadius() >= getHeight()) {
            Vector2 vel = particles[i].getVel();
            particles[i].setY(getHeight() - particles[i].getRadius());
            particles[i].setVel(Vector2Scale({vel.x, -vel.y}, (1 - collisionDampening)));
        }
        else if (particles[i].getPos().x - particles[i].getRadius() <= 0) {
            Vector2 vel = particles[i].getVel();
            particles[i].setX(particles[i].getRadius());
            particles[i].setVel(Vector2Scale({ -vel.x, vel.y }, (1 - collisionDampening)));
        }
        else if (particles[i].getPos().x + particles[i].getRadius() >= getWidth()) {
            Vector2 vel = particles[i].getVel();
            particles[i].setX(getWidth() - particles[i].getRadius());
            particles[i].setVel(Vector2Scale({ -vel.x, vel.y }, (1 - collisionDampening)));
        }
        else if (particles[i].getPos().y - particles[i].getRadius() <= 0) {
            Vector2 vel = particles[i].getVel();
            particles[i].setY(particles[i].getRadius());
            particles[i].setVel(Vector2Scale({ vel.x, -vel.y }, (1 - collisionDampening)));
        }

        particles[i].update(deltaTime);
    }
}

void Simulation::draw() {
    for (int i = 0; i < particles.getCount(); i++) {
        Vector2 pos = particles[i].getPos();
        DrawCircle(bounds.x + pos.x, bounds.y + pos.y, smoothingRadius, { 0, 121, 241, 31 });
        particles[i].draw();
    }
    DrawRectangleLines(bounds.x, bounds.y, bounds.z - bounds.x, bounds.w - bounds.y, BLACK);
}
