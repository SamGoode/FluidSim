#include "Simulation.h"
#include "raymath.h"

Simulation::Simulation(Vector4 _bounds) {
    bounds = _bounds;

    gravity = { 0, 0 };
    collisionDampening = 0.2;

<<<<<<< HEAD
<<<<<<< HEAD
    showSmoothingRadius = false;
    smoothingRadius = 40;
    pressureMultiplier = 10000;
=======
=======
>>>>>>> 64af7581f75ac825354e0cfd2274aa807470c69c
    smoothingRadius = 30;
    pressureMultiplier = 1000;
>>>>>>> 64af7581f75ac825354e0cfd2274aa807470c69c

    defaultMass = 1;
    defaultRadius = 5;
    //particles = Array<Particle>(256);
    //for (int i = 0; i < particles.getCount(); i++) {
    //    // generate random location
    //    float x = (rand() * (bounds.z - bounds.x)) / RAND_MAX;
    //    float y = (rand() * (bounds.w - bounds.y)) / RAND_MAX;
    //    Vector2 position = { x, y };

<<<<<<< HEAD
<<<<<<< HEAD
    //    // generate random velocity
    //    //float velX = ((rand() * 40) / RAND_MAX) - 20;
    //    //float velY = ((rand() * 40) / RAND_MAX) - 20;
    //    //Vector2 vel = { velX, velY };
    //    Vector2 vel = { 0, 0 };

    //    particles[i] = { defaultMass, defaultRadius, position, vel };
=======
=======
>>>>>>> 64af7581f75ac825354e0cfd2274aa807470c69c
    //    // initiate particle with random velocity
    //    float velX = ((rand() * 40) / RAND_MAX) - 20;
    //    float velY = ((rand() * 40) / RAND_MAX) - 20;
    //    particles[i].setVel({ velX, velY });
<<<<<<< HEAD
>>>>>>> 64af7581f75ac825354e0cfd2274aa807470c69c
=======
>>>>>>> 64af7581f75ac825354e0cfd2274aa807470c69c
    //}

    particles = Array<Particle>(400);
    // initiate particles in a square formation
    for (int i = 0; i < particles.getCount(); i++) {
<<<<<<< HEAD
<<<<<<< HEAD
        particles[i] = { 1, 5, {100 + (float)(i % 20) * 10, 100 + (float)(i / 20) * 10}, {0, 0} };
=======
        particles[i] = Particle(this, { 100 + (float)(i % 20) * 10, 100 + (float)(i / 20) * 10 });
>>>>>>> 64af7581f75ac825354e0cfd2274aa807470c69c
=======
        particles[i] = Particle(this, { 100 + (float)(i % 20) * 10, 100 + (float)(i / 20) * 10 });
>>>>>>> 64af7581f75ac825354e0cfd2274aa807470c69c
    }

    densities = Array<float>(particles.getCount());
}

// smoothing function
float Simulation::smoothing(float radius, float dist) {
<<<<<<< HEAD
<<<<<<< HEAD
    if (dist >= radius) {
        return 0;
    }

    float volume = pow(radius, 5) * PI * 0.1;
    float difference = std::max(0.f, radius - dist);
    float value = difference * difference * difference;
    return value / volume;
}

//(r-d)^3
float Simulation::smoothingGradient(float radius, float dist) {
    if (dist >= radius) {
        return 0;
    }
    
    float difference = radius - dist;
    float gradient = -(30 * difference * difference) / (PI * pow(radius, 5));
    return gradient;
=======
    float volume = PI * 0.5 * radius * radius * radius * radius;
    float value = std::max(0.f, radius - dist);
    return (value * value * value)/volume;
>>>>>>> 64af7581f75ac825354e0cfd2274aa807470c69c
=======
    float volume = PI * 0.5 * radius * radius * radius * radius;
    float value = std::max(0.f, radius - dist);
    return (value * value * value)/volume;
>>>>>>> 64af7581f75ac825354e0cfd2274aa807470c69c
}

// calculates density level at a particle's position
float Simulation::calculateDensity(Vector2 pos) {
    float density = 0;

    for (int i = 0; i < particles.getCount(); i++) {
        float dist = Vector2Distance(pos, particles[i].pos);
        density += smoothing(smoothingRadius, dist) * particles[i].mass;
    }

    return density;
}

// calculates gradient vector within the density field at a particle's position
Vector2 Simulation::calculateGradientVec(Vector2 pos) {
    Vector2 gradientVec = { 0, 0 };
    for (int i = 0; i < particles.getCount(); i++) {
        Vector2 pointToParticle = Vector2Subtract(particles[i].pos, pos);
        float dist = Vector2Length(pointToParticle);

        // skip if particle is right on top of location
        if (dist == 0 || dist >= smoothingRadius) { continue; }

        // normalised direction vector
        Vector2 dir = Vector2Divide(pointToParticle, dist);
        // magnitude of gradient vector
        float magnitude = smoothingGradient(smoothingRadius, dist);

        // because of how we compute gradient vector the direction automatically descends
        gradientVec = Vector2Add(gradientVec, Vector2Scale(dir, magnitude));
    }
    return gradientVec;
}

//force = change in mass * velocity / change in time
void Simulation::update(float deltaTime) {
    // update densities cache
    for (int i = 0; i < particles.getCount(); i++) {
        densities[i] = calculateDensity(particles[i].pos);
    }

    for (int i = 0; i < particles.getCount(); i++) {
        // gravity application
        //particles[i].applyForce(Vector2Scale(gravity, particles[i].getMass() * deltaTime));

        Vector2 gradientVec = calculateGradientVec(particles[i].pos);
        // invert the gradientVec because gradient vector aims in direction of maximising function
        Vector2 forceVec = Vector2Scale(gradientVec, 1);
        particles[i].vel = Vector2Add(particles[i].vel, Vector2Scale(forceVec, pressureMultiplier));

        // boundary collisions
        if (particles[i].pos.y + particles[i].radius >= getHeight()) {
            Vector2 vel = particles[i].vel;
            particles[i].pos.y = (getHeight() - particles[i].radius);
            particles[i].vel = Vector2Scale({vel.x, -vel.y}, (1 - collisionDampening));
        }
        else if (particles[i].pos.x - particles[i].radius <= 0) {
            Vector2 vel = particles[i].vel;
            particles[i].pos.x = particles[i].radius;
            particles[i].vel = Vector2Scale({ -vel.x, vel.y }, (1 - collisionDampening));
        }
        else if (particles[i].pos.x + particles[i].radius >= getWidth()) {
            Vector2 vel = particles[i].vel;
            particles[i].pos.x = (getWidth() - particles[i].radius);
            particles[i].vel = Vector2Scale({ -vel.x, vel.y }, (1 - collisionDampening));
        }
        else if (particles[i].pos.y - particles[i].radius <= 0) {
            Vector2 vel = particles[i].vel;
            particles[i].pos.y = particles[i].radius;
            particles[i].vel = Vector2Scale({ vel.x, -vel.y }, (1 - collisionDampening));
        }

        // update particle positions
        particles[i].pos = Vector2Add(particles[i].pos, Vector2Scale(particles[i].vel, deltaTime));
    }
}

void Simulation::draw() {
    for (int i = 0; i < particles.getCount(); i++) {
        Vector2 pos = particles[i].pos;
        float scaledDensity = densities[i] * 310;
        
        Color densityColor = { 255 * scaledDensity, 0, 255 * (1 - scaledDensity), 255};
        if (showSmoothingRadius) {
            DrawCircle(bounds.x + pos.x, bounds.y + pos.y, smoothingRadius, { 0, 121, 241, 31 });
        }
        DrawCircle(bounds.x + pos.x, bounds.y + pos.y, particles[i].radius, densityColor);
    }
    DrawRectangleLines(bounds.x, bounds.y, bounds.z - bounds.x, bounds.w - bounds.y, BLACK);
}
