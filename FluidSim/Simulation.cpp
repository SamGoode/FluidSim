#include "Simulation.h"
#include "raymath.h"

Simulation::Simulation(Vector4 _bounds) {
    bounds = _bounds;

    gravity = { 0, 1 };
    collisionDampening = 0.1;

    showSmoothingRadius = true;
    smoothingRadius = 10;
    targetDensity = 0.1;
    pressureMultiplier = 1000;
    timeMultiplier = 8;

    defaultMass = 10;
    defaultRadius = 2;
    //particles = Array<Particle>(400);
    //for (int i = 0; i < particles.getCount(); i++) {
    //    // generate random location
    //    float x = (rand() * (bounds.z - bounds.x)) / RAND_MAX;
    //    float y = (rand() * (bounds.w - bounds.y)) / RAND_MAX;
    //    Vector2 position = { x, y };


    //    // generate random velocity
    //    //float velX = ((rand() * 40) / RAND_MAX) - 20;
    //    //float velY = ((rand() * 40) / RAND_MAX) - 20;
    //    //Vector2 vel = { velX, velY };
    //    Vector2 vel = { 0, 0 };

    //    particles[i] = { defaultMass, defaultRadius, position, vel };
    //}

    particles = Array<Particle>(400);
    // initiate particles in a square formation
    for (int i = 0; i < particles.getCount(); i++) {
        particles[i] = { defaultMass, defaultRadius, {100 + (float)(i % 20) * 4, 100 + (float)(i / 20) * 4}, {0, 0} };
    }

    densities = Array<float>(particles.getCount());
}

// smoothing function
float Simulation::smoothing(float radius, float dist) {
    if (dist >= radius) {
        return 0;
    }

    float volume = pow(radius, 5) * PI * 0.1;
    float difference = std::max(0.f, radius - dist);
    float value = difference * difference * difference;
    return value / volume;
}

// (r-d)^3
float Simulation::smoothingGradient(float radius, float dist) {
    if (dist >= radius) {
        return 0;
    }
    
    float difference = radius - dist;
    float gradient = -(30 * difference * difference) / (PI * pow(radius, 5));
    return gradient;
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

float Simulation::convertDensityToPressure(float density) {
    float densityError = density - targetDensity;
    float pressure = densityError * pressureMultiplier;
    return pressure;
}

float Simulation::calculateSharedPressure(float densityA, float densityB) {
    float pressureA = convertDensityToPressure(densityA);
    float pressureB = convertDensityToPressure(densityB);

    return (pressureA + pressureB) / 2;
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
        float mass = particles[i].mass;
        float density = densities[i];
        float sharedPressure = calculateSharedPressure(density, densities[i]);

        // because of how we compute gradient vector the direction automatically descends
        gradientVec = Vector2Add(gradientVec, Vector2Scale(dir, sharedPressure * magnitude * (mass / density)));
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
        particles[i].vel = Vector2Add(particles[i].vel, Vector2Scale(gravity, deltaTime * timeMultiplier));

        Vector2 gradientVec = calculateGradientVec(particles[i].pos);
        Vector2 forceVec = Vector2Divide(gradientVec, particles[i].mass);
        particles[i].vel = Vector2Add(particles[i].vel, Vector2Scale(forceVec, deltaTime * timeMultiplier));

        // boundary collisions
        if (particles[i].pos.y + particles[i].radius >= getHeight()) {
            Vector2 vel = particles[i].vel;
            particles[i].pos.y = (getHeight() - particles[i].radius);
            particles[i].vel = { vel.x, -vel.y * (1 - collisionDampening) };
        }
        else if (particles[i].pos.x - particles[i].radius <= 0) {
            Vector2 vel = particles[i].vel;
            particles[i].pos.x = particles[i].radius;
            particles[i].vel = { -vel.x * (1 - collisionDampening), vel.y };
        }
        else if (particles[i].pos.x + particles[i].radius >= getWidth()) {
            Vector2 vel = particles[i].vel;
            particles[i].pos.x = (getWidth() - particles[i].radius);
            particles[i].vel = { -vel.x * (1 - collisionDampening), vel.y };
        }
        else if (particles[i].pos.y - particles[i].radius <= 0) {
            Vector2 vel = particles[i].vel;
            particles[i].pos.y = particles[i].radius;
            particles[i].vel = { vel.x, -vel.y * (1 - collisionDampening) };
        }

        // update particle positions
        particles[i].pos = Vector2Add(particles[i].pos, Vector2Scale(particles[i].vel, deltaTime * timeMultiplier));
    }
}

void Simulation::draw() {
    if (showSmoothingRadius) {
        for (int i = 0; i < particles.getCount(); i++) {
            Vector2 pos = particles[i].pos;
            DrawCircle(bounds.x + pos.x, bounds.y + pos.y, smoothingRadius, { 0, 121, 241, 31 });
        }
    }

    for (int i = 0; i < particles.getCount(); i++) {
        Vector2 pos = particles[i].pos;
        float scaledDensity = densities[i] * 4;
        
        Color densityColor = { std::min(255.f, 255 * scaledDensity), 0, std::max(0.f, 255 * (1 - scaledDensity)), 255};
        
        DrawCircle(bounds.x + pos.x, bounds.y + pos.y, particles[i].radius, densityColor);
    }
    DrawRectangleLines(bounds.x, bounds.y, bounds.z - bounds.x, bounds.w - bounds.y, BLACK);
}
