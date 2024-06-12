#include "Simulation.h"
#include "raymath.h"
#include "RaylibOverloads.h"

Simulation::Simulation(Vector4 _bounds) {
    bounds = _bounds;
    // scale is pixels/unit
    scale = 8;

    gravity = { 0, 1 };
    collisionDampening = 0.2;

    showSmoothingRadius = false;
    smoothingRadius = 1.2;
    targetDensity = smoothing(smoothingRadius, 0);
    pressureMultiplier = 100;
    timeMultiplier = 4;

    defaultMass = 1;
    defaultRadius = 0.5;
    //particles = Array<Particle>(400);
    //for (int i = 0; i < particles.getCount(); i++) {
    //    // generate random location
    //    float x = (rand() * (bounds.z - bounds.x)/scale) / RAND_MAX;
    //    float y = (rand() * (bounds.w - bounds.y)/scale) / RAND_MAX;
    //    Vector2 position = { x, y };


    //    // generate random velocity
    //    float velX = ((rand() * 40) / RAND_MAX) - 20;
    //    float velY = ((rand() * 40) / RAND_MAX) - 20;
    //    Vector2 vel = { velX, velY };
    //    //Vector2 vel = { 0, 0 };

    //    particles[i] = { defaultMass, defaultRadius, position, vel };
    //}

    particles = Array<Particle>(900);
    // initiate particles in a square formation
    for (int i = 0; i < particles.getCount(); i++) {
        particles[i] = { defaultMass, defaultRadius, {5 + (float)(i % 30) * 1, 5 + (float)(i / 30) * 1}, {0, 0} };
    }

    projectedPositions = Array<Vector2>(particles.getCount());
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

    for (int i = 0; i < projectedPositions.getCount(); i++) {
        float dist = Vector2Distance(pos, projectedPositions[i]);
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
// then scales it based on mass, density and pressure to get the pressure vector acting upon a particle
Vector2 Simulation::calculateGradientVec(Vector2 pos) {
    Vector2 gradientVec = { 0, 0 };
    float density = calculateDensity(pos);
    for (int i = 0; i < particles.getCount(); i++) {
        Vector2 pointToParticle = projectedPositions[i] - pos;
        float dist = Vector2Length(pointToParticle);

        // skip if particle is right on top of location
        if (dist == 0 || dist >= smoothingRadius) { continue; }

        // normalised direction vector
        Vector2 dir = Vector2Divide(pointToParticle, dist);
        // magnitude of gradient vector
        float magnitude = smoothingGradient(smoothingRadius, dist);
        float mass = particles[i].mass;
        
        float sharedPressure = calculateSharedPressure(density, densities[i]);

        // because of how we compute gradient vector the direction automatically descends
        gradientVec += (dir * magnitude) * (sharedPressure * (mass / density));
    }
    return gradientVec;
}

//force = change in mass * velocity / change in time
void Simulation::update(float deltaTime) {
    if (deltaTime > 0.5) {
        return;
    }

    // gravity application
    for (int i = 0; i < particles.getCount(); i++) {
        particles[i].vel += gravity * (deltaTime * timeMultiplier);
    }

    // update projected positions cache
    for (int i = 0; i < particles.getCount(); i++) {
        projectedPositions[i] = particles[i].pos + (particles[i].vel * (deltaTime * timeMultiplier));
    }

    // update densities cache
    for (int i = 0; i < particles.getCount(); i++) {
        densities[i] = calculateDensity(projectedPositions[i]);
    }

    // apply pressure vectors to particles
    for (int i = 0; i < particles.getCount(); i++) {
        Vector2 gradientVec = calculateGradientVec(projectedPositions[i]);
        Vector2 forceVec = gradientVec / particles[i].mass;
        particles[i].vel += forceVec * (deltaTime * timeMultiplier);
    }

    if (IsMouseButtonDown(MOUSE_LEFT_BUTTON)) {
        for (int i = 0; i < particles.getCount(); i++) {
            Vector2 mousePos = convertToSimPos(GetMousePosition());
            Vector2 MtoP = particles[i].pos - mousePos;
            float dist = Vector2Length(MtoP);
            if (dist <= 5) {
                Vector2 unitDir = MtoP / dist;
                Vector2 forceVec = unitDir * -4;
                particles[i].vel += forceVec * (deltaTime * timeMultiplier);
            }
        }
    }
    else if (IsMouseButtonDown(MOUSE_RIGHT_BUTTON)) {
        for (int i = 0; i < particles.getCount(); i++) {
            Vector2 mousePos = convertToSimPos(GetMousePosition());
            Vector2 MtoP = particles[i].pos - mousePos;
            float dist = Vector2Length(MtoP);
            if (dist <= 5) {
                Vector2 unitDir = MtoP / dist;
                Vector2 forceVec = unitDir * 4;
                particles[i].vel += forceVec * (deltaTime * timeMultiplier);
            }
        }
    }

    for (int i = 0; i < particles.getCount(); i++) {
        

        // boundary collisions
        if (particles[i].pos.y + particles[i].radius >= getScaledHeight()) {
            Vector2 vel = particles[i].vel;
            particles[i].pos.y = (getScaledHeight() - particles[i].radius);
            particles[i].vel.y = -vel.y * (1 - collisionDampening);
        }
        else if (particles[i].pos.y - particles[i].radius <= 0) {
            Vector2 vel = particles[i].vel;
            particles[i].pos.y = particles[i].radius;
            particles[i].vel.y = -vel.y * (1 - collisionDampening);
        }

        if (particles[i].pos.x - particles[i].radius <= 0) {
            Vector2 vel = particles[i].vel;
            particles[i].pos.x = particles[i].radius;
            particles[i].vel.x = -vel.x * (1 - collisionDampening);
        }
        else if (particles[i].pos.x + particles[i].radius >= getScaledWidth()) {
            Vector2 vel = particles[i].vel;
            particles[i].pos.x = (getScaledWidth() - particles[i].radius);
            particles[i].vel.x = -vel.x * (1 - collisionDampening);
        }
        
        // update particle positions
        particles[i].pos += particles[i].vel * (deltaTime * timeMultiplier);
    }
}

void Simulation::draw() {
    if (showSmoothingRadius) {
        for (int i = 0; i < particles.getCount(); i++) {
            Vector2 screenPos = convertToScreenPos(particles[i].pos);
            DrawCircle(screenPos.x, screenPos.y, smoothingRadius * scale, { 0, 121, 241, 31 });
        }
    }

    for (int i = 0; i < particles.getCount(); i++) {
        Vector2 screenPos = convertToScreenPos(particles[i].pos);
        float densityError = std::max(0.f, (densities[i] - targetDensity)) * 255 * 2 * 1.5;

        float r = std::max(0.f, std::min(255.f, densityError));
        float g = std::max(0.f, std::min(255.f, densityError - 255));
        float b = std::max(0.f, std::min(255.f, densityError - 511));

        Color densityColor = { r, g, b, 255};
        
        DrawCircle(screenPos.x, screenPos.y, particles[i].radius * scale, densityColor);
    }

    Vector2 mousePos = GetMousePosition();
    DrawCircleLines(mousePos.x, mousePos.y, 5 * scale, BLACK);

    DrawRectangleLines(bounds.x, bounds.y, bounds.z - bounds.x, bounds.w - bounds.y, BLACK);
}
