#include "Simulation.h"
#include "raymath.h"
#include "RaylibOverloads.h"
#include <algorithm>

Simulation::Simulation(Vector4 _bounds) {
    bounds = _bounds;
    // scale is pixels/unit
    scale = 8;

    gravity = { 0, 1 };
    collisionDampening = 0.1;

    showSmoothingRadius = false;
    smoothingRadius = 2;
    targetDensity = smoothing(smoothingRadius, 0);
    pressureMultiplier = 100;
    timeMultiplier = 4;
    mouseInteractRadius = 5;
    mouseInteractForce = 10;

    defaultMass = 1;
    defaultRadius = 0.5;
    //particles = Array<Particle>(1024, 1024);
    //for (int i = 0; i < particles.getCount(); i++) {
    //    // generate random location
    //    float x = (rand() * (bounds.z - bounds.x)/scale) / RAND_MAX;
    //    float y = (rand() * (bounds.w - bounds.y)/scale) / RAND_MAX;
    //    Vector2 position = { x, y };


    //    // generate random velocity
    //    float velX = ((rand() * 20) / RAND_MAX) - 10;
    //    float velY = ((rand() * 20) / RAND_MAX) - 10;
    //    Vector2 vel = { velX, velY };
    //    //Vector2 vel = { 0, 0 };

    //    particles[i] = { defaultMass, defaultRadius, position, vel };
    //}

    particles = Array<Particle>(1024);
    // initiate particles in a square formation
    for (int i = 0; i < particles.getCount(); i++) {
        particles[i] = { defaultMass, defaultRadius, {10 + (float)(i % 32) * defaultRadius * 2, 10 + (float)(i / 32) * defaultRadius * 2}, {0, 0} };
    }

    projectedPositions = Array<Vector2>(particles.getCount());
    densities = Array<float>(particles.getCount());

    spatialHash = SpatialHashGrid({getScaledWidth(), getScaledHeight()}, smoothingRadius, smoothingRadius);
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

    Array<int> nearbyIDs = spatialHash.findNearby(spatialHash.getCellPos(pos));
    for (int i = 0; i < nearbyIDs.getCount(); i++) {
        float dist = Vector2Distance(pos, projectedPositions[nearbyIDs[i]]);
        density += smoothing(smoothingRadius, dist) * particles[nearbyIDs[i]].mass;
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
Vector2 Simulation::calculateGradientVec(int particleID) {
    Vector2 gradientVec = { 0, 0 };
    float density = densities[particleID];
    Vector2 pos = projectedPositions[particleID];

    Array<int> nearbyIDs = spatialHash.findNearby(spatialHash.getCellPos(pos));
    for (int i = 0; i < nearbyIDs.getCount(); i++) {
        if (particleID == nearbyIDs[i]) {
            continue;
        }

        Vector2 pointToParticle = projectedPositions[nearbyIDs[i]] - pos;
        float dist = Vector2Length(pointToParticle);

        Vector2 dir;
        // generate random direction vector if other particle is right on top of this one
        if (dist == 0) { 
            float randAngle = (rand() * 2 * PI) / RAND_MAX;
            dir = { (float)cos(randAngle), (float)sin(randAngle) };
        }
        else {
            // normalised direction vector
            dir = Vector2Divide(pointToParticle, dist);
        }

        // magnitude of gradient vector
        float magnitude = smoothingGradient(smoothingRadius, dist);
        float mass = particles[nearbyIDs[i]].mass;
        
        float sharedPressure = calculateSharedPressure(density, densities[nearbyIDs[i]]);

        // because of how we compute gradient vector the direction automatically descends
        gradientVec += (dir * magnitude) * (sharedPressure * (mass / density));
    }
    return gradientVec;
}

//force = change in mass * velocity / change in time
void Simulation::update(float deltaTime) {
    if (deltaTime > 0.01) {
        deltaTime = 0.01;
    }

    // gravity application
    for (int i = 0; i < particles.getCount(); i++) {
        particles[i].vel += gravity * (deltaTime * timeMultiplier);
    }

    // update projected positions cache
    for (int i = 0; i < particles.getCount(); i++) {
        projectedPositions[i] = particles[i].pos + (particles[i].vel * ((1.f/120.f) * timeMultiplier));
    }

    spatialHash.generateHashList(projectedPositions);
    spatialHash.sortByCellHash();
    spatialHash.generateLookup();

    // update densities cache
    for (int i = 0; i < particles.getCount(); i++) {
        densities[i] = calculateDensity(projectedPositions[i]);
    }

    // apply pressure vectors to particles
    for (int particleID = 0; particleID < particles.getCount(); particleID++) {
        Vector2 gradientVec = calculateGradientVec(particleID);
        Vector2 forceVec = gradientVec / particles[particleID].mass;
        particles[particleID].vel += forceVec * (deltaTime * timeMultiplier);
    }

    if (IsMouseButtonDown(MOUSE_LEFT_BUTTON)) {
        for (int i = 0; i < particles.getCount(); i++) {
            Vector2 mousePos = convertToSimPos(GetMousePosition());
            Vector2 MtoP = particles[i].pos - mousePos;
            float dist = Vector2Length(MtoP);
            if (dist <= mouseInteractRadius) {
                Vector2 unitDir = MtoP / dist;
                Vector2 forceVec = unitDir * -mouseInteractForce;
                particles[i].vel += forceVec * (deltaTime * timeMultiplier);
            }
        }
    }
    else if (IsMouseButtonDown(MOUSE_RIGHT_BUTTON)) {
        for (int i = 0; i < particles.getCount(); i++) {
            Vector2 mousePos = convertToSimPos(GetMousePosition());
            Vector2 MtoP = particles[i].pos - mousePos;
            float dist = Vector2Length(MtoP);
            if (dist <= mouseInteractRadius) {
                Vector2 unitDir = MtoP / dist;
                Vector2 forceVec = unitDir * mouseInteractForce;
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
    Vector2 mousePos = GetMousePosition();

    if (showSmoothingRadius) {
        for (int i = 0; i < particles.getCount(); i++) {
            Vector2 screenPos = convertToScreenPos(particles[i].pos);
            DrawCircle(screenPos.x, screenPos.y, smoothingRadius * scale, { 0, 121, 241, 31 });
        }
    }

    for (int i = 0; i < particles.getCount(); i++) {
        Vector2 screenPos = convertToScreenPos(particles[i].pos);
        float densityError = std::clamp((densities[i] / targetDensity) / 2.8f, 0.f, 1.f);
        
        float r = (1 - abs(densityError - 1)) * 255;
        float g = (1 - abs(densityError - 0.5) * 2) * 0.7f * 255;
        float b = (1 - abs(densityError - 0)) * 255;

        Color densityColor = { r, g, b, 255};
        
        DrawCircle(screenPos.x, screenPos.y, particles[i].radius * scale, densityColor);
        //DrawRectangle(screenPos.x - (particles[i].radius * scale), screenPos.y - (particles[i].radius * scale), 2 * particles[i].radius * scale, 2 * particles[i].radius * scale, densityColor);
    }

    // spatial hash testing
    //spatialHash.draw({ bounds.x, bounds.y }, scale, convertToSimPos(mousePos));
    //Array<int> nearbyIDs = spatialHash.findNearby(spatialHash.getCellPos(convertToSimPos(mousePos)));
    //Array<Vector2> nearbyPositions(nearbyIDs.getCount());
    //for (int i = 0; i < nearbyIDs.getCount(); i++) {
    //    Vector2 screenPos = convertToScreenPos(projectedPositions[nearbyIDs[i]]);
    //    DrawCircle(screenPos.x, screenPos.y, defaultRadius * scale * 2, WHITE);
    //}

    DrawCircleLines(mousePos.x, mousePos.y, mouseInteractRadius * scale, BLACK);

    DrawRectangleLines(bounds.x, bounds.y, bounds.z - bounds.x, bounds.w - bounds.y, BLACK);
}
