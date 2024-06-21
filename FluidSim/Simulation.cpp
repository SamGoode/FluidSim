#include "Simulation.h"
#include "raymath.h"
#include "rlgl.h"
#include "RaylibOverloads.h"
#include <algorithm>

// must be larger than particle count
#define BUFFER_SIZE 16384
#define WORKGROUP_SIZE 512

struct SimData {
    Vector2 gravity;
    float deltaTime;
    float projectedTime;
    float timeDilation;
};

Simulation::Simulation(Vector4 _bounds) {
    bounds = _bounds;
    // scale is pixels/unit
    scale = 4;
    resolution = { getWidth(), getHeight() };
    Image whiteImage = GenImageColor(getWidth(), getHeight(), WHITE);
    texture = LoadTextureFromImage(whiteImage);
    UnloadImage(whiteImage);

    gravity = { 0, 1 };
    collisionDampening = 0.1;

    showSmoothingRadius = false;
    smoothingRadius = 1.4f;
    sqrRadius = smoothingRadius * smoothingRadius;
    targetDensity = smoothing(smoothingRadius, 0);
    pressureMultiplier = 400;
    timeDilation = 3;
    mouseInteractRadius = 10;
    mouseInteractForce = 20;

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

    particles = Array<Particle>(16384);
    // initiate particles in a square formation
    for (int i = 0; i < particles.getCount(); i++) {
        particles[i] = { defaultMass, defaultRadius, {10 + (float)(i % 128) * defaultRadius * 2, 10 + (float)(i / 128) * defaultRadius * 2}, {0, 0} };
    }

    fixedTimeStep = 0.01f;
    timePassed = 0;

    projectedPositions = Array<Vector2>(particles.getCount());
    densities = Array<float>(particles.getCount());

    spatialHash = SpatialHashGrid({getScaledWidth(), getScaledHeight()}, smoothingRadius, smoothingRadius);
    // 3x3 grid around center cell
    for (int x = 0; x < 3; x++) {
        for (int y = 0; y < 3; y++) {
            cellOffsets[x + y * 3] = { x - 1, y - 1 };
        }
    }

    char* updateParticleCode = LoadFileText("updateParticle.glsl");
    unsigned int updateParticleShader = rlCompileShader(updateParticleCode, RL_COMPUTE_SHADER);
    updateParticleProgram = rlLoadComputeShaderProgram(updateParticleShader);
    UnloadFileText(updateParticleCode);

    char* gravProjectionCode = LoadFileText("gravProjection.glsl");
    unsigned int gravProjectionShader = rlCompileShader(gravProjectionCode, RL_COMPUTE_SHADER);
    gravProjectionProgram = rlLoadComputeShaderProgram(gravProjectionShader);
    UnloadFileText(gravProjectionCode);

    renderSimShader = LoadShader(NULL, "renderSim.glsl");
    resUniformLoc = GetShaderLocation(renderSimShader, "resolution");

    particleSSBO = rlLoadShaderBuffer(particles.getCount() * sizeof(Particle), NULL, RL_DYNAMIC_COPY);
    projectedPositionSSBO = rlLoadShaderBuffer(projectedPositions.getCount() * sizeof(Vector2), NULL, RL_DYNAMIC_COPY);
    simDataSSBO = rlLoadShaderBuffer(sizeof(SimData), NULL, RL_DYNAMIC_COPY);
    textureBuffer = Array<Vector4>(getWidth() * getHeight());
    textureBuffer.clear({1, 1, 1, 0});
    textureSSBO = rlLoadShaderBuffer(getWidth() * getHeight() * sizeof(Vector4), NULL, RL_DYNAMIC_COPY);

    rlUpdateShaderBuffer(particleSSBO, particles.begin(), particles.getCount() * sizeof(Particle), 0);
}

Simulation::~Simulation() {
    rlUnloadShaderBuffer(particleSSBO);
    rlUnloadShaderBuffer(projectedPositionSSBO);
    rlUnloadShaderBuffer(simDataSSBO);
    rlUnloadShaderBuffer(textureSSBO);

    rlUnloadShaderProgram(updateParticleProgram);
    rlUnloadShaderProgram(gravProjectionProgram);

    UnloadTexture(texture);
    UnloadShader(renderSimShader);
}

// smoothing function
// dist should never exceed radius
float Simulation::smoothing(float radius, float dist) {
    //float volume = pow(radius, 5) * PI * 0.1;
    float difference = radius - dist;//std::max(0.f, radius - dist);
    float value = (difference * difference) / (radius * radius);
    return value;
}

// derivative of smoothing function
// dist should never exceed radius
float Simulation::smoothingGradient(float radius, float dist) {
    float gradient = (2 * (dist - radius)) / (radius * radius);
    return gradient;
}

// calculates density level at a particle's position
float Simulation::calculateDensity(int particleID) {
    const Array<int2>& hashList = spatialHash.getHashList();
    const Array<int2>& indexLookup = spatialHash.getIndexLookup();

    Vector2 pos = projectedPositions[particleID];
    int2 cellPos = spatialHash.getCellPos(pos);

    float density = 0;
    for (int i = 0; i < 9; i++) {
        int2 offsetCellPos = cellPos + cellOffsets[i];

        if (!spatialHash.isValidCellPos(offsetCellPos)) {
            continue;
        }

        int cellHash = spatialHash.getCellHash(offsetCellPos);

        int startIndex = indexLookup[cellHash].x;
        if (startIndex < 0) {
            continue;
        }

        int endIndex = indexLookup[cellHash].y;

        for (int n = startIndex; n < endIndex + 1; n++) {
            int otherParticleID = hashList[n].x;
            Vector2 otherPos = projectedPositions[otherParticleID];

            float sqrDist = Vector2DistanceSqr(pos, otherPos);

            if (sqrDist > sqrRadius) {
                continue;
            }

            float dist = sqrt(sqrDist);
            density += smoothing(smoothingRadius, dist) * particles[otherParticleID].mass;
        }
    }

    return density;
}

float Simulation::convertDensityToPressure(float density) {
    float densityError = std::max(0.f, density - targetDensity);
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
    const Array<int2>& hashList = spatialHash.getHashList();
    const Array<int2>& indexLookup = spatialHash.getIndexLookup();

    float density = densities[particleID];
    Vector2 pos = projectedPositions[particleID];
    int2 cellPos = spatialHash.getCellPos(pos);

    Vector2 gradientVec = { 0, 0 };
    for (int i = 0; i < 9; i++) {
        int2 offsetCellPos = cellPos + cellOffsets[i];

        if (!spatialHash.isValidCellPos(offsetCellPos)) {
            continue;
        }

        int cellHash = spatialHash.getCellHash(offsetCellPos);

        int startIndex = indexLookup[cellHash].x;
        if (startIndex < 0) {
            continue;
        }

        int endIndex = indexLookup[cellHash].y;

        for (int n = startIndex; n < endIndex + 1; n++) {
            int otherParticleID = hashList[n].x;
            if (particleID == otherParticleID) {
                continue;
            }

            Vector2 otherPos = projectedPositions[otherParticleID];
            Vector2 posToOtherPos = otherPos - pos;

            float sqrDist = Vector2LengthSqr(posToOtherPos);
            if (sqrDist > sqrRadius) {
                continue;
            }

            float dist = sqrt(sqrDist);
            Vector2 dir;

            if (dist == 0) {
                // generate random direction vector if other particle is right on top of this one
                float randAngle = (rand() * 2 * PI) / RAND_MAX;
                dir = { (float)cos(randAngle), (float)sin(randAngle) };
            }
            else {
                // unit direction vector
                dir = posToOtherPos / dist;
            }

            // magnitude of gradient vector
            float magnitude = smoothingGradient(smoothingRadius, dist);
            float otherMass = particles[otherParticleID].mass;
            float otherDensity = densities[otherParticleID];

            float sharedPressure = calculateSharedPressure(density, otherDensity);

            // because of how we compute gradient vector the direction automatically descends
            gradientVec += (dir * magnitude) * (sharedPressure * (otherMass / density));
        }
    }

    return gradientVec;
}

void Simulation::update(float deltaTime) {
    // Cap in case of step number being too high causing a low performance loop
    if (deltaTime > 0.02f) {
        deltaTime = 0.02f;
    }

    timePassed += deltaTime;

    int numberSteps = 0;
    if (timePassed >= fixedTimeStep) {
        numberSteps = (int)floor(timePassed / fixedTimeStep);
        timePassed = timePassed - (fixedTimeStep * numberSteps);
    }

    for (int step = 0; step < numberSteps; step++) {
        //SimData simData;
        //simData.gravity = gravity;
        //simData.deltaTime = deltaTime;
        //simData.projectedTime = projectedTime;
        //simData.timeDilation = timeDilation;

        //rlUpdateShaderBuffer(simDataSSBO, &simData, sizeof(SimData), 0);

        // gravity application
        for (int i = 0; i < particles.getCount(); i++) {
            particles[i].vel += gravity * (fixedTimeStep * timeDilation);
        }

        // update projected positions cache
        for (int i = 0; i < particles.getCount(); i++) {
            projectedPositions[i] = particles[i].pos + (particles[i].vel * (fixedTimeStep * timeDilation));
        }

        //// loading particle data into shader buffer
        //rlUpdateShaderBuffer(particleSSBO, particles.begin(), particles.getCount() * sizeof(Particle), 0);

        //// processing particle data using compute shader
        //rlEnableShader(gravProjectionProgram);
        //rlBindShaderBuffer(particleSSBO, 1);
        //rlBindShaderBuffer(simDataSSBO, 2);
        //rlBindShaderBuffer(projectedPositionSSBO, 3);
        //rlComputeShaderDispatch(BUFFER_SIZE / WORKGROUP_SIZE, 1, 1);
        //rlDisableShader();

        //// unloading data from shader buffer
        //rlReadShaderBuffer(particleSSBO, particles.begin(), particles.getCount() * sizeof(Particle), 0);
        //rlReadShaderBuffer(projectedPositionSSBO, projectedPositions.begin(), projectedPositions.getCount() * sizeof(Vector2), 0);

        spatialHash.generateHashList(projectedPositions);
        spatialHash.sortByCellHash();
        spatialHash.generateLookup();

        // update densities cache
        for (int particleID = 0; particleID < particles.getCount(); particleID++) {
            densities[particleID] = calculateDensity(particleID);
        }

        // calculate cumulative pressure force vector and apply it
        for (int particleID = 0; particleID < particles.getCount(); particleID++) {
            float mass = particles[particleID].mass;

            Vector2 gradientVec = calculateGradientVec(particleID);

            Vector2 forceVec = gradientVec / mass;
            particles[particleID].vel += forceVec * (fixedTimeStep * timeDilation);
        }

        // mouse interaction
        if (IsMouseButtonDown(MOUSE_LEFT_BUTTON)) {
            for (int i = 0; i < particles.getCount(); i++) {
                Vector2 mousePos = convertToSimPos(GetMousePosition());
                Vector2 MtoP = particles[i].pos - mousePos;
                float dist = Vector2Length(MtoP);
                if (dist <= mouseInteractRadius) {
                    Vector2 unitDir = MtoP / dist;
                    Vector2 forceVec = unitDir * -mouseInteractForce;
                    particles[i].vel += forceVec * (fixedTimeStep * timeDilation);
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
                    particles[i].vel += forceVec * (fixedTimeStep * timeDilation);
                }
            }
        }

        // boundary collision check
        for (int i = 0; i < particles.getCount(); i++) {
            float edgeOffset = particles[i].radius * 1.2f;

            if (particles[i].pos.y + particles[i].radius >= getScaledHeight()) {
                Vector2 vel = particles[i].vel;
                particles[i].pos.y = (getScaledHeight() - edgeOffset);
                particles[i].vel.y = -abs(vel.y) * (1 - collisionDampening);
            }
            else if (particles[i].pos.y - particles[i].radius <= 0) {
                Vector2 vel = particles[i].vel;
                particles[i].pos.y = edgeOffset;
                particles[i].vel.y = abs(vel.y) * (1 - collisionDampening);
            }

            if (particles[i].pos.x - particles[i].radius <= 0) {
                Vector2 vel = particles[i].vel;
                particles[i].pos.x = edgeOffset;
                particles[i].vel.x = abs(vel.x) * (1 - collisionDampening);
            }
            else if (particles[i].pos.x + particles[i].radius >= getScaledWidth()) {
                Vector2 vel = particles[i].vel;
                particles[i].pos.x = (getScaledWidth() - edgeOffset);
                particles[i].vel.x = -abs(vel.x) * (1 - collisionDampening);
            }
        }

        // update particle positions
        for (int i = 0; i < particles.getCount(); i++) {
            particles[i].pos += particles[i].vel * (fixedTimeStep * timeDilation);
        }

        //// loading particle data into shader buffer
        //rlUpdateShaderBuffer(particleSSBO, particles.begin(), particles.getCount() * sizeof(Particle), 0);
        //rlUpdateShaderBuffer(textureSSBO, textureBuffer.begin(), textureBuffer.getCount() * sizeof(Vector4), 0);

        //// processing particle data using compute shader
        //rlEnableShader(updateParticleProgram);
        //rlBindShaderBuffer(particleSSBO, 1);
        //rlBindShaderBuffer(simDataSSBO, 2);
        //rlBindShaderBuffer(textureSSBO, 3);
        //rlComputeShaderDispatch(BUFFER_SIZE/WORKGROUP_SIZE, 1, 1);
        //rlDisableShader();

        //rlReadShaderBuffer(particleSSBO, particles.begin(), particles.getCount() * sizeof(Particle), 0);
    }

    // Iterated through all fixed updates
}

void Simulation::draw() {
    //rlBindShaderBuffer(textureSSBO, 1);
    //SetShaderValue(renderSimShader, resUniformLoc, &resolution, SHADER_UNIFORM_VEC2);

    //BeginShaderMode(renderSimShader);
    //DrawTexture(texture, bounds.x, bounds.y, WHITE);
    //EndShaderMode();

    Vector2 mousePos = GetMousePosition();

    if (showSmoothingRadius) {
        for (int i = 0; i < particles.getCount(); i++) {
            Vector2 screenPos = convertToScreenPos(particles[i].pos);
            DrawCircle(screenPos.x, screenPos.y, smoothingRadius * scale, { 0, 121, 241, 31 });
        }
    }

    for (int i = 0; i < particles.getCount(); i++) {
        Vector2 screenPos = convertToScreenPos(particles[i].pos);
        float densityError = std::clamp((densities[i] / targetDensity) / 2.5f, 0.f, 1.f);
        
        float r = (1 - abs(densityError - 1)) * 255;
        float g = (1 - abs(densityError - 0.5) * 2) * 255;
        float b = (1 - abs(densityError - 0)) * 255;

        //float densityError = std::max((densities[i] - targetDensity), 0.f) * 255 * 2;

        //float r = std::max(0.f, std::min(255.f, densityError));
        //float g = std::max(0.f, std::min(255.f, densityError - 255));
        //float b = std::max(0.f, std::min(255.f, densityError - 511));

        Color densityColor = { r, g, b, 255 };
        
        //DrawCircle(screenPos.x, screenPos.y, particles[i].radius * scale, densityColor);
        DrawRectangle(screenPos.x - (particles[i].radius * scale), screenPos.y - (particles[i].radius * scale), 2 * particles[i].radius * scale, 2 * particles[i].radius * scale, densityColor);
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
