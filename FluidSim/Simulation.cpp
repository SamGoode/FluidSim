#include "Simulation.h"
#include "raymath.h"
#include "rlgl.h"
#include "RaylibOverloads.h"
#include <algorithm>

// must be larger than particle count
#define MAX_PARTICLE_COUNT 8192
#define WORKGROUP_SIZE 512

Simulation::Simulation(Vector4 _bounds) {
    bounds = _bounds;
    // scale is pixels/unit
    scale = 6;
    resolution = { getWidth(), getHeight() };
    Image whiteImage = GenImageColor(getWidth(), getHeight(), WHITE);
    texture = LoadTextureFromImage(whiteImage);
    UnloadImage(whiteImage);

    gravity = { 0, 2 };
    collisionDampening = 0.1;

    showSmoothingRadius = false;
    smoothingRadius = 1.3f;
    sqrRadius = smoothingRadius * smoothingRadius;
    particleRadius = 0.5;
    targetDensity = smoothing(smoothingRadius, 0);
    pressureMultiplier = 400;
    timeDilation = 2;

    mouseInteractRadius = 8;
    mouseInteractForce = 50;

    airflows = Array<AirflowSpace>(0);
    //airflows[0] = {
    //    {0, 30},
    //    100,
    //    40,
    //    {2, 0}
    //};
    //airflows[1] = {
    //    {30, 0},
    //    103.33f,
    //    20,
    //    {-6, 0}
    //};
    //airflows[2] = {
    //    {30, 80},
    //    103.33f,
    //    20,
    //    {-6, 0}
    //};
    //airflows[3] = {
    //    {100, 20},
    //    33.33f,
    //    30,
    //    {0, -4}
    //};
    //airflows[4] = {
    //    {100, 50},
    //    33.33f,
    //    30,
    //    {0, 4}
    //};
    //airflows[5] = {
    //    {0, 0},
    //    30,
    //    30,
    //    {0, 4}
    //};
    //airflows[6] = {
    //    {0, 70},
    //    30,
    //    30,
    //    {0, -4}
    //};

    ball = {
        {80, 50},
        5
    };


    defaultMass = 1;
    particles = Array<Particle>(MAX_PARTICLE_COUNT);
    objectPool = Array<int>(MAX_PARTICLE_COUNT);
    for (int i = 0; i < MAX_PARTICLE_COUNT; i++) {
        objectPool[i] = i;
    }
    activeCount = 8192;
    
    // initiate particles at random positions and with random velocities
    for (int i = 0; i < activeCount; i++) {
        // generate random location
        float x = (rand() * (bounds.z - bounds.x)/scale) / RAND_MAX;
        float y = (rand() * (bounds.w - bounds.y)/scale) / RAND_MAX;
        Vector2 position = { x, y };

        // generate random velocity
        float velX = ((rand() * 10) / RAND_MAX) - 5;
        float velY = ((rand() * 10) / RAND_MAX) - 5;
        //Vector2 vel = { velX, velY };
        Vector2 vel = { 0, 0 };

        particles[i] = { 1, defaultMass, position, vel };
    }
    
    // initiate particles in a square formation
    //for (int i = 0; i < particles.getCount(); i++) {
    //    particles[i] = { defaultMass, defaultRadius, {10 + (float)(i % 128) * defaultRadius * 1.2f, 5 + (float)(i / 128) * defaultRadius * 1.2f}, {0, 0} };
    //}

    fixedTimeStep = 0.02f;
    timePassed = 0;

    densities = Array<float>(particles.getCount());
    previousPositions = Array<Vector2>(particles.getCount());
    for (int i = 0; i < particles.getCount(); i++) {
        previousPositions[i] = particles[i].pos;
    }

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

    char* clearTextureBufferCode = LoadFileText("clearTextureBuffer.glsl");
    unsigned int clearTextureBufferShader = rlCompileShader(clearTextureBufferCode, RL_COMPUTE_SHADER);
    clearTextureBufferProgram = rlLoadComputeShaderProgram(clearTextureBufferShader);
    UnloadFileText(clearTextureBufferCode);

    char* updateTextureBufferCode = LoadFileText("updateTextureBuffer.glsl");
    unsigned int updateTextureBufferShader = rlCompileShader(updateTextureBufferCode, RL_COMPUTE_SHADER);
    updateTextureBufferProgram = rlLoadComputeShaderProgram(updateTextureBufferShader);
    UnloadFileText(updateTextureBufferCode);

    renderSimShader = LoadShader(NULL, "renderSim.glsl");
    resUniformLoc = GetShaderLocation(renderSimShader, "resolution");

    simDataSSBO = rlLoadShaderBuffer(sizeof(SimData), NULL, RL_DYNAMIC_COPY);
    particleSSBO = rlLoadShaderBuffer(particles.getCount() * sizeof(Particle), NULL, RL_DYNAMIC_COPY);
    densitySSBO = rlLoadShaderBuffer(densities.getCount() * sizeof(float), NULL, RL_DYNAMIC_COPY);
    textureSSBO = rlLoadShaderBuffer(getWidth() * getHeight() * sizeof(Vector4), NULL, RL_DYNAMIC_COPY);

    simData = {
        gravity,
        targetDensity,
        fixedTimeStep
    };
    rlUpdateShaderBuffer(simDataSSBO, &simData, sizeof(simData), 0);
}

Simulation::~Simulation() {
    rlUnloadShaderBuffer(simDataSSBO);
    rlUnloadShaderBuffer(particleSSBO);
    rlUnloadShaderBuffer(projectedPositionSSBO);
    rlUnloadShaderBuffer(densitySSBO);
    rlUnloadShaderBuffer(textureSSBO);

    rlUnloadShaderProgram(updateParticleProgram);
    rlUnloadShaderProgram(gravProjectionProgram);
    rlUnloadShaderProgram(clearTextureBufferProgram);
    rlUnloadShaderProgram(updateTextureBufferProgram);

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

    //Vector2 pos = projectedPositions[particleID];
    Vector2 pos = particles[particleID].pos;
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
            //Vector2 otherPos = projectedPositions[otherParticleID];
            Vector2 otherPos = particles[otherParticleID].pos;

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
    //const Array<int>& hashOffsets = spatialHash.getHashOffsets();

    const Array<int2>& hashList = spatialHash.getHashList();
    const Array<int2>& indexLookup = spatialHash.getIndexLookup();

    float density = densities[particleID];
    Vector2 pos = particles[particleID].pos;
    int2 cellPos = spatialHash.getCellPos(pos);
    //int cellHash = spatialHash.getCellHash(cellPos);

    Vector2 gradientVec = { 0, 0 };
    for (int i = 0; i < 9; i++) {
        int2 offsetCellPos = cellPos + cellOffsets[i];

        if (!spatialHash.isValidCellPos(offsetCellPos)) {
            continue;
        }

        int offsetCellHash = spatialHash.getCellHash(offsetCellPos);
        //int offsetCellHash = cellHash + hashOffsets[i];
        //if (!spatialHash.isValidCellHash(offsetCellHash)) {
        //    continue;
        //}

        int startIndex = indexLookup[offsetCellHash].x;
        if (startIndex < 0) {
            continue;
        }

        int endIndex = indexLookup[offsetCellHash].y;

        for (int n = startIndex; n < endIndex + 1; n++) {
            int otherParticleID = hashList[n].x;
            if (particleID == otherParticleID) {
                continue;
            }

            //Vector2 otherPos = projectedPositions[otherParticleID];
            Vector2 otherPos = particles[otherParticleID].pos;
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

void Simulation::spawnParticle(float mass, Vector2 pos, Vector2 vel) {
    if (activeCount >= objectPool.getCount()) {
        return;
    }

    previousPositions[objectPool[activeCount]] = pos;
    particles[objectPool[activeCount]] = { 1, mass, pos, vel };
    activeCount++;
}

void Simulation::despawnParticle(int poolIndex) {
    int particleID = objectPool[poolIndex];
    particles[particleID].isActive = 0;

    objectPool[poolIndex] = objectPool[activeCount - 1];
    objectPool[activeCount - 1] = particleID;
    activeCount--;
}

void Simulation::update(float deltaTime) {
    timePassed += deltaTime * timeDilation;

    int numberSteps = 0;
    if (timePassed >= fixedTimeStep) {
        numberSteps = (int)floor(timePassed / fixedTimeStep);
        timePassed = timePassed - (fixedTimeStep * numberSteps);
    }

    // Cap in case of step number being too high causing a low performance loop
    if (numberSteps > 4) {
        numberSteps = 4;
    }

    for (int step = 0; step < numberSteps; step++) {
        // gravity application
        for (int poolIndex = 0; poolIndex < activeCount; poolIndex++) {
            int particleID = objectPool[poolIndex];
            particles[particleID].vel += gravity * fixedTimeStep;
        }

        // save current positions and shift particle positions to projected positions
        for (int poolIndex = 0; poolIndex < activeCount; poolIndex++) {
            int particleID = objectPool[poolIndex];
            previousPositions[particleID] = particles[particleID].pos;
            particles[particleID].pos += particles[particleID].vel * fixedTimeStep;
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

        spatialHash.generateHashList(particles, objectPool, activeCount);
        spatialHash.sortByCellHash();
        spatialHash.generateLookup();

        // update densities cache
        for (int poolIndex = 0; poolIndex < activeCount; poolIndex++) {
            int particleID = objectPool[poolIndex];
            densities[particleID] = calculateDensity(particleID);
        }

        // calculate cumulative pressure force vector and apply it
        for (int poolIndex = 0; poolIndex < activeCount; poolIndex++) {
            int particleID = objectPool[poolIndex];

            float mass = particles[particleID].mass;
            Vector2 pressureForce = calculateGradientVec(particleID);
            particles[particleID].pos += pressureForce * ((fixedTimeStep * fixedTimeStep) / mass);
        }

        // airflow spaces
        for (int i = 0; i < airflows.getCount(); i++) {
            for (int poolIndex = 0; poolIndex < activeCount; poolIndex++) {
                int particleID = objectPool[poolIndex];

                Vector2 particlePos = particles[particleID].pos;
                float mass = particles[particleID].mass;

                if (particlePos.x < airflows[i].pos.x || particlePos.x > airflows[i].pos.x + airflows[i].width || particlePos.y < airflows[i].pos.y || particlePos.y > airflows[i].pos.y + airflows[i].height) {
                    continue;
                }

                particles[particleID].pos += airflows[i].force * (fixedTimeStep * fixedTimeStep / mass);
            }
        }

        // ball collisions
        for (int poolIndex = 0; poolIndex < activeCount; poolIndex++) {
            int particleID = objectPool[poolIndex];
            
            Vector2 BtoP = particles[particleID].pos - ball.pos;
            float sqrDist = Vector2LengthSqr(BtoP);
            if (sqrDist >= ball.radius * ball.radius) {
                continue;
            }

            float dist = sqrt(sqrDist);
            Vector2 unitNormal = BtoP / dist;

            Vector2 velocity = (particles[particleID].pos - previousPositions[particleID]) / fixedTimeStep;
            float velNormalMag = Vector2DotProduct(velocity, unitNormal);
            Vector2 velNormal = unitNormal * velNormalMag;
            Vector2 velTangent = velocity - velNormal;

            float frictionCoefficient = 0.8f;
            Vector2 impulse = (velNormal - (velTangent * (1 - frictionCoefficient))) * -1;
            particles[particleID].pos += impulse * fixedTimeStep * fixedTimeStep;

            BtoP = particles[particleID].pos - ball.pos;
            sqrDist = Vector2LengthSqr(BtoP);
            if (sqrDist >= ball.radius * ball.radius) {
                continue;
            }

            dist = sqrt(sqrDist);
            unitNormal = BtoP / dist;
            particles[particleID].pos += unitNormal * (ball.radius - dist);
        }

        // mouse interaction
        if (IsMouseButtonDown(MOUSE_LEFT_BUTTON) || IsMouseButtonDown(MOUSE_RIGHT_BUTTON)) {
            float interactForceMag = mouseInteractForce;
            if (IsMouseButtonDown(MOUSE_LEFT_BUTTON)) {
                interactForceMag = -interactForceMag;
            }

            for (int poolIndex = 0; poolIndex < activeCount; poolIndex++) {
                int particleID = objectPool[poolIndex];

                Vector2 mousePos = convertToSimPos(GetMousePosition());
                Vector2 MtoP = particles[particleID].pos - mousePos;

                float dist = Vector2Length(MtoP);
                if (dist <= mouseInteractRadius) {
                    Vector2 unitDir = MtoP / dist;
                    Vector2 forceVec = unitDir * interactForceMag;
                    float mass = particles[particleID].mass;
                    particles[particleID].pos += (forceVec * fixedTimeStep * fixedTimeStep) / mass;
                }
            }
        }

        if (IsKeyDown(KEY_ONE)) {
            Vector2 mousePos = convertToSimPos(GetMousePosition());

            float randAngle = (rand() * 2 * PI) / RAND_MAX;
            float randMag = (rand() * mouseInteractRadius) / RAND_MAX;
            Vector2 randOffset = { (float)cos(randAngle) * randMag, (float)sin(randAngle) * randMag };

            spawnParticle(defaultMass, mousePos + randOffset, { 0, 0 });
        }
        else if (IsKeyDown(KEY_TWO)) {
            for (int poolIndex = 0; poolIndex < activeCount; poolIndex++) {
                int particleID = objectPool[poolIndex];

                Vector2 mousePos = convertToSimPos(GetMousePosition());
                Vector2 MtoP = particles[particleID].pos - mousePos;

                float dist = Vector2Length(MtoP);
                if (dist <= mouseInteractRadius) {
                    despawnParticle(poolIndex);
                }
            }
        }
        

        // boundary collision check
        for (int poolIndex = 0; poolIndex < activeCount; poolIndex++) {
            int particleID = objectPool[poolIndex];
            
            if (particles[particleID].pos.y + particleRadius >= getScaledHeight()) {
                particles[particleID].pos.y = (getScaledHeight() - particleRadius);
            }
            else if (particles[particleID].pos.y - particleRadius <= 0) {
                particles[particleID].pos.y = particleRadius;
            }

            if (particles[particleID].pos.x - particleRadius <= 0) {
                particles[particleID].pos.x = particleRadius;
            }
            else if (particles[particleID].pos.x + particleRadius >= getScaledWidth()) {
                particles[particleID].pos.x = (getScaledWidth() - particleRadius);
            }
        }

        // compute implicit velocity
        for (int poolIndex = 0; poolIndex < activeCount; poolIndex++) {
            int particleID = objectPool[poolIndex];

            particles[particleID].vel = (particles[particleID].pos - previousPositions[particleID]) / fixedTimeStep;
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

        // shader logic
        rlEnableShader(clearTextureBufferProgram);
        rlBindShaderBuffer(textureSSBO, 1);
        rlComputeShaderDispatch((int)ceil((resolution.x * resolution.y) / WORKGROUP_SIZE), 1, 1);
        rlDisableShader();

        rlUpdateShaderBuffer(particleSSBO, particles.begin(), particles.getCount() * sizeof(Particle), 0);
        rlUpdateShaderBuffer(densitySSBO, densities.begin(), densities.getCount() * sizeof(float), 0);

        rlEnableShader(updateTextureBufferProgram);
        rlBindShaderBuffer(simDataSSBO, 1);
        rlBindShaderBuffer(particleSSBO, 2);
        rlBindShaderBuffer(densitySSBO, 3);
        rlBindShaderBuffer(textureSSBO, 4);
        int dispatches = MAX_PARTICLE_COUNT / WORKGROUP_SIZE;
        if (MAX_PARTICLE_COUNT % WORKGROUP_SIZE != 0) {
            dispatches++;
        }
        rlComputeShaderDispatch(dispatches, 1, 1);
        rlDisableShader();
    }

    // Iterated through all fixed updates
}

void Simulation::draw() {
    rlBindShaderBuffer(textureSSBO, 1);
    SetShaderValue(renderSimShader, resUniformLoc, &resolution, SHADER_UNIFORM_VEC2);

    BeginShaderMode(renderSimShader);
    DrawTexture(texture, bounds.x, bounds.y, WHITE);
    EndShaderMode();

    Vector2 mousePos = GetMousePosition();

    if (showSmoothingRadius) {
        for (int i = 0; i < particles.getCount(); i++) {
            Vector2 screenPos = convertToScreenPos(particles[i].pos);
            DrawCircle(screenPos.x, screenPos.y, smoothingRadius * scale, { 0, 121, 241, 31 });
        }
    }

    //for (int i = 0; i < particles.getCount(); i++) {
    //    Vector2 screenPos = convertToScreenPos(particles[i].pos);
    //    float densityError = targetDensity / densities[i];
    //    
    //    float r = (1 - abs(0.4 - densityError)) * 255;
    //    float g = (1 - abs(0.7 - densityError)) * 255;
    //    float b = densityError * 255;

    //    //float densityError = std::max((densities[i] - targetDensity), 0.f) * 255 * 2;

    //    //float r = std::max(0.f, std::min(255.f, densityError));
    //    //float g = std::max(0.f, std::min(255.f, densityError - 255));
    //    //float b = std::max(0.f, std::min(255.f, densityError - 511));

    //    Color densityColor = { r, g, b, 255 };
    //    
    //    //DrawCircle(screenPos.x, screenPos.y, particles[i].radius * scale, densityColor);
    //    DrawRectangle(screenPos.x - (particles[i].radius * scale), screenPos.y - (particles[i].radius * scale), 2 * particles[i].radius * scale, 2 * particles[i].radius * scale, densityColor);
    //}

    // spatial hash testing
    //spatialHash.draw({ bounds.x, bounds.y }, scale, convertToSimPos(mousePos));
    //Array<int> nearbyIDs = spatialHash.findNearby(spatialHash.getCellPos(convertToSimPos(mousePos)));
    //Array<Vector2> nearbyPositions(nearbyIDs.getCount());
    //for (int i = 0; i < nearbyIDs.getCount(); i++) {
    //    Vector2 screenPos = convertToScreenPos(projectedPositions[nearbyIDs[i]]);
    //    DrawCircle(screenPos.x, screenPos.y, defaultRadius * scale * 2, WHITE);
    //}

    for (int i = 0; i < airflows.getCount(); i++) {
        Vector2 airflowScreenPos = convertToScreenPos(airflows[i].pos);
        DrawRectangleLines(airflowScreenPos.x, airflowScreenPos.y, airflows[i].width * scale, airflows[i].height * scale, BLACK);
    }

    Vector2 ballScreenPos = convertToScreenPos(ball.pos);
    DrawCircleLines(ballScreenPos.x, ballScreenPos.y, ball.radius * scale, BLACK);

    DrawCircleLines(mousePos.x, mousePos.y, mouseInteractRadius * scale, BLACK);

    DrawRectangleLines(bounds.x, bounds.y, bounds.z - bounds.x, bounds.w - bounds.y, BLACK);
}
