#include "Simulation.h"
#include "raymath.h"
#include "rlgl.h"
#include "RaylibOverloads.h"
#include <algorithm>

#define MAX_PARTICLE_COUNT 16384
#define WORKGROUP_SIZE 512

#define SIM_WIDTH 800
#define SIM_HEIGHT 600

Simulation::Simulation(Vector4 _bounds) {
    bounds = _bounds;
    // scale is pixels/unit
    scale = 6;
    resolution = { getWidth(), getHeight() };
    Image whiteImage = GenImageColor(getWidth(), getHeight(), WHITE);
    texture = LoadTextureFromImage(whiteImage);
    UnloadImage(whiteImage);

    gravity = { 0, 1 };
    frictionCoefficient = 0.5f;
    stickyDist = 1.f;
    stickyCoefficient = 0.1f;

    viscLinear = 0.f;
    viscQuad = 0.2f;

    showSmoothingRadius = false;
    smoothingRadius = 2.f;
    sqrRadius = smoothingRadius * smoothingRadius;
    particleRadius = 0.5f;
    targetDensity = 1.5f;
    pressureMultiplier = 10;
    nearPressureMultiplier = 20;
    timeDilation = 2.f;

    mouseInteractRadius = 8;
    mouseInteractForce = 10;

    spawnAmount = 5;
    spawnArea = {
        {0, 0, 2, 100},
        {8, 0}
    };

    despawnArea = {
        {131.33f, 0, 133.33f, 100}
    };

    airflows = Array<AirflowSpace>(0);

    balls = Array<Ball>(0);
    //balls[0] = {
    //    {67, 51},
    //    6
    //};
    //balls[1] = {
    //    {69, 51},
    //    7
    //};
    //balls[2] = {
    //    {72, 50.5},
    //    8
    //};
    //balls[3] = {
    //    {76, 50},
    //    9.5
    //};
    //balls[4] = {
    //    {80, 50},
    //    10
    //};
    //balls[5] = {
    //    {85, 50},
    //    10
    //};
    //balls[6] = {
    //    {89, 51},
    //    9
    //};
    //balls[7] = {
    //    {92, 52},
    //    8
    //};
    //balls[8] = {
    //    {95, 53},
    //    7
    //};
    //balls[9] = {
    //    {98, 54},
    //    6
    //};
    //balls[10] = {
    //    {100, 55},
    //    5
    //};

    fixedTimeStep = 0.02f;
    timePassed = 0;
    timeStepUpperBound = 0.1f;
    timeStepLowerBound = 0.015f;

    defaultMass = 1;
    densities = Array<float>(MAX_PARTICLE_COUNT);
    nearDensities = Array<float>(MAX_PARTICLE_COUNT);
    previousPositions = Array<Vector2>(MAX_PARTICLE_COUNT);
    positions = Array<Vector2>(MAX_PARTICLE_COUNT);
    velocities = Array<Vector2>(MAX_PARTICLE_COUNT);
    masses = Array<float>(MAX_PARTICLE_COUNT);

    objectPool = Array<int>(MAX_PARTICLE_COUNT);
    for (int i = 0; i < MAX_PARTICLE_COUNT; i++) {
        objectPool[i] = i;
    }
    activeCount = 4000;
    
    // initiate particles at random positions and with random velocities
    for (int i = 0; i < activeCount; i++) {
        // generate random location
        float x = (rand() * (bounds.z - bounds.x)/scale) / RAND_MAX;
        float y = (rand() * (bounds.w - bounds.y)/scale) / RAND_MAX;
        Vector2 position = { x, y };

        positions[i] = { x, y };

        // generate random velocity
        float velX = ((rand() * 10) / RAND_MAX) - 5;
        float velY = ((rand() * 10) / RAND_MAX) - 5;
        //Vector2 vel = { velX, velY };
        Vector2 vel = { 0, 0 };

        velocities[i] = { 0, 0 };
        masses[i] = defaultMass;

    }

    maxVelocity = 0;
    for (int i = 0; i < activeCount; i++) {
        float velMag = Vector2Length(velocities[i]);
        if (maxVelocity < velMag) {
            maxVelocity = velMag;
        }
    }
    
    // initiate particles in a square formation
    //for (int i = 0; i < particles.getCount(); i++) {
    //    particles[i] = { defaultMass, defaultRadius, {10 + (float)(i % 128) * defaultRadius * 1.2f, 5 + (float)(i / 128) * defaultRadius * 1.2f}, {0, 0} };
    //}

    for (int i = 0; i < activeCount; i++) {
        previousPositions[i] = positions[i];
    }

    spatialHash = SpatialHashGrid({getScaledWidth(), getScaledHeight()}, smoothingRadius, smoothingRadius);
    // 3x3 grid around center cell
    for (int x = 0; x < 3; x++) {
        for (int y = 0; y < 3; y++) {
            cellOffsets[x + y * 3] = { x - 1, y - 1 };
        }
    }

    spatialRendering = SpatialHashGrid({ getWidth(), getHeight() }, 8.f, 8.f);
    unscaledPositions = Array<Vector2>(MAX_PARTICLE_COUNT);

    char* updateTextureBufferCode = LoadFileText("updateTextureBuffer.glsl");
    unsigned int updateTextureBufferShader = rlCompileShader(updateTextureBufferCode, RL_COMPUTE_SHADER);
    updateTextureBufferProgram = rlLoadComputeShaderProgram(updateTextureBufferShader);
    UnloadFileText(updateTextureBufferCode);

    renderSimShader = LoadShader(NULL, "renderSim.glsl");
    resUniformLoc = GetShaderLocation(renderSimShader, "resolution");

    simDataSSBO = rlLoadShaderBuffer(sizeof(SimData), NULL, RL_DYNAMIC_COPY);
    poolSSBO = rlLoadShaderBuffer(MAX_PARTICLE_COUNT * sizeof(int), NULL, RL_DYNAMIC_COPY);
    positionSSBO = rlLoadShaderBuffer(MAX_PARTICLE_COUNT * sizeof(Vector2), NULL, RL_DYNAMIC_COPY);
    densitySSBO = rlLoadShaderBuffer(densities.getCount() * sizeof(float), NULL, RL_DYNAMIC_COPY);

    int pixelCount = getWidth() * getHeight();
    Array<Vector4> initialTextureData = Array<Vector4>(pixelCount);
    for (int i = 0; i < pixelCount; i++) {
        initialTextureData[i] = Vector4(0, 0, 0, 0);
    }

    defaultTextureSSBO = rlLoadShaderBuffer(getWidth() * getHeight() * sizeof(Vector4), initialTextureData.begin(), RL_DYNAMIC_COPY);
    textureSSBO = rlLoadShaderBuffer(getWidth() * getHeight() * sizeof(Vector4), NULL, RL_DYNAMIC_COPY);
    
    hashListSSBO = rlLoadShaderBuffer(MAX_PARTICLE_COUNT * sizeof(int2), NULL, RL_DYNAMIC_COPY);
    lookupSSBO = rlLoadShaderBuffer((SIM_WIDTH/8) * (SIM_HEIGHT/8) * sizeof(int2), NULL, RL_DYNAMIC_COPY);

    simData = {
        gravity,
        targetDensity,
        fixedTimeStep,
        activeCount,
        particleRadius
    };
    rlUpdateShaderBuffer(simDataSSBO, &simData, sizeof(simData), 0);
}

Simulation::~Simulation() {
    rlUnloadShaderBuffer(simDataSSBO);
    rlUnloadShaderBuffer(positionSSBO);
    rlUnloadShaderBuffer(poolSSBO);

    rlUnloadShaderBuffer(hashListSSBO);
    rlUnloadShaderBuffer(lookupSSBO);
    
    rlUnloadShaderBuffer(projectedPositionSSBO);
    rlUnloadShaderBuffer(densitySSBO);
    rlUnloadShaderBuffer(textureSSBO);

    rlUnloadShaderProgram(updateTextureBufferProgram);

    UnloadTexture(texture);
    UnloadShader(renderSimShader);
}

// density kernel
// dist should never exceed radius
float Simulation::densityKernel(float radius, float dist) {
    float value = 1 - (dist / radius);
    return value * value;
}

float Simulation::nearDensityKernel(float radius, float dist) {
    float value = 1 - (dist / radius);
    return value * value * value;
}

// calculates density level at a particle's position
// (density, nearDensity)
std::pair<float, float> Simulation::calculateDensity(int particleID) {
    const Array<int2>& hashList = spatialHash.getHashList();
    const Array<int2>& indexLookup = spatialHash.getIndexLookup();

    Vector2 pos = positions[particleID];
    int2 cellPos = spatialHash.getCellPos(pos);

    float density = 0;
    float nearDensity = 0;
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
            Vector2 otherPos = positions[otherParticleID];

            float sqrDist = Vector2DistanceSqr(pos, otherPos);

            if (sqrDist > sqrRadius) {
                continue;
            }

            float dist = sqrt(sqrDist);
            density += densityKernel(smoothingRadius, dist) * masses[otherParticleID];
            nearDensity += nearDensityKernel(smoothingRadius, dist) * masses[otherParticleID];
        }
    }

    return { density, nearDensity };
}

float Simulation::convertDensityToPressure(float density) {
    float densityError = density - targetDensity;
    float pressure = densityError * pressureMultiplier;
    return pressure;
}

Vector2 Simulation::calculatePressureForce(int particleID) {
    const Array<int2>& hashList = spatialHash.getHashList();
    const Array<int2>& indexLookup = spatialHash.getIndexLookup();

    Vector2 pos = positions[particleID];
    int2 cellPos = spatialHash.getCellPos(pos);

    Vector2 pressureForce = { 0, 0 };
    for (int i = 0; i < 9; i++) {
        int2 offsetCellPos = cellPos + cellOffsets[i];

        if (!spatialHash.isValidCellPos(offsetCellPos)) {
            continue;
        }

        int offsetCellHash = spatialHash.getCellHash(offsetCellPos);

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

            Vector2 otherPos = positions[otherParticleID];
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

            float pressure = (densities[otherParticleID] - targetDensity) * pressureMultiplier;
            float nearPressure = nearDensities[otherParticleID] * nearPressureMultiplier;

            float weight = 1 - (dist / smoothingRadius);
            
            // BIG PROBLEM HERE TO FIX, THIS DOESN'T CONSERVE LINEAR MOMENTUM
            pressureForce -= dir * (pressure * weight + nearPressure * weight * weight) / 2;
        }
    }

    return pressureForce;
}

void Simulation::applyPressureDisplacements(int particleID, float deltaTime) {
    const Array<int2>& hashList = spatialHash.getHashList();
    const Array<int2>& indexLookup = spatialHash.getIndexLookup();

    Vector2 pos = positions[particleID];
    int2 cellPos = spatialHash.getCellPos(pos);

    Vector2 pressureForceSum = { 0, 0 };
    for (int i = 0; i < 9; i++) {
        int2 offsetCellPos = cellPos + cellOffsets[i];

        if (!spatialHash.isValidCellPos(offsetCellPos)) {
            continue;
        }

        int offsetCellHash = spatialHash.getCellHash(offsetCellPos);

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

            Vector2 otherPos = positions[otherParticleID];
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

            float pressure = (densities[otherParticleID] - targetDensity) * pressureMultiplier;
            float nearPressure = nearDensities[otherParticleID] * nearPressureMultiplier;

            float weight = 1 - (dist / smoothingRadius);

            float sharedPressure = (pressure * weight + nearPressure * weight * weight);
            float pressureForce = sharedPressure / 2;

            // BIG PROBLEM HERE TO FIX, THIS DOESN'T CONSERVE LINEAR MOMENTUM
            //pressureForce -= dir * (pressure * weight + nearPressure * weight * weight) / 2;
            positions[otherParticleID] += dir * (deltaTime * deltaTime * pressureForce / masses[otherParticleID]);
            pressureForceSum -= dir * pressureForce;
        }
    }

    positions[particleID] += pressureForceSum * (deltaTime * deltaTime / masses[particleID]);
}

Vector2 Simulation::calculateViscosityImpulse(int particleID) {
    const Array<int2>& hashList = spatialHash.getHashList();
    const Array<int2>& indexLookup = spatialHash.getIndexLookup();

    Vector2 pos = positions[particleID];
    int2 cellPos = spatialHash.getCellPos(pos);

    Vector2 viscosityImpulse = { 0, 0 };
    for (int i = 0; i < 9; i++) {
        int2 offsetCellPos = cellPos + cellOffsets[i];

        if (!spatialHash.isValidCellPos(offsetCellPos)) {
            continue;
        }

        int offsetCellHash = spatialHash.getCellHash(offsetCellPos);

        int startIndex = indexLookup[offsetCellHash].x;
        if (startIndex < 0) {
            continue;
        }

        int endIndex = indexLookup[offsetCellHash].y;

        for (int n = startIndex; n < endIndex + 1; n++) {
            int otherParticleID = hashList[n].x;
            //if (particleID >= otherParticleID) {
            //    continue;
            //}

            Vector2 otherPos = positions[otherParticleID];

            Vector2 toOther = otherPos - pos;
            float sqrDist = Vector2LengthSqr(toOther);

            if (sqrDist == 0 || sqrDist > sqrRadius) {
                continue;
            }

            float dist = sqrt(sqrDist);
            Vector2 normal = toOther / dist;

            // inward radial velocity
            float inRadVel = Vector2DotProduct((velocities[particleID] - velocities[otherParticleID]), normal);
            if (inRadVel <= 0) {
                continue;
            }

            Vector2 impulse = normal * ((1 - (dist / smoothingRadius)) * ((viscLinear * inRadVel) + (viscQuad * (inRadVel * inRadVel))));

            velocities[particleID] -= impulse / 2;
            //velocities[otherParticleID] += impulse / 2;
        }
    }

    return viscosityImpulse;
}

void Simulation::spawnParticle(float mass, Vector2 pos, Vector2 vel) {
    if (activeCount >= objectPool.getCount()) {
        return;
    }

    activeCount++;
    int particleID = objectPool[activeCount - 1];

    previousPositions[particleID] = pos;
    masses[particleID] = mass;
    positions[particleID] = pos;
    velocities[particleID] = vel;
}

void Simulation::despawnParticle(int poolIndex) {
    int particleID = objectPool[poolIndex];

    activeCount--;

    objectPool[poolIndex] = objectPool[activeCount];
    objectPool[activeCount] = particleID;
}

// time step = scaling * (particle diameter / magnitude of highest particle velocity)
float Simulation::calculateTimeStepSize() {
    if (maxVelocity > 0) {
        float maxTimeStep = (2 * particleRadius) / maxVelocity;
        float scaling = 0.9f;
        float value = maxTimeStep * scaling;

        return std::clamp(value, timeStepLowerBound, timeStepLowerBound);
    }

    return timeStepLowerBound;
}

void Simulation::update(float deltaTime) {
    if (deltaTime > 0.05f) {
        deltaTime = 0;
    }

    timePassed += deltaTime * timeDilation;


    float nextTimeStepSize = calculateTimeStepSize();
    while (timePassed > nextTimeStepSize) {
        stepForward();
        timePassed -= nextTimeStepSize;
        nextTimeStepSize = calculateTimeStepSize();
    }
}

void Simulation::stepForward() {
    float timeStep = calculateTimeStepSize();

    // apply gravity
    for (int poolIndex = 0; poolIndex < activeCount; poolIndex++) {
        int particleID = objectPool[poolIndex];
        velocities[particleID] += gravity * timeStep;
    }

    // spatial hash current positions
    spatialHash.generateHashList(positions, objectPool, activeCount);
    spatialHash.sortByCellHash();
    spatialHash.generateLookup();

    // apply viscosity
    for (int poolIndex = 0; poolIndex < activeCount; poolIndex++) {
        int particleID = objectPool[poolIndex];

        Vector2 viscImpulse = calculateViscosityImpulse(particleID);
        velocities[particleID] += viscImpulse * timeStep;
    }

    // save current positions and shift particle positions to projected positions
    for (int poolIndex = 0; poolIndex < activeCount; poolIndex++) {
        int particleID = objectPool[poolIndex];
        previousPositions[particleID] = positions[particleID];
        positions[particleID] += velocities[particleID] * timeStep;
    }

    // spatial hash current positions
    spatialHash.generateHashList(positions, objectPool, activeCount);
    spatialHash.sortByCellHash();
    spatialHash.generateLookup();

    // update densities and near densities cache
    for (int poolIndex = 0; poolIndex < activeCount; poolIndex++) {
        int particleID = objectPool[poolIndex];
        std::pair<float, float> densityPair = calculateDensity(particleID);
        densities[particleID] = densityPair.first;
        nearDensities[particleID] = densityPair.second;
    }

    //// calculate cumulative pressure forces and apply them
    //for (int poolIndex = 0; poolIndex < activeCount; poolIndex++) {
    //    int particleID = objectPool[poolIndex];

    //    float mass = masses[particleID];
    //    Vector2 pressureForce = calculatePressureForce(particleID);
    //    positions[particleID] += pressureForce * ((timeStep * timeStep) / mass);
    //}

    // calculate cumulative pressure forces and apply them
    for (int poolIndex = 0; poolIndex < activeCount; poolIndex++) {
        int particleID = objectPool[poolIndex];

        applyPressureDisplacements(particleID, timeStep);
    }

    //// particle spawning
    //for (int i = 0; i < spawnAmount; i++) {
    //    Vector2 spawnPos = { spawnArea.bounds.x + ((rand() * spawnArea.getWidth()) / RAND_MAX), spawnArea.bounds.y + ((rand() * spawnArea.getHeight()) / RAND_MAX) };

    //    spawnParticle(defaultMass, spawnPos, spawnArea.initVel);
    //}

    //// particle despawning
    //for (int poolIndex = 0; poolIndex < activeCount; poolIndex++) {
    //    int particleID = objectPool[poolIndex];
    //    Vector2 particlePos = positions[particleID];

    //    if (particlePos.x < despawnArea.bounds.x || particlePos.x > despawnArea.bounds.z || particlePos.y < despawnArea.bounds.y || particlePos.y > despawnArea.bounds.w) {
    //        continue;
    //    }

    //    despawnParticle(poolIndex);
    //}

    // airflow spaces
    for (int i = 0; i < airflows.getCount(); i++) {
        for (int poolIndex = 0; poolIndex < activeCount; poolIndex++) {
            int particleID = objectPool[poolIndex];

            Vector2 particlePos = positions[particleID];
            float mass = masses[particleID];

            if (particlePos.x < airflows[i].pos.x || particlePos.x > airflows[i].pos.x + airflows[i].width || particlePos.y < airflows[i].pos.y || particlePos.y > airflows[i].pos.y + airflows[i].height) {
                continue;
            }

            positions[particleID] += airflows[i].force * (timeStep * timeStep / mass);
        }
    }

    // ball collisions
    for (int i = 0; i < balls.getCount(); i++) {
        for (int poolIndex = 0; poolIndex < activeCount; poolIndex++) {
            int particleID = objectPool[poolIndex];

            Vector2 BtoP = positions[particleID] - balls[i].pos;
            float sqrDist = Vector2LengthSqr(BtoP);

            float dist = sqrt(sqrDist);
            Vector2 unitNormal = BtoP / dist;

            float surfDist = dist - balls[i].radius;
            if (surfDist < stickyDist) {
                Vector2 impulseStick = unitNormal * -stickyCoefficient * surfDist * (1 - (surfDist / stickyDist));
                positions[particleID] += impulseStick * timeStep;
            }

            if (sqrDist >= balls[i].radius * balls[i].radius) {
                continue;
            }

            //float dist = sqrt(sqrDist);
            //Vector2 unitNormal = BtoP / dist;

            Vector2 velocity = (positions[particleID] - previousPositions[particleID]) / timeStep;
            float velNormalMag = Vector2DotProduct(velocity, unitNormal);
            Vector2 velNormal = unitNormal * velNormalMag;
            Vector2 velTangent = velocity - velNormal;

            Vector2 impulse = (velNormal - (velTangent * (1 - frictionCoefficient))) * -1;
            positions[particleID] += impulse * timeStep * timeStep;

            BtoP = positions[particleID] - balls[i].pos;
            sqrDist = Vector2LengthSqr(BtoP);
            if (sqrDist >= balls[i].radius * balls[i].radius) {
                continue;
            }

            dist = sqrt(sqrDist);
            unitNormal = BtoP / dist;
            positions[particleID] += unitNormal * (balls[i].radius - dist);
        }
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
            Vector2 MtoP = positions[particleID] - mousePos;

            float dist = Vector2Length(MtoP);
            if (dist <= mouseInteractRadius) {
                Vector2 unitDir = MtoP / dist;
                Vector2 forceVec = unitDir * interactForceMag;
                float mass = masses[particleID];
                positions[particleID] += forceVec * (timeStep * timeStep / mass);
            }
        }
    }

    // mouse spawning and despawning
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
            Vector2 MtoP = positions[particleID] - mousePos;

            float dist = Vector2Length(MtoP);
            if (dist <= mouseInteractRadius) {
                despawnParticle(poolIndex);
            }
        }
    }

    // boundary collision check
    for (int poolIndex = 0; poolIndex < activeCount; poolIndex++) {
        int particleID = objectPool[poolIndex];

        if (positions[particleID].y + smoothingRadius >= getScaledHeight()) {
            //positions[particleID].y = (getScaledHeight() - smoothingRadius);
            
            float dist = getScaledHeight() - positions[particleID].y;
            float value = (1 - (dist / smoothingRadius));
            float pressure = value * value;
            positions[particleID].y -= pressure * value * timeStep;
        }
        else if (positions[particleID].y - smoothingRadius <= 0) {
            //positions[particleID].y = smoothingRadius;

            float dist = positions[particleID].y;
            float value = (1 - (dist / smoothingRadius));
            float pressure = value * value;
            positions[particleID].y += pressure * value * timeStep;
        }

        if (positions[particleID].x - smoothingRadius <= 0) {
            //positions[particleID].x = smoothingRadius;

            float dist = positions[particleID].x;
            float value = (1 - (dist / smoothingRadius));
            float pressure = value * value;
            positions[particleID].x += pressure * value * timeStep;
        }
        else if (positions[particleID].x + smoothingRadius >= getScaledWidth()) {
            //positions[particleID].x = (getScaledWidth() - smoothingRadius);

            float dist = getScaledWidth() - positions[particleID].x;
            float value = (1 - (dist / smoothingRadius));
            float pressure = value * value;
            positions[particleID].x -= pressure * value * timeStep;
        }
    }

    // compute implicit velocity
    maxVelocity = 0;
    for (int poolIndex = 0; poolIndex < activeCount; poolIndex++) {
        int particleID = objectPool[poolIndex];

        Vector2 vel = (positions[particleID] - previousPositions[particleID]) / timeStep;
        velocities[particleID] = vel;
        
        float velMag = Vector2Length(vel);
        if (maxVelocity < velMag) {
            maxVelocity = velMag;
        }
    }

    // reset texture buffer
    rlCopyShaderBuffer(textureSSBO, defaultTextureSSBO, 0, 0, getWidth() * getHeight() * sizeof(Vector4));

    // subdivide the pixel space into 8x8 squares of pixels assuming each contains a max of 16 particles

    // unscaling first
    for (int i = 0; i < MAX_PARTICLE_COUNT; i++) {
        Vector2 pos = positions[i];
        unscaledPositions[i] = pos * scale;
    }

    spatialRendering.generateHashList(unscaledPositions, objectPool, activeCount);
    spatialRendering.sortByCellHash();
    spatialRendering.generateLookup();

    simData.activeCount = activeCount;
    rlUpdateShaderBuffer(simDataSSBO, &simData, sizeof(simData), 0);
    rlUpdateShaderBuffer(poolSSBO, objectPool.begin(), MAX_PARTICLE_COUNT * sizeof(int), 0);
    rlUpdateShaderBuffer(positionSSBO, positions.begin(), MAX_PARTICLE_COUNT * sizeof(Vector2), 0);
    rlUpdateShaderBuffer(densitySSBO, densities.begin(), MAX_PARTICLE_COUNT * sizeof(float), 0);

    rlUpdateShaderBuffer(hashListSSBO, spatialRendering.getHashList().begin(), activeCount * sizeof(int2), 0);
    rlUpdateShaderBuffer(lookupSSBO, spatialRendering.getIndexLookup().begin(), (SIM_WIDTH / 8)* (SIM_HEIGHT / 8) * sizeof(int2), 0);

    rlEnableShader(updateTextureBufferProgram);
    rlBindShaderBuffer(simDataSSBO, 1);
    rlBindShaderBuffer(positionSSBO, 2);
    rlBindShaderBuffer(densitySSBO, 3);
    rlBindShaderBuffer(textureSSBO, 4);
    rlBindShaderBuffer(poolSSBO, 5);

    rlBindShaderBuffer(hashListSSBO, 6);
    rlBindShaderBuffer(lookupSSBO, 7);

    rlComputeShaderDispatch(SIM_WIDTH / 8, SIM_HEIGHT / 8, 9);

    rlDisableShader();
}

void Simulation::draw() {
    rlBindShaderBuffer(textureSSBO, 1);
    SetShaderValue(renderSimShader, resUniformLoc, &resolution, SHADER_UNIFORM_VEC2);

    BeginShaderMode(renderSimShader);
    DrawTexture(texture, bounds.x, bounds.y, WHITE);
    EndShaderMode();

    Vector2 mousePos = GetMousePosition();

    //if (showSmoothingRadius) {
    //    for (int i = 0; i < particles.getCount(); i++) {
    //        Vector2 screenPos = convertToScreenPos(particles[i].pos);
    //        DrawCircle(screenPos.x, screenPos.y, smoothingRadius * scale, { 0, 121, 241, 31 });
    //    }
    //}

    Vector2 spawnScreenPos = convertToScreenPos(spawnArea.getStartPos());
    DrawRectangleLines(spawnScreenPos.x, spawnScreenPos.y, spawnArea.getWidth() * scale, spawnArea.getHeight() * scale, BLACK);

    Vector2 despawnScreenPos = convertToScreenPos(despawnArea.getStartPos());
    DrawRectangleLines(despawnScreenPos.x, despawnScreenPos.y, despawnArea.getWidth() * scale, despawnArea.getHeight() * scale, BLACK);

    for (int i = 0; i < airflows.getCount(); i++) {
        Vector2 airflowScreenPos = convertToScreenPos(airflows[i].pos);
        DrawRectangleLines(airflowScreenPos.x, airflowScreenPos.y, airflows[i].width * scale, airflows[i].height * scale, BLACK);
    }

    for (int i = 0; i < balls.getCount(); i++) {
        Vector2 ballScreenPos = convertToScreenPos(balls[i].pos);
        //DrawCircleLines(ballScreenPos.x, ballScreenPos.y, balls[i].radius * scale, BLACK);
        DrawCircle(ballScreenPos.x, ballScreenPos.y, balls[i].radius * scale, LIGHTGRAY);
    }

    DrawCircleLines(mousePos.x, mousePos.y, mouseInteractRadius * scale, BLACK);

    DrawRectangleLines(bounds.x, bounds.y, bounds.z - bounds.x, bounds.w - bounds.y, BLACK);
}
