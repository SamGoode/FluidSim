#version 430

// updating texture buffer

#define MAX_PARTICLE_COUNT 16384

#define SIM_WIDTH 1000
#define SIM_HEIGHT 800
#define SCALE 4

layout (local_size_x = 8, local_size_y = 8, local_size_z = 16) in;

layout(std430, binding = 1) buffer SimData {
    vec2 gravity;
    float targetDensity;
    float fixedTimeStep;
    int activeCount;
    float particleRadius;
} simData;

layout(std430, binding = 2) buffer PositionBuffer {
    vec2 positions[MAX_PARTICLE_COUNT];
} positionBuffer;

layout(std430, binding = 3) buffer DensityBuffer {
    float densities[MAX_PARTICLE_COUNT];
} densityBuffer;

layout(std430, binding = 4) buffer TextureBuffer {
    vec4 pixels[SIM_WIDTH * SIM_HEIGHT];
} textureBuffer;

layout(std430, binding = 5) buffer PoolBuffer {
    int particleIDs[MAX_PARTICLE_COUNT];
} poolBuffer;

layout(std430, binding = 6) buffer HashListBuffer {
    ivec2 hashList[MAX_PARTICLE_COUNT];
} hashListBuffer;

layout(std430, binding = 7) buffer LookupBuffer {
    ivec2 lookup[(SIM_WIDTH / 8) * (SIM_HEIGHT / 8)];
} lookupBuffer;

void main() {
    ivec2 offsets[9] = {
        ivec2(-1, -1), ivec2(0, -1), ivec2(1, -1),
        ivec2(-1, 0), ivec2(0, 0), ivec2(1, 0),
        ivec2(-1, 1), ivec2(0, 1), ivec2(1, 1)
    };

    ivec2 cellPos = ivec2(gl_WorkGroupID.xy) + offsets[gl_WorkGroupID.z];

    if(cellPos.x < 0 || cellPos.x >= gl_NumWorkGroups.x || cellPos.y < 0 || cellPos.y >= gl_NumWorkGroups.y) {
        return;
    }

    //uint cellIndex = gl_WorkGroupID.x + gl_WorkGroupID.y * gl_NumWorkGroups.x;
    uint cellIndex = cellPos.x + cellPos.y * gl_NumWorkGroups.x;
    
    ivec2 pixelPos = ivec2(gl_GlobalInvocationID.xy);
    int localParticleIndex = int(gl_LocalInvocationID.z);

    int startIndex = lookupBuffer.lookup[cellIndex].x;
    int endIndex = lookupBuffer.lookup[cellIndex].y;
    int particleCount = (endIndex - startIndex) + 1;

    if(startIndex < 0) {
        return;
    }

    if(localParticleIndex < particleCount) {
        int particleIndex = hashListBuffer.hashList[startIndex + localParticleIndex].x;
        
        vec2 particlePos = positionBuffer.positions[particleIndex] * SCALE;

        vec2 particleToPixel = vec2(pixelPos) - particlePos;

        float dist = sqrt((particleToPixel.x * particleToPixel.x) + (particleToPixel.y * particleToPixel.y));

        if(dist < simData.particleRadius * SCALE) { // Replace with particle radius later
            float densityError = 1.5f / densityBuffer.densities[particleIndex];
            float r = 1 - abs(0.4 - densityError);
            float g = 1 - abs(0.7 - densityError);
            float b = densityError;
            vec4 color = vec4(r, g, b, 1);

            textureBuffer.pixels[pixelPos.x + pixelPos.y * SIM_WIDTH] = color;
        }
    }
}