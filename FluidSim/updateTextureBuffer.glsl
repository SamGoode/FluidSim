#version 430

// updating texture buffer

#define MAX_PARTICLE_COUNT 16384
// process 512 particles at a time
#define WORKGROUP_SIZE 512

#define SIM_WIDTH 800
#define SIM_HEIGHT 600
#define SCALE 8

layout (local_size_x = WORKGROUP_SIZE, local_size_y = 1, local_size_z = 1) in;

struct Particle {
    int isActive;
    float mass;
    vec2 pos;
    vec2 vel;
};

layout(std430, binding = 1) buffer SimData {
    vec2 gravity;
    float targetDensity;
    float fixedTimeStep;
    int activeCount;
} simData;

//layout(std430, binding = 2) buffer ParticleBuffer {
//    Particle particles[MAX_PARTICLE_COUNT];
//} particleBuffer;

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

void main() {
    uint poolIndex = gl_GlobalInvocationID.x;
    
    if(poolIndex < simData.activeCount) {
        uint index = poolBuffer.particleIDs[poolIndex];

        float densityError = simData.targetDensity / densityBuffer.densities[index];
        float r = 1 - abs(0.4 - densityError);
        float g = 1 - abs(0.7 - densityError);
        float b = densityError;
        vec4 color = vec4(r, g, b, 1);

        ivec2 pos = ivec2(positionBuffer.positions[index] * SCALE);
        textureBuffer.pixels[(pos.x - 1) + (pos.y - 1) * SIM_WIDTH] = color;
        textureBuffer.pixels[(pos.x + 0) + (pos.y - 1) * SIM_WIDTH] = color;
        textureBuffer.pixels[(pos.x + 1) + (pos.y - 1) * SIM_WIDTH] = color;
        textureBuffer.pixels[(pos.x - 1) + (pos.y + 0) * SIM_WIDTH] = color;
        textureBuffer.pixels[(pos.x + 0) + (pos.y + 0) * SIM_WIDTH] = color;
        textureBuffer.pixels[(pos.x + 1) + (pos.y + 0) * SIM_WIDTH] = color;
        textureBuffer.pixels[(pos.x - 1) + (pos.y + 1) * SIM_WIDTH] = color;
        textureBuffer.pixels[(pos.x + 0) + (pos.y + 1) * SIM_WIDTH] = color;
        textureBuffer.pixels[(pos.x + 1) + (pos.y + 1) * SIM_WIDTH] = color;
    }
}