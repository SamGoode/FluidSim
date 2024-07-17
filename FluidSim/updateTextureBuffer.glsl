#version 430

// updating texture buffer

#define MAX_PARTICLE_COUNT 4096
// process 512 particles at a time
#define WORKGROUP_SIZE 512

#define SIM_WIDTH 800
#define SIM_HEIGHT 600
#define SCALE 8

layout (local_size_x = 8, local_size_y = 8, local_size_z = 16) in;

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

layout(std430, binding = 6) buffer HashListBuffer {
    ivec2 hashList[MAX_PARTICLE_COUNT];
} hashListBuffer;

layout(std430, binding = 7) buffer LookupBuffer {
    ivec2 lookup[(SIM_WIDTH / 8) * (SIM_HEIGHT / 8)];
} lookupBuffer;

void main() {
    uint cellIndex = gl_WorkGroupID.x + gl_WorkGroupID.y * gl_NumWorkGroups.x;

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
        
        ivec2 particlePos = ivec2(positionBuffer.positions[particleIndex] * SCALE);

        ivec2 particleToPixel = pixelPos - particlePos;

        float dist = sqrt((particleToPixel.x * particleToPixel.x) + (particleToPixel.y * particleToPixel.y));

        if(dist < 0.5 * SCALE) { // Replace with particle radius later
            float densityError = 1.5f / densityBuffer.densities[particleIndex];
            float r = 1 - abs(0.4 - densityError);
            float g = 1 - abs(0.7 - densityError);
            float b = densityError;
            vec4 color = vec4(r, g, b, 1);

            textureBuffer.pixels[pixelPos.x + pixelPos.y * SIM_WIDTH] = color;
        }
    }

    //float value = float(particleCount) / 3;
    //vec4 testColor = vec4(value, value, value, 1);
    //textureBuffer.pixels[pixelPos.x + pixelPos.y * SIM_WIDTH] = testColor;

    //if(cellX == 2 && cellY == 1) {
        //textureBuffer.pixels[pixelPos.x + pixelPos.y * SIM_WIDTH] = vec4(0, 0, 0, 1);
    //}
    //float red = float(cellX)/gl_NumWorkGroups.x;
    //float green = float(cellY)/gl_NumWorkGroups.y;

    //vec4 testColor = vec4(red, green, 1, 1);
    //textureBuffer.pixels[pixelPos.x + pixelPos.y * SIM_WIDTH] = testColor;

    //uint poolIndex = gl_GlobalInvocationID.x;
    //
    //if (poolIndex < simData.activeCount) {
    //    uint index = poolBuffer.particleIDs[poolIndex];
    //
    //    float densityError = 1.5f / densityBuffer.densities[index];
    //    float r = 1 - abs(0.4 - densityError);
    //    float g = 1 - abs(0.7 - densityError);
    //    float b = densityError;
    //    vec4 color = vec4(r, g, b, 1);
    //
    //    ivec2 pos = ivec2(positionBuffer.positions[index] * SCALE);
    //    textureBuffer.pixels[(pos.x - 1) + (pos.y - 1) * SIM_WIDTH] = color;
    //    textureBuffer.pixels[(pos.x + 0) + (pos.y - 1) * SIM_WIDTH] = color;
    //    textureBuffer.pixels[(pos.x + 1) + (pos.y - 1) * SIM_WIDTH] = color;
    //    textureBuffer.pixels[(pos.x - 1) + (pos.y + 0) * SIM_WIDTH] = color;
    //    textureBuffer.pixels[(pos.x + 0) + (pos.y + 0) * SIM_WIDTH] = color;
    //    textureBuffer.pixels[(pos.x + 1) + (pos.y + 0) * SIM_WIDTH] = color;
    //    textureBuffer.pixels[(pos.x - 1) + (pos.y + 1) * SIM_WIDTH] = color;
    //    textureBuffer.pixels[(pos.x + 0) + (pos.y + 1) * SIM_WIDTH] = color;
    //    textureBuffer.pixels[(pos.x + 1) + (pos.y + 1) * SIM_WIDTH] = color;
    //}
}