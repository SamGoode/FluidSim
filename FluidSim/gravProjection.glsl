#version 430

// particle gravity application and projected positions

#define BUFFER_SIZE 16384
// process 512 particles at a time
#define WORKGROUP_SIZE 512

layout (local_size_x = WORKGROUP_SIZE, local_size_y = 1, local_size_z = 1) in;

struct Particle {
    float mass;
    float radius;
    vec2 pos;
    vec2 vel;
};

layout(std430, binding = 1) buffer ParticleBuffer {
    Particle particles[BUFFER_SIZE];
} particleBuffer;

layout(std430, binding = 2) readonly restrict buffer SimData {
    vec2 gravity;
    float deltaTime;
    float projectedTime;
    float timeDilation;
} simData;

layout(std430, binding = 3) writeonly restrict buffer ProjectedPositionBuffer {
    vec2 projectedPositions[BUFFER_SIZE];
} projectedPositionBuffer;

void main() {
    uint index = gl_GlobalInvocationID.x;

    particleBuffer.particles[index].vel = particleBuffer.particles[index].vel + simData.gravity * simData.deltaTime * simData.timeDilation;
    projectedPositionBuffer.projectedPositions[index] = particleBuffer.particles[index].pos + particleBuffer.particles[index].vel * simData.projectedTime * simData.timeDilation;
}