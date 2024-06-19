#version 430

// particle physics update

#define BUFFER_SIZE 16384
// process 1024 particles at a time
#define WORKGROUP_SIZE 1024

layout (local_size_x = WORKGROUP_SIZE, local_size_y = 1, local_size_z = 1) in;

struct particleInfo {
    vec2 pos;
    vec2 vel;
};

layout(std430, binding = 1) buffer ParticleBuffer {
    particleInfo particles[BUFFER_SIZE];
    vec2 gravity;
    float deltaTime;
    float projectedTime;
    float timeMultiplier;
} particleBuffer;

void main() {
    uint index = gl_GlobalInvocationID.x;

    particleBuffer.particles[index].vel = particleBuffer.particles[index].vel + particleBuffer.gravity * particleBuffer.deltaTime * particleBuffer.timeMultiplier;
    particleBuffer.particles[index].pos = particleBuffer.particles[index].pos + particleBuffer.particles[index].vel * particleBuffer.projectedTime * particleBuffer.timeMultiplier;
}