#version 430

// particle physics update

#define BUFFER_SIZE 16384
// process 512 particles at a time
#define WORKGROUP_SIZE 512

#define SIM_WIDTH 800
#define SIM_HEIGHT 600
#define SCALE 4

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

layout(std430, binding = 3) buffer TextureBuffer {
    vec4 pixels[SIM_WIDTH * SIM_HEIGHT];
} textureBuffer;

void main() {
    uint index = gl_GlobalInvocationID.x;

    Particle particle = particleBuffer.particles[index];
    particleBuffer.particles[index].pos = particle.pos + particle.vel * simData.deltaTime * simData.timeDilation;
    //particleBuffer.particles[index].pos = particleBuffer.particles[index].pos + particleBuffer.particles[index].vel * simData.deltaTime * simData.timeDilation;

    ivec2 pos = ivec2(particleBuffer.particles[index].pos * SCALE);
    textureBuffer.pixels[(pos.x - 1) + (pos.y - 1) * SIM_WIDTH] = vec4(0, 0.5, 1, 1);
    textureBuffer.pixels[(pos.x + 0) + (pos.y - 1) * SIM_WIDTH] = vec4(0, 0.5, 1, 1);
    textureBuffer.pixels[(pos.x + 1) + (pos.y - 1) * SIM_WIDTH] = vec4(0, 0.5, 1, 1);
    textureBuffer.pixels[(pos.x - 1) + (pos.y + 0) * SIM_WIDTH] = vec4(0, 0.5, 1, 1);
    textureBuffer.pixels[(pos.x + 0) + (pos.y + 0) * SIM_WIDTH] = vec4(0, 0.5, 1, 1);
    textureBuffer.pixels[(pos.x + 1) + (pos.y + 0) * SIM_WIDTH] = vec4(0, 0.5, 1, 1);
    textureBuffer.pixels[(pos.x - 1) + (pos.y + 1) * SIM_WIDTH] = vec4(0, 0.5, 1, 1);
    textureBuffer.pixels[(pos.x + 0) + (pos.y + 1) * SIM_WIDTH] = vec4(0, 0.5, 1, 1);
    textureBuffer.pixels[(pos.x + 1) + (pos.y + 1) * SIM_WIDTH] = vec4(0, 0.5, 1, 1);
}