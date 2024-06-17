#version 430

// particle physics update

#define BUFFER_SIZE 256

// update 64 particles at a time
layout (local_size_x = 1, local_size_y = 1, local_size_z = 1) in;

struct particleInfo {
    vec2 pos;
    vec2 vel;
    float deltaTime;
};

layout(std430, binding = 1) readonly restrict buffer bufferLayout {
    particleInfo inBuffer[];
};

layout(std430, binding = 2) writeonly restrict buffer bufferLayout2 {
    particleInfo outBuffer[];
};

void main() {
    uint id = gl_GlobalInvocationID.x;

    outBuffer[id].pos = inBuffer[id].pos + inBuffer[id].vel * inBuffer[id].deltaTime;
    outBuffer[id].vel = vec2(0, 0);
    outBuffer[id].deltaTime = 0;
}