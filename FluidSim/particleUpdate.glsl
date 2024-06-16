#version 430

// Updating particle positions

struct updateInfo {
    vec2 pos;
    vec2 vel;
    float deltaTime;
    float timeMultiplier;
}

layout (local_size_x = 64, local_size_y = 1, local_size_z = 1) in;

layout(std430, binding = 1) readonly restrict buffer particleBufferLayout {
    updateInfo particleBuffer[];
};

layout(std430, binding = 2) writeonly restrict buffer particleBufferLayout2 {
    updateInfo particleBufferDest[];
};

void main() {
    uint id = gl_GlobalInvocationID.x;
    vec2 newPos = particleBuffer[id].pos + particleBuffer[id].vel * particleBuffer[id].deltaTime * particleBuffer[id].timeMultiplier;
    particleBufferDest[id] = newPos;
}