#version 430

// clearing texture buffer

// process 512 pixels at a time
#define WORKGROUP_SIZE 512

#define SIM_WIDTH 800
#define SIM_HEIGHT 600

layout (local_size_x = WORKGROUP_SIZE, local_size_y = 1, local_size_z = 1) in;

layout(std430, binding = 1) buffer TextureBuffer {
    vec4 pixels[SIM_WIDTH * SIM_HEIGHT];
} textureBuffer;

void main() {
    uint index = gl_GlobalInvocationID.x;

    textureBuffer.pixels[index] = vec4(0, 0, 0, 0);
}