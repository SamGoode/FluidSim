#version 430

// Simulation rendering shader

#define SIM_WIDTH 800
#define SIM_HEIGHT 600

// Input vertex attributes
in vec2 fragTexCoord;

// Output fragment color
out vec4 finalColor;

// Input texture buffer
layout(std430, binding = 1) readonly buffer TextureBuffer {
    vec4 pixels[SIM_WIDTH * SIM_HEIGHT];
} textureBuffer;

// Output resolution
uniform vec2 resolution;

void main() {
    ivec2 coords = ivec2(fragTexCoord * resolution);
    
    //finalColor = vec4(1, 1, 0, 1);
    //finalColor = vec4(coords.x/resolution.x, 0, coords.y/resolution.y, 1);
    //vec4 test = textureBuffer.pixels[coords.x + coords.y * uvec2(resolution).x];
    finalColor = textureBuffer.pixels[coords.x + coords.y * uvec2(resolution).x];
}