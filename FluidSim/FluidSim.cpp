#include "raylib.h"
#include "rlgl.h"
#include "Simulation.h"
#include "Particle.h"

#define BUFFER_SIZE 64

typedef struct particleInfo {
    Vector2 pos;
    Vector2 vel;
    float deltaTime;
    float timeMultiplier;
};

//class ComputeShader {
//public:
//    unsigned int ID;
//
//public:
//    ComputeShader(const char* computePath) {
//        unsigned int compute;
//        compute = glCreateShader(GL_COMPUTE_SHADER);
//    }
//};

int main() {
    char* particleUpdateCode = LoadFileText("particleUpdate.glsl");
    unsigned int particleUpdateShader = rlCompileShader(particleUpdateCode, RL_COMPUTE_SHADER);
    unsigned int particleUpdateProgram = rlLoadComputeShaderProgram(particleUpdateShader);
    UnloadFileText(particleUpdateCode);

    unsigned int ssbo = rlLoadShaderBuffer(BUFFER_SIZE * sizeof(particleInfo), NULL, RL_DYNAMIC_COPY);

    int screenWidth = 1600;
    int screenHeight = 800;

    InitWindow(screenWidth, screenHeight, "Smooth Particle Hydrodynamics Sim");

    //SetTargetFPS(120);

    int seed = 844134593;
    srand(seed);

    Simulation sim({ 500, 200, (float)screenWidth - 500, (float)screenHeight - 200 });
    sim.update(0.000001);

    bool paused = true;

    while (!WindowShouldClose()) {
        // Updates
        float delta = GetFrameTime();

        if (IsKeyReleased(KEY_SPACE)) {
            paused = !paused;
        }

        Vector2 mousePos = GetMousePosition();
        Vector4 bounds = sim.getBounds();
        Vector2 mouseSimPos = { mousePos.x - bounds.x, mousePos.y - bounds.y };

        if (!paused) {
            sim.update(delta);
        }

        // Drawing
        BeginDrawing();

        ClearBackground(RAYWHITE);

        sim.draw();

        DrawFPS(10, 10);

        std::string mouseSimInfo = std::to_string(mouseSimPos.x) + "," + std::to_string(mouseSimPos.y);
        DrawText(mouseSimInfo.c_str(), 10, 50, 20, BLUE);
        std::string smoothing = std::to_string(sim.smoothing(50, 0));
        DrawText(smoothing.c_str(), 10, 80, 20, BLUE);
        //std::string density = std::to_string(sim.calculateDensity(mouseSimPos));
        //DrawText(density.c_str(), 10, 110, 20, BLUE);
        //Vector2 gradientVec = sim.calculateGradientVec(mouseSimPos);
        //std::string gradient = std::to_string(gradientVec.x) + "," + std::to_string(gradientVec.y);
        //DrawText(gradient.c_str(), 10, 130, 20, BLUE);

        EndDrawing();
    }

    //rlUnloadShaderProgram(particleUpdateProgram);
}
