#include "raylib.h"
#include "rlgl.h"
#include "Simulation.h"

int main() {
    int screenWidth = 1600;
    int screenHeight = 1000;

    InitWindow(screenWidth, screenHeight, "Smooth Particle Hydrodynamics Sim");

    //SetTargetFPS(240);

    int seed = 844134593;
    srand(seed);

    Simulation sim({ 300, 100, (float)screenWidth - 300, (float)screenHeight - 100 });
    // update once to initialise all values
    sim.stepForward();

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
        std::string msText = std::to_string(delta * 1000) + " ms";
        DrawText(msText.c_str(), 10, 30, 20, DARKGREEN);

        std::string mouseSimInfo = std::to_string(mouseSimPos.x) + "," + std::to_string(mouseSimPos.y);
        DrawText(mouseSimInfo.c_str(), 10, 60, 20, BLUE);
        std::string smoothing = std::to_string(sim.densityKernel(50, 0));
        DrawText(smoothing.c_str(), 10, 80, 20, BLUE);
        //std::string density = std::to_string(sim.calculateDensity(mouseSimPos));
        //DrawText(density.c_str(), 10, 110, 20, BLUE);
        //Vector2 gradientVec = sim.calculateGradientVec(mouseSimPos);
        //std::string gradient = std::to_string(gradientVec.x) + "," + std::to_string(gradientVec.y);
        //DrawText(gradient.c_str(), 10, 130, 20, BLUE);

        EndDrawing();
    }
}
