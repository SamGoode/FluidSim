#include "raylib.h"
#include "Simulation.h"
#include "Particle.h"

int main() {
    int screenWidth = 1600;
    int screenHeight = 800;

    InitWindow(screenWidth, screenHeight, "Smooth Particle Hydrodynamics Sim");

    SetTargetFPS(240);

    srand(78925311);

    Simulation sim({ 200, 100, (float)screenWidth - 200, (float)screenHeight - 100 });

    while (!WindowShouldClose()) {
        // Updates
        float delta = GetFrameTime();


        // Drawing
        BeginDrawing();

        ClearBackground(RAYWHITE);

        sim.draw();

        DrawFPS(10, 10);

        EndDrawing();
    }
}
