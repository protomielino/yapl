#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>

#include <raylib.h>
#define RAYMATH_IMPLEMENTATION
#include <raymath.h>

#ifndef STB_DS_IMPLEMENTATION
#define STB_DS_IMPLEMENTATION
#endif
#include "stb_ds.h"

#include "maths.h"

#include "kd.h"
#include "sim.h"

const int WIDTH  = 700;
const int HEIGHT = 700;

Vector2 mouse_screen;
Vector2 mouse_world;

Vector2 offset;
Vector2 scale;

Vector2 startPan;

int bx, tx, dx;
int by, ty, dy;
int textDim;
char text[1024] = {};

size_t frameCounter;

bool show_help;
bool show_HUD;
bool show_tree;

void showHUD()
{
    if(show_HUD) {
        sprintf(text, "counter: %ld", frameCounter);
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "offset: {%d, %d} [CLICK+DRAG]", (int)offset.x, (int)offset.y);
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "scale: {%.2f, %.2f} [WHEEL]", scale.x, scale.y);
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        ty += dy;
    }
}

void showHelp()
{
    if(show_help) {
        sprintf(text, "[R] reset pan & zoom");
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "[MMB] reset pan");
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "[d] show/hide tree");
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        ty += dy;
    }
}

void processInputs()
{
    if(IsKeyDown(KEY_R) && IsKeyDown(KEY_LEFT_SHIFT)) {
        offset = Vector2Zero(); //(Vector2){-WIDTH/2.0f, -HEIGHT/2.0f};
        offset = Vector2Divide(offset, scale);
        scale = Vector2One();
    }
    if (IsMouseButtonPressed(MOUSE_BUTTON_MIDDLE)) {
        offset = Vector2Zero(); //(Vector2){-WIDTH/2.0f, -HEIGHT/2.0f};
        offset = Vector2Divide(offset, scale);
    }
    if(IsKeyPressed(KEY_D)) {
        show_tree = !show_tree;
    }
    if(IsKeyPressed(KEY_H) && IsKeyDown(KEY_LEFT_SHIFT)) {
        show_HUD = !show_HUD;
    }
    if(IsKeyPressed(KEY_H) && !IsKeyDown(KEY_LEFT_SHIFT)) {
        show_help = !show_help;
    }
}

void updatePanZoom()
{
    //----------------------------------------------------------------------------------
    // pan-zoom
    //----------------------------------------------------------------------------------
    // Just grab a copy of mouse coordinates for convenience
    mouse_screen = (Vector2){ (float)GetMouseX(), (float)GetMouseY() };
    mouse_world = ScreenToWorld(mouse_screen);

    // For panning, we need to capture the screen location when the user starts
    // to pan...
    if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
        startPan = mouse_screen;
    }

    // ...as the mouse moves, the screen location changes. Convert this screen
    // coordinate change into world coordinates to implement the pan. Simples.
    if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
        offset = Vector2Subtract(
                offset,
                Vector2Divide(
                        Vector2Subtract(mouse_screen, startPan),
                                scale));

        // Start "new" pan for next epoch
        startPan = mouse_screen;
    }

    // For zoom, we need to extract the location of the cursor before and after the
    // scale is changed. Here we get the cursor and translate into world space...
    Vector2 mouseWorld_BeforeZoom = {};
    mouseWorld_BeforeZoom = ScreenToWorld(mouse_screen);

    // ...change the scale as required...
    if (GetMouseWheelMove() > 0) {
        scale.x *= 1.1f;
        scale.y *= 1.1f;
    } else if (GetMouseWheelMove() < 0) {
        scale.x *= 0.9f;
        scale.y *= 0.9f;
    }

    // ...now get the location of the cursor in world space again - It will have changed
    // because the scale has changed, but we can offset our world now to fix the zoom
    // location in screen space, because we know how much it changed laterally between
    // the two spatial scales. Neat huh? ;-)
    Vector2 mouseWorld_AfterZoom = {};;
    mouseWorld_AfterZoom = ScreenToWorld(mouse_screen);
    offset = Vector2Add(
            offset,
            Vector2Subtract(
                    mouseWorld_BeforeZoom,
                    mouseWorld_AfterZoom));
}

void initPanZoom(Vector2 pan, Vector2 sca)
{
    startPan = pan;

    offset = Vector2Zero(); //(Vector2){-WIDTH/2.0f, -HEIGHT/2.0f};
    scale = sca;
}

void initHUD()
{
    textDim = 20;

    by = 50;
    ty = by;
    dy = textDim;
    bx = 10;
    tx = bx;
    dx = textDim;
}

// -------------------- configuration (match JS defaults) --------------------
const int DEFAULT_N = 1000;
const float DEFAULT_DT = 0.02f;
const float DEFAULT_FRICTION_HALF_LIFE = 0.040f;
const float DEFAULT_RMAX = 0.25f;
const int DEFAULT_M = 3;
const float FORCE_FACTOR = 10.0f;

// add periodic world size (positions are in [0,1] space; world wraps at 0/1)
const int USE_PERIODIC = 1; // 1 = wrap (periodic), 0 = no wrap (original clamp behavior)
const float WORLD_SIZE = 1.0f; // domain size in both x and y (positions remain in [0,WORLD_SIZE))

int KD_CURRENT_POINT_COUNT = 0;

// pointer set by sim_update before queries to allow kd_query_radius to access mark buffer
sim *CURRENT_SIM_INSTANCE = NULL;

// -------------------- drawing API --------------------
void draw_line_world(float x1, float y1, float x2, float y2)
{
    Vector2 v1 = { x1 * WIDTH, y1 * HEIGHT };
    Vector2 v2 = { x2 * WIDTH, y2 * HEIGHT };
    v1 = WorldToScreen(v1);
    v2 = WorldToScreen(v2);
    Color c = WHITE; c.a = 32;
    DrawLine(v1.x, v1.y, v2.x, v2.y, c);
}
void draw_particle_world(float x, float y, float size, int color)
{
    Color c = ColorFromHSV(360.0f / DEFAULT_M * color, 1.0f, 1.0f);
    Vector2 v = { x * WIDTH, y * HEIGHT };
    v = WorldToScreen(v);
    DrawPixelV(v, c);
//    DrawCircleV(v, size, c);
}

unsigned int rnd_state;

int main(int argc, char *argv[])
{
//    srand48(clock());
    srand((unsigned)clock());
    frameCounter = 0;

    // Initialization

    sim* s = sim_create(DEFAULT_N, DEFAULT_M, DEFAULT_DT, DEFAULT_FRICTION_HALF_LIFE, DEFAULT_RMAX, FORCE_FACTOR);

    // Initialise offset so 0,0 top left of the screen
    initPanZoom(Vector2Zero(), Vector2One());
    initHUD();

    //---------------------------------------------------------------------------------------

    InitWindow(WIDTH, HEIGHT, "Particle Life");
    //SetWindowState(FLAG_FULLSCREEN_MODE);
    SetTargetFPS(30);
    //--------------------------------------------------------------------------------------

    // Main game loop
    while (!WindowShouldClose())    // Detect window close button or ESC key
    {
        // Update
        //----------------------------------------------------------------------------------

        updatePanZoom();
        processInputs();

        frameCounter++;
        //----------------------------------------------------------------------------------

        // Draw
        //----------------------------------------------------------------------------------
        BeginDrawing(); {
            ClearBackground(BLACK);
//            draw_clear_background();

            sim_update(s);
            sim_draw_frame(s);

            Vector2 cpos = Vector2Zero();
            DrawLineV(
                    WorldToScreen((Vector2){cpos.x, 0.0f}),
                    WorldToScreen((Vector2){cpos.x, HEIGHT/2}), GREEN);
            DrawLineV(
                    WorldToScreen((Vector2){0.0f, cpos.y}),
                    WorldToScreen((Vector2){WIDTH/2, cpos.y}), RED);
            DrawCircleLinesV(WorldToScreen(cpos), 2.5*Vector2Length(scale), WHITE);

            sprintf(text, "[0, 0]");
            cpos = WorldToScreen(cpos);
            DrawText(text, cpos.x+5, cpos.y+5, textDim/3*Vector2Length(scale), WHITE);

            ty = by;
            showHUD();
            showHelp();
            sprintf(text, "help, HUD: [h,H]");
            DrawText(text, 10, ty, textDim, WHITE);
            ty += dy;

            DrawFPS(20, 20);
        } EndDrawing();
        //----------------------------------------------------------------------------------
    }

    // De-Initialization
    sim_free(s);
    //--------------------------------------------------------------------------------------
    CloseWindow();        // Close window and OpenGL context
    //--------------------------------------------------------------------------------------

    return 0;
}
