#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>

#include <raylib.h>
#define RAYMATH_IMPLEMENTATION
#include <raymath.h>

#define STB_DS_IMPLEMENTATION
#include "stb_ds.h"

#include "maths.h"

#define WIDTH  1360
#define HEIGHT 768

Vector2 mouse_scren;
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
bool pause;
bool display;

int n;
float dt;
float frictionHalfLife;
float rMax;
int m;
float **matrix;
float forceFactor;

float frictionFactor;

float real_time = 0.0f;


void randomizeParticles();
float** makeRandomMatrix();

void showHUD()
{
    if(show_HUD) {
        sprintf(text, "spf: {%.4f} [s]", GetFrameTime());
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "counter: %ld", frameCounter);
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "offset: {%d, %d} [CLICK+DRAG]", (int)offset.x, (int)offset.y);
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "scale: {%.2f, %.2f} [WHEEL]", scale.x, scale.y);
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "dt: {%.4f} [ms]", dt);
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "t: {%.4f} [s]", real_time);
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;

        ty += dy;
    }
}

void processInputs()
{
    if(IsKeyDown(KEY_R) && IsKeyDown(KEY_LEFT_SHIFT)) {
        offset = Vector2Zero();
        scale = Vector2One();
    }
    if(IsKeyPressed(KEY_R) && !IsKeyDown(KEY_LEFT_SHIFT)) {
        randomizeParticles();
    }
    if (IsMouseButtonPressed(MOUSE_BUTTON_MIDDLE)) {
        offset = (Vector2){ 0, 0 };
//            scale = Vector2One();
    }
    if(IsKeyPressed(KEY_H) && IsKeyDown(KEY_LEFT_SHIFT)) {
        show_HUD = !show_HUD;
    }
    if(IsKeyPressed(KEY_SPACE) && !IsKeyDown(KEY_LEFT_SHIFT)) {
        matrix = makeRandomMatrix();
//        dt = 0.02f;
    }
    if(IsKeyPressed(KEY_P)) {
        pause = !pause;
    }
    if(IsKeyPressed(KEY_D)) {
        display = !display;
    }
    if(IsKeyPressed(KEY_T) && IsKeyDown(KEY_LEFT_SHIFT)) {
        dt -= 0.001f;
    }
    if(IsKeyPressed(KEY_T) && !IsKeyDown(KEY_LEFT_SHIFT)) {
        dt += 0.001f;
    }
    if(IsKeyPressed(KEY_Z)) {
        dt = 0.02f;
    }
}

void updatePanZoom()
{
    //----------------------------------------------------------------------------------
    // pan-zoom
    //----------------------------------------------------------------------------------
    // Just grab a copy of mouse coordinates for convenience
    mouse_scren = (Vector2){ (float)GetMouseX(), (float)GetMouseY() };
    mouse_world = ScreenToWorld(mouse_scren);

    // For panning, we need to capture the screen location when the user starts
    // to pan...
    if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
        startPan = mouse_scren;
    }

    // ...as the mouse moves, the screen location changes. Convert this screen
    // coordinate change into world coordinates to implement the pan. Simples.
    if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
        offset = Vector2Subtract(
                offset,
                Vector2Divide(
                        Vector2Subtract(
                                mouse_scren,
                                startPan),
                                scale));

        // Start "new" pan for next epoch
        startPan = mouse_scren;
    }

    // For zoom, we need to extract the location of the cursor before and after the
    // scale is changed. Here we get the cursor and translate into world space...
    Vector2 mouseWorld_BeforeZoom = {};
    mouseWorld_BeforeZoom = ScreenToWorld(mouse_scren);

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
    mouseWorld_AfterZoom = ScreenToWorld(mouse_scren);
    offset = Vector2Add(
            offset,
            Vector2Subtract(
                    mouseWorld_BeforeZoom,
                    mouseWorld_AfterZoom));
}

void initPanZoom(Vector2 pan, Vector2 sca)
{
    startPan = Vector2Zero();

    // Initialise offset so 0,0 top left of the screen
    offset = pan;
    scale = sca;
}

void iniHUD()
{
    textDim = 20;

    by = 50;
    ty = by;
    dy = textDim;
    bx = 10;
    tx = bx;
    dx = textDim;
}

float** makeRandomMatrix()
{
    if (arrlen(matrix) != 0) {
        for (int i = 0; i < arrlen(matrix); i++) {
            arrfree(matrix[i]);
        }
    }
    arrfree(matrix);

    float **rows = NULL;
    arrsetlen(rows, m);
    for (int i = 0; i < m; i++) {
        rows[i] = NULL;
        arrsetlen(rows[i], m);
        for (int j = 0; j < m; j++) {
            rows[i][j] = drand48() * 2.0f - 1.0f;
        }
    }
    return rows;
}

unsigned int *colors = NULL;
float *positionsX = NULL;
float *positionsY  = NULL;
float *velocitiesX = NULL;
float *velocitiesY = NULL;

void randomizeParticles()
{
    for(int i = 0; i < n; i++) {
        colors[i] = floorf(drand48() * m);
        positionsX[i] = drand48()*WIDTH;
        positionsY[i] = drand48()*HEIGHT;
        velocitiesX[i] = 0;
        velocitiesY[i] = 0;
    }
}

float force(float r, float a)
{
    float beta = 0.3;
    if (r < beta) {
        return r / beta - 1;
    } else if (beta < r && r < 1) {
        return a * (1 - fabsf(2 * r - 1 - beta) / (1 - beta));
    } else {
        return 0;
    }
}

void updateParticles()
{
    float totalVelocity_len_max = 0.0f;
    float totalVelocity_len_min = FLT_MAX;
    float totalVelocity_len = 0.0f;
    float DT = 0.0f;

    float totalForceX = 0.0f;
    float totalForceY = 0.0f;
    float totalVelocityX = 0.0f;
    float totalVelocityY = 0.0f;

    for (int i = 0; i < n; i++) {
        totalForceX = 0.0f;
        totalForceY = 0.0f;
        totalVelocityX = 0.0f;
        totalVelocityY = 0.0f;

        for (int j = 0; j < n; j++) {
            if (j == i) continue;
            float rx = positionsX[j] - positionsX[i];
            float ry = positionsY[j] - positionsY[i];
            if(rx > 0.5f * WIDTH) {
                rx -= WIDTH;
            }
            if(rx < -0.5f * WIDTH) {
                rx += WIDTH;
            }
            if(ry > 0.5f * HEIGHT) {
                ry -= HEIGHT;
            }
            if(ry < -0.5f * HEIGHT) {
                ry += HEIGHT;
            }

            float r = sqrt(rx*rx + ry*ry);
            if (r <= rMax) {
                if (r > 0 && r < rMax) {
                    float f = force(r / rMax, matrix[colors[i]][colors[j]]);
                    totalForceX += rx / r * f;
                    totalForceY += ry / r * f;
                }
            }
        }

        totalForceX *= rMax * forceFactor;
        totalForceY *= rMax * forceFactor;

        velocitiesX[i] *= frictionFactor;
        velocitiesY[i] *= frictionFactor;

        velocitiesX[i] += totalForceX * dt;
        velocitiesY[i] += totalForceY * dt;

        totalVelocityX += velocitiesX[i];
        totalVelocityY += velocitiesY[i];

#if 0
        Vector2 v = {totalVelocityX, totalVelocityY};
        totalVelocity_len = Vector2Length(v);
        totalVelocity_len_max =
                totalVelocity_len_max < totalVelocity_len ?
                        totalVelocity_len :
                        totalVelocity_len_max;
        totalVelocity_len_min =
                totalVelocity_len_min > totalVelocity_len ?
                        totalVelocity_len :
                        totalVelocity_len_min;
    }

    float min_vel = 0.0f;
    float max_vel = 750.0f;

    DT = map(totalVelocity_len, min_vel, max_vel, 1.0f, 0.1f);

    if (DT == 0.0f) {
        DT = 0.999;
    }

    if (totalVelocity_len_max < max_vel && totalVelocity_len_min > min_vel) {
        if (dt > 0.02f) {
            dt *= 0.95;
            dt *= DT;
        }
        if (dt < 0.02f) {
            dt *= 1.05;
            dt /= DT;
        }
    }
    if (totalVelocity_len_max > max_vel) {
        dt /= DT;
        dt *= 0.999;
        if (dt < 0.005f) {
            dt = 0.005f;
        }
    }
#else
    }
#endif
    real_time += dt;


//    if (v_len_min < min_vel) {
//        dt *= 1.005f;
//        if (dt > 2.0f) {
//            dt = 2.0f;
//        }
//    }

    for (int i = 0; i < n; i++) {
        positionsX[i] += velocitiesX[i] * dt;
        positionsY[i] += velocitiesY[i] * dt;

        // wrap around - toroidal world
        positionsX[i] = fmodf(positionsX[i] + (float)WIDTH, (float)WIDTH);
        positionsY[i] = fmodf(positionsY[i] + (float)HEIGHT, (float)HEIGHT);
    }
}

void drawParticles()
{
    //ClearBackground(BLACK);

    for (int i = 0; i < n; i++) {
        float screenX = positionsX[i];
        float screenY = positionsY[i];
        Color col = ColorFromHSV(360.0f * ((float)colors[i] / (float)m), 1.0f, 1.0);
        Vector2 world = {screenX, screenY};
//        DrawPixelV(WorldToScreen(world), col);
        DrawCircleV(WorldToScreen(world), 2, col);
    }
}

int main(int argc, char *argv[])
{
//    srand48(clock());
    frameCounter = 0;

    // Initialization
    show_help = false;
    show_HUD = true;
    pause = false;
    display = true;

    n = 2000;
    dt = 0.02f;
    frictionHalfLife = 0.040f;
    rMax = 200.0f;
    m = 10;
    matrix = makeRandomMatrix();
    forceFactor = 1.0f;

    frictionFactor = powf(0.5, dt / frictionHalfLife);

    arrsetlen(colors, n);
    arrsetlen(positionsX, n);
    arrsetlen(positionsY, n);
    arrsetlen(velocitiesX, n);
    arrsetlen(velocitiesY, n);
    randomizeParticles();

    initPanZoom(Vector2Zero(), Vector2One());
    iniHUD();

    //---------------------------------------------------------------------------------------

    InitWindow(WIDTH, HEIGHT, "Particle Life");
    SetWindowState(FLAG_FULLSCREEN_MODE);
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

        if (!pause)
            updateParticles();

        //----------------------------------------------------------------------------------

        // Draw
        //----------------------------------------------------------------------------------
        BeginDrawing(); {
            ClearBackground(BLACK);

            Vector2 cpos = Vector2Zero();
            DrawLineV(
                    WorldToScreen((Vector2){cpos.x, 0}),
                    WorldToScreen((Vector2){cpos.x, HEIGHT}), GREEN);
            DrawLineV(
                    WorldToScreen((Vector2){0, cpos.y}),
                    WorldToScreen((Vector2){WIDTH, cpos.y}), RED);
            DrawLineV(
                    WorldToScreen((Vector2){WIDTH, 0}),
                    WorldToScreen((Vector2){WIDTH, HEIGHT}), WHITE);
            DrawLineV(
                    WorldToScreen((Vector2){0, HEIGHT}),
                    WorldToScreen((Vector2){WIDTH, HEIGHT}), WHITE);
            DrawCircleLinesV(WorldToScreen(cpos), 2.5*Vector2Length(scale), WHITE);

            if (display)
                drawParticles();

            sprintf(text, "[0, 0]");
            cpos = WorldToScreen(cpos);
            DrawText(text, cpos.x+5, cpos.y+5, textDim/3*Vector2Length(scale), WHITE);

            ty = by;
            showHUD();
            sprintf(text, "help, HUD: [h,H]");
            DrawText(text, 10, ty, textDim, WHITE);
            ty += dy;

            DrawFPS(20, 20);
        } EndDrawing();
        //----------------------------------------------------------------------------------
    }

    // De-Initialization
    //--------------------------------------------------------------------------------------
    CloseWindow();        // Close window and OpenGL context
    //--------------------------------------------------------------------------------------

    return 0;
}
