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

#define WIDTH  700
#define HEIGHT 768

int worldWidth;
int worldHeight;

Vector2 mouse_scren;
Vector2 mouse_world;

Vector2 offset;
Vector2 scale;

Vector2 startPan;

int bx, tx, dx;
int by, ty, dy;
float textDim;
char text[1024] = {};

size_t frameCounter;

bool pause;
bool display;
bool showHelp;
bool showHUD;
bool showForcesDBG;
bool showMinDistancesDBG;
bool showMassesDBG;
bool showRadiiDBG;

int numParticles;
int numTypes;
float rMax;
float dt;
float **forces;
float forceFactor;
float frictionHalfLife;

float hueStep;
float satStep;
float valStep;

float frictionFactor;

float real_time = 0.0f;

unsigned int *particleType = NULL;
float *positionsX = NULL;
float *positionsY  = NULL;
float *velocitiesX = NULL;
float *velocitiesY = NULL;

void randomizeParticles();
float** makeRandomMatrix();

float cell_size;
int cell_row_hover;
int cell_col_hover;
int cell_row_selected;
int cell_col_selected;

void shift_particles_up()
{
    for (int i = 0; i < numParticles; ++i) {
        positionsY[i] -= 10.0f;
    }
}
void shift_particles_down()
{
    for (int i = 0; i < numParticles; ++i) {
        positionsY[i] += 10.0f;
    }
}
void shift_particles_left()
{
    for (int i = 0; i < numParticles; ++i) {
        positionsX[i] -= 10.0f;
    }
}
void shift_particles_right()
{
    for (int i = 0; i < numParticles; ++i) {
        positionsX[i] += 10.0f;
    }
}

void zeroParticlesForces()
{
    for(int i = 0; i < numTypes; i++) {
        for(int j = 0; j < numTypes; j++) {
            forces[i][j] = 0.0f;
        }
    }
}

void show_HUD()
{
    if(showHUD) {
        sprintf(text, "spf: {%.4f} [s]", GetFrameTime());
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "offset: {%d, %d} [CLICK+DRAG]", (int)offset.x, (int)offset.y);
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "scale: {%.2f, %.2f} [WHEEL]", scale.x, scale.y);
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "numParticles: {%d}", numParticles);
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "dt: {%.4f} [ms]", dt);
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "t: {%.4f} [s]", real_time);
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "counter: %ld", frameCounter);
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;

        ty += dy;
    }
}

void show_DBG()
{
    float textDim_tmp = textDim;
    float dx_tmp = dx;
    float dy_tmp = dy;
//    textDim /= 1.5;
    dx = dx*3 + 2;
    //dy = dy/2;


    if(showForcesDBG) {
        for (int col = 0; col < numTypes; ++col) {
            Vector2 pos = {bx + cell_size + (cell_size+2)*col + 5, ty};
            Vector2 siz = {cell_size, cell_size};
            DrawRectangleV(pos, siz, ColorFromHSV((float)col*hueStep, 0.8f, 1.0f));
        }
        for (int row = 0; row < numTypes; ++row) {
            Vector2 pos = {bx, ty + cell_size + (cell_size+2)*row + 5};
            Vector2 siz = {cell_size, cell_size};
            DrawRectangleV(pos, siz, ColorFromHSV((float)row*hueStep, 0.8f, 1.0f));
        }
        for (int row = 0; row < numTypes; ++row) {
            for (int col = 0; col < numTypes; ++col) {
                Vector2 pos = {bx + cell_size + (cell_size+2)*col + 5, ty + cell_size + (cell_size+2)*row + 5};
                Vector2 siz = {cell_size, cell_size};
                Color c = {};
                if (forces[row][col] <= 0.0f) {
                    c = ColorFromHSV(360.0f, 1.0f, 1.0f - forces[row][col]);
                } else {
                    c = ColorFromHSV(360.0f / 3.0f, 1.0f, forces[row][col]);
                }
                if (forces[row][col] == 0.0f) {
                    c = BLACK;
                }
                if (forces[row][col] == -0.0f) {
                    c = BLACK;
                }
                DrawRectangleV(pos, siz, c);
                if (cell_row_hover == row && cell_col_hover == col) {
                    DrawRectangleLines(pos.x-1, pos.y-1, siz.x+2, siz.y+2, WHITE);
                }
                if (cell_row_selected == row && cell_col_selected == col) {
                    DrawRectangleLines(pos.x-2, pos.y-2, siz.x+4, siz.y+4, YELLOW);
                }

                sprintf(text, "%+.1f", forces[row][col]);
                DrawText(text, pos.x+cell_size/2.0f - textDim, pos.y + cell_size/2.0f - textDim/3, textDim+3, BLACK);
                DrawText(text, 2+ pos.x+cell_size/2.0f - textDim, 1+ pos.y + cell_size/2.0f - textDim/3, textDim, WHITE);
            }
        }
#if 0
        float col_tdx = dx;
        float col_tdy = dy;
        float col_tbx = bx;
        float col_tby = ty;
        float col_tx = col_tbx + col_tdx/2;
        float col_ty = col_tby + col_tdy*2.0f + textDim/2.0f;

        for (int row = 0; row < numTypes; ++row) {
            DrawCircle(col_tx, col_ty, 5,ColorFromHSV((float)row*hueStep, 0.8f, 1.0f));
            col_tx += col_tdx;
        }
        col_tx = bx;
        for (int col = 0; col < numTypes; ++col) {
            DrawCircle(col_tx, col_ty+col_tdy, 5,ColorFromHSV((float)col*hueStep, 0.8f, 1.0f));
            col_ty += col_tdy;
        }

        ty += dy;

        sprintf(text, "forces:");
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        ty += dy;
        for (int i = 0; i < numTypes; i++) {
            tx = bx+dx/4;

            for (int j = 0; j < numTypes; j++) {
                sprintf(text, "[       ]");
                DrawText(text, tx, ty, textDim,
                            ColorFromHSV((float)i*hueStep,
                                         0.8f,
                                         1.0f));
                sprintf(text, "  %+.1f  ", forces[i][j]);
                DrawText(text, tx, ty, textDim,
                            ColorFromHSV((float)j*hueStep,
                                         0.8f,
                                         1.0f));
                tx += dx;
            }
            ty += dy;
        }
        ty += dy;
#endif
    }
    if(showMinDistancesDBG) {
//        sprintf(text, "minDistances:");
//        DrawText(text, 10, ty, textDim, WHITE);
//        ty += dy;
//        for (int i = 0; i < numTypes; i++) {
//            tx = bx;
//            for (int j = 0; j < numTypes; j++) {
//                sprintf(text, "[         ]");
//                DrawText(text, tx, ty, textDim,
//                            ColorFromHSV((float)i*hueStep,
//                                         0.8f,
//                                         1.0f));
//                sprintf(text, "  %+.1f  ", minDistances[i][j]);
//                DrawText(text, tx, ty, textDim,
//                            ColorFromHSV((float)j*hueStep,
//                                         0.8f,
//                                         1.0f));
//                tx += dx;
//            }
//            ty += dy;
//        }
//        ty += dy;
    }
    if(showMassesDBG) {
//        sprintf(text, "masses:");
//        DrawText(text, 10, ty, textDim, WHITE);
//        ty += dy;
//        tx = bx;
//        for (int i = 0; i < numTypes; i++) {
//            sprintf(text, "[       ]");
//            DrawText(text, tx, ty, textDim,
//                        ColorFromHSV((float)i*hueStep,
//                                     0.8f,
//                                     1.0f));
//            sprintf(text, "  %+.1f  ", masses[i]);
//            DrawText(text, tx, ty, textDim,
//                        ColorFromHSV((float)i*hueStep,
//                                     0.8f,
//                                     1.0f));
//            tx += dx;
//        }
//        ty += dy;
//        ty += dy;
    }
    if(showRadiiDBG) {
//        sprintf(text, "radii:");
//        DrawText(text, 10, ty, textDim, WHITE);
//        ty += dy;
//        for (int i = 0; i < numTypes; i++) {
//            tx = bx;
//            for (int j = 0; j < numTypes; j++) {
//                sprintf(text, "[          ]");
//                DrawText(text, tx, ty, textDim,
//                            ColorFromHSV((float)i*hueStep,
//                                         0.8f,
//                                         1.0f));
//                sprintf(text, "  %+.1f  ", radii[i][j]);
//                DrawText(text, tx, ty, textDim,
//                            ColorFromHSV((float)j*hueStep,
//                                         0.8f,
//                                         1.0f));
//                tx += dx;
//            }
//            ty += dy;
//        }
//        ty += dy;
    }

    dx = dx_tmp;
    dy = dy_tmp;
    textDim = textDim_tmp;

    if(showHelp) {
        tx = bx;
        sprintf(text, "[r] : randomize positions");
        DrawText(text, tx, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "[R] : reset pan-zoom");
        DrawText(text, tx, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "[SPACE] : randomize parameters");
        DrawText(text, tx, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "[z,x,c,v] : forces, minDistances, masses, radii");
        DrawText(text, tx, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "[LSHIFT+DIR_KEYS] : shift field");
        DrawText(text, tx, ty, textDim, WHITE);
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
    if(IsKeyPressed(KEY_RIGHT_CONTROL)) {
        zeroParticlesForces();
    }
    if (IsMouseButtonPressed(MOUSE_BUTTON_MIDDLE)) {
        offset = (Vector2){ 0, 0 };
//            scale = Vector2One();
    }
    if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT) && IsKeyDown(KEY_LEFT_CONTROL)) {
            cell_row_selected = cell_row_hover;
            cell_col_selected = cell_col_hover;
    }
    if(IsKeyPressed(KEY_H) && IsKeyDown(KEY_LEFT_SHIFT)) {
        showHUD = !showHUD;
    }
    if(IsKeyPressed(KEY_H) && !IsKeyDown(KEY_LEFT_SHIFT)) {
        showHelp = !showHelp;
    }
    if(IsKeyPressed(KEY_SPACE) && !IsKeyDown(KEY_LEFT_SHIFT)) {
        forces = makeRandomMatrix();
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
        showForcesDBG = !showForcesDBG;
        showMinDistancesDBG = false;
        showMassesDBG = false;
        showRadiiDBG = false;
    }
    if(IsKeyPressed(KEY_X)) {
        showForcesDBG = false;
        showMinDistancesDBG = !showMinDistancesDBG;
        showMassesDBG = false;
        showRadiiDBG = false;
    }
    if(IsKeyPressed(KEY_C)) {
        showForcesDBG = false;
        showMinDistancesDBG = false;
        showMassesDBG = !showMassesDBG;
        showRadiiDBG = false;
    }
    if(IsKeyPressed(KEY_V)) {
        showForcesDBG = false;
        showMinDistancesDBG = false;
        showMassesDBG = false;
        showRadiiDBG = !showRadiiDBG;
    }
    if(IsKeyDown(KEY_UP) && IsKeyDown(KEY_LEFT_SHIFT)) {
        shift_particles_down();
    }
    if(IsKeyDown(KEY_DOWN) && IsKeyDown(KEY_LEFT_SHIFT)) {
        shift_particles_up();
    }
    if(IsKeyDown(KEY_LEFT) && IsKeyDown(KEY_LEFT_SHIFT)) {
        shift_particles_right();
    }
    if(IsKeyDown(KEY_RIGHT) && IsKeyDown(KEY_LEFT_SHIFT)) {
        shift_particles_left();
    }
    if(IsKeyPressed(KEY_UP)) {
        forces[cell_row_hover][cell_col_hover] += 0.1f;
    }
    if(IsKeyPressed(KEY_DOWN)) {
        forces[cell_row_hover][cell_col_hover] -= 0.1f;
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
    if (arrlen(forces) != 0) {
        for (int i = 0; i < arrlen(forces); i++) {
            arrfree(forces[i]);
        }
    }
    arrfree(forces);

    float **rows = NULL;
    arrsetlen(rows, numTypes);
    for (int i = 0; i < numTypes; i++) {
        rows[i] = NULL;
        arrsetlen(rows[i], numTypes);
        for (int j = 0; j < numTypes; j++) {
            rows[i][j] = drand48() * 2.0f - 1.0f;
        }
    }
    return rows;
}

void randomizeParticles()
{
    for(int i = 0; i < numParticles; i++) {
        particleType[i] = floorf(drand48() * numTypes);
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

#undef FOECR_EASE // experimental
void updateParticles()
{
#ifdef FOECR_EASE
    float totalVelocity_len_max = 0.0f;
    float totalVelocity_len_min = FLT_MAX;
    float totalVelocity_len = 0.0f;
    float DT = 0.0f;
#endif /* FOECR_EASE */

    float totalForceX = 0.0f;
    float totalForceY = 0.0f;
    float totalVelocityX = 0.0f;
    float totalVelocityY = 0.0f;

    for (int i = 0; i < numParticles; i++) {
        totalForceX = 0.0f;
        totalForceY = 0.0f;
        totalVelocityX = 0.0f;
        totalVelocityY = 0.0f;

        for (int j = 0; j < numParticles; j++) {
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
                    float f = force(r / rMax, forces[particleType[i]][particleType[j]]);
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

#ifdef FOECR_EASE
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
//    if (v_len_min < min_vel) {
//        dt *= 1.005f;
//        if (dt > 2.0f) {
//            dt = 2.0f;
//        }
//    }
#else
    }
#endif /* FOECR_EASE */

    for (int i = 0; i < numParticles; i++) {
        positionsX[i] += velocitiesX[i] * dt;
        positionsY[i] += velocitiesY[i] * dt;

        // wrap around - toroidal world
        positionsX[i] = fmodf(positionsX[i] + (float)WIDTH, (float)WIDTH);
        positionsY[i] = fmodf(positionsY[i] + (float)HEIGHT, (float)HEIGHT);
    }

    real_time += dt;
}

void drawParticles()
{
    //ClearBackground(BLACK);

    for (int i = 0; i < numParticles; i++) {
        float screenX = positionsX[i];
        float screenY = positionsY[i];
//        Color col = ColorFromHSV(360.0f * ((float)particleType[i] / (float)numTypes), 1.0f, 1.0);
        Color col = ColorFromHSV(particleType[i]*hueStep, 1.0f, 1.0);
        Vector2 world = {screenX, screenY};
//        DrawPixelV(WorldToScreen(world), col);
        float r = 2.0f * Vector2Length(scale);
        r = r > 1.0f ? r : 1.0f;
        DrawCircleV(WorldToScreen(world), r, col);
    }
}

int main(int argc, char *argv[])
{
    srand48(clock());
    frameCounter = 0;

    // Initialization
    showHelp = false;
    showHUD = true;
    showForcesDBG = false;
    showMinDistancesDBG = false;
    showMassesDBG = false;
    showRadiiDBG = false;

    pause = false;
    display = true;

    numParticles = 1000;
    numTypes = 20;
    dt = 0.02f;
    frictionHalfLife = 0.040f;
    rMax = 150.0f;
    forces = makeRandomMatrix();
    forceFactor = 1.0f;

    frictionFactor = powf(0.5, dt / frictionHalfLife);

    worldWidth = WIDTH;
    worldHeight = HEIGHT;
    scale = Vector2One();

    startPan = (Vector2){ 0.0f, 0.0f };

    hueStep = 360.0f / (float)numTypes;
    satStep = 1.0f / (float)numTypes;
    valStep = 1.0f / (float)numTypes;

    textDim = 20;

    by = 50;
    ty = by;
    dy = textDim;
    bx = 10;
    tx = bx;
    dx = textDim;

    cell_size = 50.0f;
    cell_row_hover = 0;
    cell_col_hover = 0;
    cell_row_selected = -1;
    cell_col_selected = -1;


    arrsetlen(particleType, numParticles);
    arrsetlen(positionsX, numParticles);
    arrsetlen(positionsY, numParticles);
    arrsetlen(velocitiesX, numParticles);
    arrsetlen(velocitiesY, numParticles);
    randomizeParticles();

    initPanZoom(Vector2Zero(), Vector2One());
    iniHUD();

    //---------------------------------------------------------------------------------------

    InitWindow(WIDTH, HEIGHT, "Particle Life");
//    SetWindowState(FLAG_FULLSCREEN_MODE);
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

        if (showForcesDBG) {
            Vector2 mouse_pos = GetMousePosition();

            cell_row_hover = (mouse_pos.y - ty - cell_size - 5) / (cell_size+2);
            cell_col_hover = (mouse_pos.x - bx - cell_size - 5) / (cell_size+2);

            cell_row_hover =
                    cell_row_hover < 0 ? 0 :
                            cell_row_hover > numTypes ? numTypes :
                                    cell_row_hover;
            cell_col_hover =
                    cell_col_hover < 0 ? 0 :
                            cell_col_hover > numTypes ? numTypes :
                                    cell_col_hover;
        }

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
            show_HUD();
            sprintf(text, "help, HUD: [h,H]");
            DrawText(text, 10, ty, textDim, WHITE);
            ty += dy;

            show_DBG();

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
