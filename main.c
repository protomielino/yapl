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

#include "particle.h"

#define WIDTH  700
#define HEIGHT 700

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
float minRadius;
float maxRadius;
float dt;
float **forces;
float **minDistances;
float *masses;
float **radii;
float forceFactor;
float frictionHalfLife;

float minParticleMass, maxParticleMass;

particle *swarm;
int swarm_len;

float hueStep;
float satStep;
float valStep;

float frictionFactor;

float real_time = 0.0f;

particle *swarm = NULL;

void randomizeParticles();
float** makeRandomMatrix();

float cell_size;
int cell_row_hover;
int cell_col_hover;
int cell_row_selected;
int cell_col_selected;

void shift_particles_up(particle *this)
{
    for (int i = 0; i < numParticles; ++i) {
        this[i].position.y -= 10.0f;
    }
}
void shift_particles_down(particle *this)
{
    for (int i = 0; i < numParticles; ++i) {
        this[i].position.y += 10.0f;
    }
}
void shift_particles_left(particle *this)
{
    for (int i = 0; i < numParticles; ++i) {
        this[i].position.x -= 10.0f;
    }
}
void shift_particles_right(particle *this)
{
    for (int i = 0; i < numParticles; ++i) {
        this[i].position.x += 10.0f;
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

void zeroParticlesMasses()
{
    for(int i = 0; i < numTypes; i++) {
            masses[i] = 0.0f;
    }
}

void zeroParticlesMinDistances()
{
    for(int i = 0; i < numTypes; i++) {
        for(int j = 0; j < numTypes; j++) {
            minDistances[i][j] = 0.0f;
        }
    }
}
void zeroParticlesRadii()
{
    for(int i = 0; i < numTypes; i++) {
        for(int j = 0; j < numTypes; j++) {
            radii[i][j] = 0.0f;
        }
    }
}

void randomizeMasses()
{
    for(int i = 0; i < numTypes; i++) {
        masses[i] = random_range(minParticleMass, maxParticleMass);
    }
}

void randomizeMinDistances()
{
    for (int row = 0; row < numTypes; ++row) {
        for (int col = 0; col < numTypes; ++col) {
            minDistances[row][col] = map(random_range(0.0f, minRadius), 0.0f, minRadius, 0.0f, minRadius/maxRadius);
        }
    }
}

void randomizeParticles(particle *this)
{
    for(int i = 0; i < numParticles; i++) {
        this[i].type = floorf(drand48() * numTypes);
        this[i].position = (Vector2){drand48()*WIDTH, drand48()*HEIGHT};
        this[i].velocity = Vector2Zero();
        this[i].mass = masses[this[i].type];
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
    dx = dx*3 + 2;

    if(showHelp) {
        ty += dy;

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

    if(showForcesDBG) {
        ty += dy;

        sprintf(text, "forces:");
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "min: %+.2f", -1.0f);
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "max: %+.2f", 1.0f);
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        ty += dy;


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
                c.a = 128;
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
    }

    if(showMinDistancesDBG) {
        ty += dy;

        sprintf(text, "minDistances:");
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "min: %.2f", 0.0f);
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "max: %.2f", minRadius);
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        ty += dy;

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
                if (minDistances[row][col] <= 0.0f) {
                    c = ColorFromHSV(360.0f, 1.0f, 1.0f - minDistances[row][col]);
                } else {
                    c = ColorFromHSV(360.0f / 3.0f, 1.0f, minDistances[row][col]);
                }
                if (minDistances[row][col] == 0.0f) {
                    c = BLACK;
                }
                if (minDistances[row][col] == -0.0f) {
                    c = BLACK;
                }
                c.a = 128;
                DrawRectangleV(pos, siz, c);
                if (cell_row_hover == row && cell_col_hover == col) {
                    DrawRectangleLines(pos.x-1, pos.y-1, siz.x+2, siz.y+2, WHITE);
                }
                if (cell_row_selected == row && cell_col_selected == col) {
                    DrawRectangleLines(pos.x-2, pos.y-2, siz.x+4, siz.y+4, YELLOW);
                }

                sprintf(text, "%+.1f", minDistances[row][col]);
                DrawText(text, pos.x+cell_size/2.0f - textDim, pos.y + cell_size/2.0f - textDim/3, textDim+3, BLACK);
                DrawText(text, 2+ pos.x+cell_size/2.0f - textDim, 1+ pos.y + cell_size/2.0f - textDim/3, textDim, WHITE);
            }
        }
    }

    if(showMassesDBG) {
        ty += dy;

        sprintf(text, "masses:");
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "min: %.1f", minParticleMass);
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "max: %.1f", maxParticleMass);
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        ty += dy;

        for (int col = 0; col < numTypes; ++col) {
            Vector2 pos = {bx + cell_size + (cell_size+2)*col + 5, ty};
            Vector2 siz = {cell_size, cell_size};
            DrawRectangleV(pos, siz, ColorFromHSV((float)col*hueStep, 0.8f, 1.0f));
        }
        for (int col = 0; col < numTypes; ++col) {
            int row = 0;
            Vector2 pos = {bx + cell_size + (cell_size+2)*col + 5, ty + cell_size + (cell_size+2)*row + 5};
            Vector2 siz = {cell_size, cell_size};
            Color c = {};
            c = ColorFromHSV(360.0f / 3.0f, 1.0f, map(masses[col], minParticleMass, maxParticleMass, 0.0f, 1.0f));
            c.a = 128;
            DrawRectangleV(pos, siz, c);
            if (cell_row_hover == row && cell_col_hover == col) {
                DrawRectangleLines(pos.x-1, pos.y-1, siz.x+2, siz.y+2, WHITE);
            }
            if (cell_row_selected == row && cell_col_selected == col) {
                DrawRectangleLines(pos.x-2, pos.y-2, siz.x+4, siz.y+4, YELLOW);
            }

            sprintf(text, " %.1f", masses[col]);
            DrawText(text, pos.x+cell_size/2.0f - textDim, pos.y + cell_size/2.0f - textDim/3, textDim+3, BLACK);
            DrawText(text, 2+ pos.x+cell_size/2.0f - textDim, 1+ pos.y + cell_size/2.0f - textDim/3, textDim, WHITE);
        }
    }

    if(showRadiiDBG) {
        ty += dy;

        sprintf(text, "radii:");
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "min: %.2f", 0.0f);
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "max: %.2f", 0.3f);
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        ty += dy;

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
                if (radii[row][col] <= 0.0f) {
                    c = ColorFromHSV(360.0f, 1.0f, 1.0f - radii[row][col]);
                } else {
                    c = ColorFromHSV(360.0f / 3.0f, 1.0f, radii[row][col]);
                }
                if (radii[row][col] == 0.0f) {
                    c = BLACK;
                }
                if (radii[row][col] == -0.0f) {
                    c = BLACK;
                }
                c.a = 128;
                DrawRectangleV(pos, siz, c);
                if (cell_row_hover == row && cell_col_hover == col) {
                    DrawRectangleLines(pos.x-1, pos.y-1, siz.x+2, siz.y+2, WHITE);
                }
                if (cell_row_selected == row && cell_col_selected == col) {
                    DrawRectangleLines(pos.x-2, pos.y-2, siz.x+4, siz.y+4, YELLOW);
                }

                sprintf(text, "%+.1f", radii[row][col]);
                DrawText(text, pos.x+cell_size/2.0f - textDim, pos.y + cell_size/2.0f - textDim/3, textDim+3, BLACK);
                DrawText(text, 2+ pos.x+cell_size/2.0f - textDim, 1+ pos.y + cell_size/2.0f - textDim/3, textDim, WHITE);
            }
        }
    }
}

void processInputs()
{
    if(IsKeyDown(KEY_R) && IsKeyDown(KEY_LEFT_SHIFT)) {
        offset = Vector2Zero();
        scale = Vector2One();
    }
    if(IsKeyPressed(KEY_R) && !IsKeyDown(KEY_LEFT_SHIFT)) {
        randomizeParticles(swarm);
    }
    if(IsKeyPressed(KEY_RIGHT_CONTROL)) {
        if (showForcesDBG)
            zeroParticlesForces();
        if (showMassesDBG)
            zeroParticlesMasses();
        if (showMinDistancesDBG)
            zeroParticlesMinDistances();
        if (showRadiiDBG)
            zeroParticlesRadii();
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
    if(IsKeyPressed(KEY_SPACE)) {
        forces = makeRandomMatrix();
        randomizeMasses();
        randomizeMinDistances();
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
        shift_particles_down(swarm);
    }
    if(IsKeyDown(KEY_DOWN) && IsKeyDown(KEY_LEFT_SHIFT)) {
        shift_particles_up(swarm);
    }
    if(IsKeyDown(KEY_LEFT) && IsKeyDown(KEY_LEFT_SHIFT)) {
        shift_particles_right(swarm);
    }
    if(IsKeyDown(KEY_RIGHT) && IsKeyDown(KEY_LEFT_SHIFT)) {
        shift_particles_left(swarm);
    }
    if(IsKeyPressed(KEY_UP)) {
        if (showForcesDBG)
            forces[cell_row_hover][cell_col_hover] += 0.1f;
        if (showMinDistancesDBG)
            minDistances[cell_row_hover][cell_col_hover] += 0.1f;
        if (showRadiiDBG)
            radii[cell_row_hover][cell_col_hover] += 0.1f;
        if (showMassesDBG)
            masses[cell_col_hover] += 0.1f;
    }
    if(IsKeyPressed(KEY_DOWN)) {
        if (showForcesDBG)
            forces[cell_row_hover][cell_col_hover] -= 0.1f;
        if (showMinDistancesDBG)
            minDistances[cell_row_hover][cell_col_hover] -= 0.1f;
        if (showRadiiDBG)
            radii[cell_row_hover][cell_col_hover] -= 0.1f;
        if (showMassesDBG)
            masses[cell_col_hover] -= 0.1f;
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

float force(float beta, float r, float a)
{
    if (r < beta) {
        return r / beta - 1.0f;
    } else if (beta < r && r < 1.0f) {
        return a * (1.0f - fabsf(2.0f * r - 1.0f - beta) / (1.0f - beta));
    } else {
        return 0.0f;
    }
}

#undef FORCE_EASE // experimental
void updateParticles(particle *this)
{
#ifdef FORCE_EASE
    float totalVelocity_len_max = 0.0f;
    float totalVelocity_len_min = FLT_MAX;
    float totalVelocity_len = 0.0f;
    float DT = 0.0f;
#endif /* FORCE_EASE */

    Vector2 totalForce = Vector2Zero();
    Vector2 totalVelocity = Vector2Zero();

    for (int i = 0; i < numParticles; i++) {
        totalForce = Vector2Zero();
        totalVelocity = Vector2Zero();

        for (int j = 0; j < numParticles; j++) {
            if (j == i) continue;
            Vector2 rv = Vector2Subtract(this[j].position, this[i].position);

            if(rv.x > 0.5f * WIDTH) {
                rv.x -= WIDTH;
            }
            if(rv.x < -0.5f * WIDTH) {
                rv.x += WIDTH;
            }
            if(rv.y > 0.5f * HEIGHT) {
                rv.y -= HEIGHT;
            }
            if(rv.y < -0.5f * HEIGHT) {
                rv.y += HEIGHT;
            }

            float r = sqrt(rv.x*rv.x + rv.y*rv.y);
            if (r <= maxRadius) {
                if (r > 0 && r < maxRadius) {
                    float f = force(minDistances[this[i].type][this[j].type], r / maxRadius, forces[this[i].type][this[j].type]);
                    //totalForce += rv / r * f;
                    Vector2 v = Vector2Scale(rv, 1.0f / r);
                    v = Vector2Scale(v, f);
                    totalForce = Vector2Add(totalForce, v);
                }
            }
        }

        totalForce = Vector2Scale(totalForce, maxRadius * forceFactor);

        this[i].velocity = Vector2Scale(this[i].velocity, frictionFactor);

        Vector2 v = Vector2Scale(totalForce, dt);
        v = Vector2Scale(v, 1/this[i].mass);
        this[i].velocity = Vector2Add(this[i].velocity, v);

        totalVelocity = Vector2Add(totalVelocity, this[i].velocity);

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
        Vector2 v = Vector2Scale(this[i].velocity, dt);
        this[i].position = Vector2Add(this[i].position, v);

        // wrap around - toroidal world
        this[i].position.x = fmodf(this[i].position.x + (float)WIDTH, (float)WIDTH);
        this[i].position.y = fmodf(this[i].position.y + (float)HEIGHT, (float)HEIGHT);
    }

    real_time += dt;
}

void drawParticles(particle *this)
{
    //ClearBackground(BLACK);

    for (int i = 0; i < numParticles; i++) {
        float screenX = this[i].position.x;
        float screenY = this[i].position.y;
//        Color col = ColorFromHSV(360.0f * ((float)particleType[i] / (float)numTypes), 1.0f, 1.0);
        Color col = ColorFromHSV(this[i].type*hueStep, 1.0f, 1.0);
        Vector2 world = {screenX, screenY};
//        DrawPixelV(WorldToScreen(world), col);
        float r = 2.0f * Vector2Length(scale);
        r = r > 1.0f ? r : 1.0f;
        world = WorldToScreen(world);
        DrawRectangle(world.x, world.y, r, r, col);
        //DrawCircleV(WorldToScreen(world), r, col);
    }
}

int main(int argc, char *argv[])
{
    srand48(clock());
    frameCounter = 0;

    // Initialization
    showHelp = false;
    showHUD = false;
    showForcesDBG = false;
    showMinDistancesDBG = false;
    showMassesDBG = false;
    showRadiiDBG = false;

    pause = false;
    display = true;

    numParticles = 1000;
    numTypes = 6;
    dt = 0.02f;
    frictionHalfLife = 0.040f;
    minRadius = 50;
    maxRadius = 150;
    minParticleMass = 0.001f;
    maxParticleMass = 2.0f;
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

    minDistances = NULL;
    arrsetlen(minDistances, numTypes);
    for (int i = 0; i < numTypes; ++i) {
        minDistances[i] = NULL;
        arrsetlen(minDistances[i], numTypes);
        for (int j = 0; j < numTypes; ++j) {
            minDistances[i][j] = map(random_range(0.0f, minRadius), 0.0f, minRadius, 0.0f, minRadius/maxRadius);
        }
    }
    radii = NULL;
    arrsetlen(radii, numTypes);
    for (int i = 0; i < numTypes; ++i) {
        radii[i] = NULL;
        arrsetlen(radii[i], numTypes);
        for (int j = 0; j < numTypes; ++j) {
            radii[i][j] = random_range(0.3f, 1.0f);
        }
    }
    arrsetlen(masses, numTypes);
    for (int i = 0; i < numTypes; ++i) {
        masses[i] = random_range(minParticleMass, maxParticleMass);
    }

    arrsetlen(swarm, numParticles);
    swarm_len = numParticles;

    for (int i = 0; i < swarm_len; ++i) {
        swarm[i].mass = masses[swarm[i].type];
    }

    randomizeParticles(swarm);

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
            updateParticles(swarm);

        if (showForcesDBG || showMinDistancesDBG || showRadiiDBG || showMassesDBG) {
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
                drawParticles(swarm);

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
