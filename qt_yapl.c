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
bool show_qt;
bool drawPixel = false;
bool pause = true;

// --------------------------------------------------------

// global simulation parameters (some set from original)
int n = 2000;
float dt = 0.02f;
float frictionHalfLife = 0.040f;
float frictionFactor;
float forceFactor = 10.0f;
float base_rMax = 0.1f;
float rMax; // may adapt
int m_colors = 6;
float beta = 0.3f;
float cutoff_factor = 0.005f; // forces smaller than cutoff_factor * rMax ignored
int seed_random = 12345;

// adaptive quadtree params (reasonable defaults; adjusted at runtime)
int base_capacity = 8;
int base_maxDepth = 16;

void showHUD()
{
    if(show_HUD) {
        sprintf(text, "frame time: %.3f", GetFrameTime());
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
        sprintf(text, "n particles: %d", n);
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "qt base capacity: %d [c,C]", base_capacity);
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "qt max tree depth: %d [d,D]", base_maxDepth);
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        ty += dy;
    }
}

void showHelp()
{
    if(show_help) {
        sprintf(text, "[R]\t reset pan & zoom");
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "[MMB] reset pan");
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        Color c = {};
        if (pause) c = RED;
        else c = WHITE;
        sprintf(text, "[p] un/pause");
        DrawText(text, 10, ty, textDim, c);
        ty += dy;
        if (show_qt) c = RED;
        else c = WHITE;
        sprintf(text, "[q] show tree");
        DrawText(text, 10, ty, textDim, c);
        ty += dy;
        c = WHITE;
        sprintf(text, "[w] base_capacity 8 <-> n");
        DrawText(text, 10, ty, textDim, c);
        ty += dy;
        if (drawPixel) c = RED;
        else c = WHITE;
        sprintf(text, "[e] pixel <-> circles");
        DrawText(text, 10, ty, textDim, c);
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
    if(IsKeyPressed(KEY_C) && !IsKeyDown(KEY_LEFT_SHIFT)) {
        base_capacity --;
        if (base_capacity < 1)
            base_capacity = 1;
    }
    if(IsKeyPressed(KEY_C) && IsKeyDown(KEY_LEFT_SHIFT)) {
        base_capacity ++;
    }
    if(IsKeyPressed(KEY_D) && !IsKeyDown(KEY_LEFT_SHIFT)) {
        base_maxDepth --;
        if (base_maxDepth < 0)
            base_maxDepth = 0;
    }
    if(IsKeyPressed(KEY_D) && IsKeyDown(KEY_LEFT_SHIFT)) {
        base_maxDepth ++;
    }
    if(IsKeyPressed(KEY_Q)) {
        show_qt = !show_qt;
    }
    if(IsKeyPressed(KEY_W)) {
        if (base_capacity == n)
            base_capacity = 8;
        else
            base_capacity = n;
    }
    if(IsKeyPressed(KEY_E)) {
        drawPixel = !drawPixel;
    }
    if(IsKeyPressed(KEY_P)) {
        pause = !pause;
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

// --------------------------------------------------------

typedef struct qtNode_s
{
    // square region [x0,x1) x [y0,y1)
    float x0, y0;
    float x1, y1;
    // center of mass and total mass
    float com_x, com_y;
    float mass;
    // children: 4 quads TL TR BL BR
    struct qtNode_s *children[4];
    // indices of particles directly stored (if leaf)
    int *indices; // dynamic array via stb_ds
    int capacity; // leaf capacity threshold (for adaptive)
    int depth;
    int is_leaf;
} qtNode;

// quadtree helpers
qtNode *qt_create(float x0, float y0, float x1, float y1, int capacity, int depth);
void qt_free(qtNode *q);
int qt_get_quadrant(qtNode *q, float x, float y);
void qt_subdivide(qtNode *q);
void qt_insert_index(qtNode *q, int idx, float x, float y, int capacity, int maxDepth, particle *particles);
void qt_insert(qtNode *q, int idx, particle *particles, int capacity, int maxDepth);
void qt_compute_com(qtNode *q, particle *particles);
void qt_accumulate_force(qtNode *q, particle *particles, int i, float rMax_local, float forceFactorLocal, float cutoff_abs, float *fx, float *fy, float *matrix, int m_colors, float theta);
void qt_draw(qtNode *q, int depth_limit);

// mock drawing functions
void draw_line(float x0, float y0, float x1, float y1);   // coordinate in [0,1)
void draw_rect(float x0, float y0, float x1, float y1);   // bordo rettangolo
// clamp wrap into [0,1)
float wrap01(float v);

// toroidal distance: compute dx,dy minimal with wrap on unit square
static inline void toroidal_delta(float x1, float y1, float x2, float y2, float *dx, float *dy)
{
    float rx = x2 - x1;
    float ry = y2 - y1;
    if (rx > 0.5f)
        rx -= 1.0f;
    else if (rx < -0.5f)
        rx += 1.0f;
    if (ry > 0.5f)
        ry -= 1.0f;
    else if (ry < -0.5f)
        ry += 1.0f;
    *dx = rx;
    *dy = ry;
}

// force function: r is normalized (r/rMax), a is matrix entry
static inline float force_func(float rnorm, float a)
{
    // consistent comparisons. use <= and >=
    if (rnorm <= beta) {
        return rnorm / beta - 1.0f;
    } else if (rnorm > beta && rnorm < 1.0f) {
        // triangular shaped scaled by a
        return a * (1.0f - fabsf(2.0f * rnorm - 1.0f - beta) / (1.0f - beta));
    } else {
        return 0.0f;
    }
}

// quadtree implementation
qtNode *qt_create(float x0, float y0, float x1, float y1, int capacity, int depth)
{
    qtNode *q = (qtNode*)malloc(sizeof(qtNode));
    q->x0 = x0; q->y0 = y0;
    q->x1 = x1; q->y1 = y1;
    q->com_x = 0.0f;
    q->com_y = 0.0f;
    q->mass = 0.0f;
    for (int i = 0; i < 4; i++) {
        q->children[i] = NULL;
    }
    q->indices = NULL;
    q->capacity = capacity;
    q->depth = depth;
    q->is_leaf = 1;
    return q;
}

void qt_free(qtNode *q)
{
    if (!q) // recursion exit condition
        return;
    for (int i = 0; i < 4; i++) {
        if (q->children[i]) {
            qt_free(q->children[i]);
            free(q->children[i]);
            q->children[i] = NULL;
        }
    }
    if (arrlen(q->indices) > 0)
        arrfree(q->indices);
}

int qt_get_quadrant(qtNode *q, float x, float y)
{
    float mx = 0.5f*(q->x0 + q->x1);
    float my = 0.5f*(q->y0 + q->y1);
    int right = x >= mx;
    int top = y >= my;
    // indexing: 0: BL,1: BR,2: TL,3: TR  (arbitrary)
    if (!right && !top)
        return 0; // BL
    if (right && !top)
        return 1;  // BR
    if (!right && top)
        return 2;  // TL
    return 3; // TR
}

void qt_subdivide(qtNode *q)
{
    float mx = 0.5f*(q->x0 + q->x1);
    float my = 0.5f*(q->y0 + q->y1);
    int child_cap = q->capacity;
    int child_depth = q->depth + 1;

    q->children[0] = qt_create(q->x0, q->y0, mx, my, child_cap, child_depth); // BL
    q->children[1] = qt_create(mx, q->y0, q->x1, my, child_cap, child_depth); // BR
    q->children[2] = qt_create(q->x0, my, mx, q->y1, child_cap, child_depth); // TL
    q->children[3] = qt_create(mx, my, q->x1, q->y1, child_cap, child_depth); // TR

    q->is_leaf = 0;
    // move existing indices into children
    for (int i = 0; i < arrlen(q->indices); i++) {
        // don't have particle array pointer here; insertion handled outside for simplicity
//        int idx = q->indices[i];
    }
    // keep indices empty in parent; caller will re-insert
    if (arrlen(q->indices) > 0)
        arrfree(q->indices);
    q->indices = NULL;
}

// insert index into quadtree leaf (caller ensures coordinates known)
void qt_insert_index(qtNode *q, int idx, float x, float y, int capacity, int maxDepth, particle *particles)
{
    // If this node isn't leaf, delegate to appropriate child
    if (!q->is_leaf) { // recursion exit condition
        int quad = qt_get_quadrant(q,x,y);
        qt_insert_index(q->children[quad], idx, x, y, capacity, maxDepth, particles);
        return;
    }
    // leaf
    arrpush(q->indices, idx);
    // if exceeds capacity and not at max depth, subdivide
    if (arrlen(q->indices) > q->capacity && q->depth < maxDepth) {
        // subdivide and reinsert
        qt_subdivide(q);
        // reinsert stored indices
//        int cnt = arrlen(q->indices); // but we freed indices in qt_subdivide; workaround: we captured earlier
        // instead: we will re-insert using local copy -- simpler approach: when subdividing, we must have saved indices
        // to avoid complex code, change approach: don't free indices inside qt_subdivide; instead move logic here.
    }
}

// simpler robust insertion implementation. recursive with reinsertion handled here
void qt_insert(qtNode *q, int idx, particle *particles, int capacity, int maxDepth)
{
    float x = particles[idx].position.x;
    float y = particles[idx].position.y;
    // if not leaf, delegate
    if (!q->is_leaf) { // recursion exit condition
        int quad = qt_get_quadrant(q, x, y);
        qt_insert(q->children[quad], idx, particles, capacity, maxDepth);
        return;
    }
    // leaf, store
    arrpush(q->indices, idx);
    // subdivide if needed
    if (arrlen(q->indices) > q->capacity && q->depth < maxDepth) {
        // save current indices
        int *saved = NULL;
        for (int i = 0; i < arrlen(q->indices); i++)
            arrpush(saved, q->indices[i]);
        if (arrlen(q->indices) > 0)
            arrfree(q->indices);
        // subdivide
        float mx = 0.5f*(q->x0 + q->x1);
        float my = 0.5f*(q->y0 + q->y1);
        int child_cap = q->capacity;
        int child_depth = q->depth + 1;
        q->children[0] = qt_create(q->x0, q->y0, mx, my, child_cap, child_depth);
        q->children[1] = qt_create(mx, q->y0, q->x1, my, child_cap, child_depth);
        q->children[2] = qt_create(q->x0, my, mx, q->y1, child_cap, child_depth);
        q->children[3] = qt_create(mx, my, q->x1, q->y1, child_cap, child_depth);
        q->is_leaf = 0;
        // reinsert saved
        for (int i = 0; i < arrlen(saved); i++) {
            qt_insert(q, saved[i], particles, capacity, maxDepth);
        }
        if (arrlen(saved) > 0)
            arrfree(saved);
    }
}

// build center of mass recursively
void qt_compute_com(qtNode *q, particle *particles)
{
    q->mass = 0.0f;
    q->com_x = 0.0f;
    q->com_y = 0.0f;
    if (q->is_leaf) {
        for (int i = 0; i < arrlen(q->indices); i++) {
            particle *p = &particles[q->indices[i]];
            q->mass += p->mass;
            q->com_x += p->position.x * p->mass;
            q->com_y += p->position.y * p->mass;
        }
    } else {
        for (int c = 0; c < 4; c++) {
            if (q->children[c]) {
                qt_compute_com(q->children[c], particles);
                q->mass += q->children[c]->mass;
                q->com_x += q->children[c]->com_x * q->children[c]->mass;
                q->com_y += q->children[c]->com_y * q->children[c]->mass;
            }
        }
    }
    if (q->mass > 0.0f) {
        q->com_x /= q->mass;
        q->com_y /= q->mass;
        // ensure com in [0,1)
        q->com_x = wrap01(q->com_x);
        q->com_y = wrap01(q->com_y);
    } else {
        q->com_x = 0.5f*(q->x0+q->x1);
        q->com_y = 0.5f*(q->y0+q->y1);
    }
}

// compute interaction between particle i and a node q, accumulating force into fx,fy
void qt_accumulate_force(qtNode *q, particle *particles, int i, float rMax_local, float forceFactorLocal, float cutoff_abs, float *fx, float *fy, float *matrix, int m_colors, float theta)
{
    if (!q || q->mass <= 0.0f)
        return;
    particle *pi = &particles[i];
    // compute toroidal delta between pi and node COM
    float dx, dy;
    toroidal_delta(pi->position.x, pi->position.y, q->com_x, q->com_y, &dx, &dy);
    float dist = hypotf(dx, dy);
    if (dist == 0.0f) {
        // avoid self-coupling; if node contains only the particle itself, skip
        if (q->is_leaf && arrlen(q->indices) == 1 && q->indices[0] == i)
            return;
        // otherwise, perturb slightly
        dist = 1e-6f;
    }
    // if node is sufficiently far relative to its size, approximate
    float width = q->x1 - q->x0;
    if (!q->is_leaf && (width / dist) < theta) {
        // approximate using COM as single particle with aggregated mass and average color?
        // for color-dependent interactions we need to choose how to handle mixing: we approximate using weighted average of colors index -> choose nearest color by weighted counts would be expensive.
        // simpler approach: compute interaction using aggregate mass and average color index weighted by mass to pick matrix row/col: use nearest integer color.
        // to keep compatibility, we'll compute average color weighted by mass (floating) and clamp to [0,m-1].
        // but we don't have aggregated per-color counts; accept approximation by using color of the nearest particle in node if available.
        int color_j = 0;
        // find a representative color: search a leaf descendant for first particle
        qtNode *scan = q;
        while (!scan->is_leaf) {
            int found = 0;
            for (int c = 0; c < 4; c++) {
                if (scan->children[c] && scan->children[c]->mass > 0.0f) {
                    scan = scan->children[c];
                    found = 1;
                    break;
                }
            }
            if (!found)
                break;
        }
        if (scan->is_leaf && arrlen(scan->indices)>0) {
            color_j = particles[scan->indices[0]].type;
        }
        // normalized distance
        float rnorm = dist / rMax_local;
        if (rnorm >= 1.0f)
            return;
        float a = matrix[pi->type * m_colors + color_j];
        float fval = force_func(rnorm, a);
        // apply cutoff on absolute magnitude
        if (fabsf(fval) < cutoff_abs)
            return;
        // direction from i to j: dx,dy (already from i to com)
        float invd = 1.0f / dist;
        float fx_inc = (dx * invd) * fval * rMax_local * forceFactorLocal * (q->mass); // scale by node mass
        float fy_inc = (dy * invd) * fval * rMax_local * forceFactorLocal * (q->mass);
        // divide by particle mass (F = m_node * f -> acceleration depends on particle mass, we will later divide by pi->mass)
        *fx += fx_inc;
        *fy += fy_inc;
        return;
    }
    // if leaf, iterate over contained particles individually
    if (q->is_leaf) {
        for (int idx_i=0; idx_i < arrlen(q->indices); idx_i++) {
            int j = q->indices[idx_i];
            if (j == i)
                continue;
            particle *pj = &particles[j];
            float rx, ry;
            toroidal_delta(pi->position.x, pi->position.y, pj->position.x, pj->position.y, &rx, &ry);
            float r = hypotf(rx, ry);
            if (r <= 0.0f || r >= rMax_local)
                continue;
            float rnorm = r / rMax_local;
            float a = matrix[pi->type * m_colors + pj->type];
            float fval = force_func(rnorm, a);
            if (fabsf(fval) < cutoff_abs)
                continue;
            float invr = 1.0f / r;
            float fx_inc = (rx * invr) * fval * rMax_local * forceFactorLocal * pj->mass;
            float fy_inc = (ry * invr) * fval * rMax_local * forceFactorLocal * pj->mass;
            *fx += fx_inc;
            *fy += fy_inc;
        }
        return;
    }
    // otherwise recurse into children
    for (int c = 0; c < 4; c++) {
        if (q->children[c] && q->children[c]->mass > 0.0f) {
            qt_accumulate_force(q->children[c], particles, i, rMax_local, forceFactorLocal, cutoff_abs, fx, fy, matrix, m_colors, theta);
        }
    }
}

// draw tree; depth_limit = -1 means no limit
void qt_draw(qtNode *q, int depth_limit)
{
    if (!q) // recursion exit condition
        return;
    // disegna il rettangolo del nodo corrente
    draw_rect(q->x0, q->y0, q->x1, q->y1);

    if (!q->is_leaf && (depth_limit < 0 || q->depth < depth_limit)) {
        for (int k = 0; k < 4; ++k) {
            if (q->children[k]) {
                qt_draw(q->children[k], depth_limit);
            }
        }
    }
}

unsigned int rnd_state;

static void init_globals()
{
    rnd_state = seed_random;
    frictionFactor = powf(0.5f, dt / frictionHalfLife);
    rMax = base_rMax;
}


void draw_line(float x0, float y0, float x1, float y1)
{
    float screenX0 = x0 * WIDTH;
    float screenY0 = y0 * HEIGHT;
    float screenX1 = x1 * WIDTH;
    float screenY1 = y1 * HEIGHT;
    Color c = WHITE; c.a = 64;
    Vector2 v0 = {screenX0, screenY0};
    Vector2 v1 = {screenX1, screenY1};
    v0 = WorldToScreen(v0);
    v1 = WorldToScreen(v1);

    DrawLine(v0.x, v0.y, v1.x, v1.y, c);
}
void draw_rect(float x0, float y0, float x1, float y1)
{
    draw_line(x0, y0, x1, y0);
    draw_line(x1, y0, x1, y1);
    draw_line(x1, y1, x0, y1);
    draw_line(x0, y1, x0, y0);
}

// create random interaction matrix [-1,1]
static float* makeRandomMatrix(int m)
{
    float *mat = malloc(sizeof(float)*m*m);
    for (int i=0;i<m*m;i++)
        mat[i] = frand()*2.0f - 1.0f;
    return mat;
}

int main(int argc, char *argv[])
{
//    srand48(clock());
//    srand(clock());
    seed_random = clock();
    frameCounter = 0;

    // Initialization

    init_globals();

    // allow override n via arg
    if (argc > 1)
        n = atoi(argv[1]);

    // allocate particles
    particle *particles = malloc(sizeof(particle) * n);
    if (!particles) { fprintf(stderr,"alloc fail\n"); return 1; }

    // randomize initial particles: positions uniform, small random vel, random color
    for (int i = 0; i < n; ++i) {
        particles[i].position.x = frand();
        particles[i].position.y = frand();
        particles[i].velocity.x = (frand() - 0.5f) * 0.01f;
        particles[i].velocity.y = (frand() - 0.5f) * 0.01f;
        particles[i].type = (int)(frand() * m_colors);
        particles[i].mass = 0.5f + frand() * 1.5f; // [0.5, 2.0]
    }

    // create random interaction matrix m_colors x m_colors in [-1,1]
    float *matrix = makeRandomMatrix(m_colors);

    // Initialise offset so 0,0 top left of the screen
    initPanZoom(Vector2Zero(), Vector2One());
    initHUD();

    //---------------------------------------------------------------------------------------

    InitWindow(WIDTH, HEIGHT, "Particle Life");
    //SetWindowState(FLAG_FULLSCREEN_MODE);
//    SetTargetFPS(30);
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

        // create root quadtree covering unit square
        qtNode *root = qt_create(0.0f, 0.0f, 1.0f, 1.0f, base_capacity, 0);

        // insert all particles into quadtree
        for (int i = 0; i < n; ++i) {
            qt_insert(root, i, particles, base_capacity, base_maxDepth);
        }

        // compute centers of mass
        qt_compute_com(root, particles);

        // prepare per-particle force accumulation and integrate one step
        float cutoff_abs = cutoff_factor * rMax;
        float theta = 0.5f; // Barnes-Hut opening angle, example value

        if (!pause) {
            for (int i = 0; i < n; ++i) {
                particle *pi = &particles[i];
                float fx = 0.0f, fy = 0.0f;
                qt_accumulate_force(root, particles, i, rMax, forceFactor, cutoff_abs, &fx, &fy, matrix, m_colors, theta);

                // apply friction
                pi->velocity.x = pi->velocity.x * frictionFactor + dt * fx / pi->mass;
                pi->velocity.y = pi->velocity.y * frictionFactor + dt * fy / pi->mass;

                // integrate positions with toroidal wrap
                pi->position.x = wrap01(pi->position.x + dt * pi->velocity.x);
                pi->position.y = wrap01(pi->position.y + dt * pi->velocity.y);
            }
        }

        // Draw
        //----------------------------------------------------------------------------------
        BeginDrawing(); {
            ClearBackground(BLACK);

            for (int i = 0; i < n; i++) {
                float screenX = particles[i].position.x * WIDTH;
                float screenY = particles[i].position.y * HEIGHT;
                Color c = ColorFromHSV(360.0f * ((float)particles[i].type / (float)m_colors), 1.0f, 1.0f);
                Vector2 v = {screenX, screenY};
                v = WorldToScreen(v);
                if (drawPixel) {
                    DrawPixel(v.x, v.y, c);
                } else {
                    float radius = fmaxf(0.5, particles[i].mass * 1.5); // tweak
                    radius *= scale.x;
                    radius = radius >= 1.0f ? radius*0.75f : 0.75f;
                    DrawCircle(v.x, v.y, radius, c);
                }
            }

            if (show_qt)
                qt_draw(root, -1); // -1 iterates up to the leaf.

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
        // cleanup
        //----------------------------------------------------------------------------------
        qt_free(root);
        free(root);
    }

    // De-Initialization
    free(particles);
    free(matrix);
    //--------------------------------------------------------------------------------------
    CloseWindow();        // Close window and OpenGL context
    //--------------------------------------------------------------------------------------

    return 0;
}
