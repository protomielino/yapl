#include <stdlib.h>
#include <math.h>

#include "stb_ds.h"

#include "maths.h"
#include "kd.h"
#include "sim.h"

extern bool show_tree;
extern const int USE_PERIODIC;
extern const float WORLD_SIZE;
extern int KD_CURRENT_POINT_COUNT;
extern sim *CURRENT_SIM_INSTANCE;
extern const int WIDTH;
extern const int HEIGHT;
extern const int DEFAULT_M;

// -------------------- drawing API --------------------
extern void draw_particle_world(float x, float y, float size, int color);

static inline float force_law(float r, float a)
{
    const float beta = 0.3f;
    if (r < beta) {
        return r / beta - 1.0f;
    } else if (beta < r && r < 1.0f) {
        return a * (1.0f - fabsf(2.0f * r - 1.0f - beta) / (1.0f - beta));
    } else {
        return 0.0f;
    }
}

static inline float wrap_coord(float v)
{
    if (!USE_PERIODIC)
        return v;
    // modulo in [0,WORLD_SIZE)
    if (v >= 0.0f && v < WORLD_SIZE)
        return v;
    // fmod can be negative, so adjust
    float r = fmodf(v, WORLD_SIZE);
    if (r < 0.0f)
        r += WORLD_SIZE;
    return r;
}

// return displacement from a to b using minimum-image convention on domain [0, WORLD_SIZE)
static inline void min_image_disp(float ax, float ay, float bx, float by, float *out_dx, float *out_dy)
{
    float dx = bx - ax;
    float dy = by - ay;
    if (!USE_PERIODIC) {
        *out_dx = dx;
        *out_dy = dy;
        return;
    }
    // wrap dx into [-WORLD_SIZE/2, WORLD_SIZE/2)
    float half = 0.5f * WORLD_SIZE;
    if (dx > half)
        dx -= WORLD_SIZE;
    else if (dx < -half)
        dx += WORLD_SIZE;
    if (dy > half)
        dy -= WORLD_SIZE;
    else if (dy < -half)
        dy += WORLD_SIZE;
    *out_dx = dx; *out_dy = dy;
}

sim* sim_create(int n, int m, float dt, float friction_half_life, float rmax, float forceFactor)
{
    sim* s = (sim*)malloc(sizeof(sim));
    s->n = n;
    s->m = m;
    s->dt = dt;
    s->friction_half_life = friction_half_life;
    s->rmax = rmax;
    s->forceFactor = forceFactor;
    s->frictionFactor = powf(0.5f, dt / friction_half_life);
    s->colors = NULL;
    s->posx = NULL;
    s->posy = NULL;
    s->velx = NULL;
    s->vely = NULL;
    s->matrix = (float*)malloc(sizeof(float) * m * m);

    arrsetlen(s->colors, n);
    arrsetlen(s->posx, n);
    arrsetlen(s->posy, n);
    arrsetlen(s->velx, n);
    arrsetlen(s->vely, n);

    if (USE_PERIODIC) {
        s->neighbor_mark = (int*)calloc(n, sizeof(int));
    } else {
        s->neighbor_mark = NULL;
    }

    for (int i = 0; i < n; ++i) {
        s->colors[i] = (int)(randf() * m);
        s->posx[i] = randf() * WORLD_SIZE;
        s->posy[i] = randf() * WORLD_SIZE;
        s->velx[i] = 0.0f;
        s->vely[i] = 0.0f;
    }
    for (int i = 0; i < m*m; ++i)
        s->matrix[i] = randf() * 2.0f - 1.0f;
    return s;
}

void sim_free(sim* s)
{
    if (!s)
        return;
    arrfree(s->colors);
    arrfree(s->posx);
    arrfree(s->posy);
    arrfree(s->velx);
    arrfree(s->vely);
    if (s->matrix)
        free(s->matrix);
    if (s->neighbor_mark)
        free(s->neighbor_mark);
    free(s);
}

void sim_update(sim* s)
{
    if (!s || s->n <= 0)
        return;

    // prepare points
    kd_point *points = NULL;
    arrsetlen(points, s->n);
    for (int i = 0; i < s->n; ++i) {
        points[i].x = s->posx[i];
        points[i].y = s->posy[i];
        points[i].index = i;
    }

    CURRENT_SIM_INSTANCE = s;
    KD_CURRENT_POINT_COUNT = s->n;

    kd_node *tree = kd_build(points, s->n);

    // debug dump (call once or periodically)
    if (show_tree)
        kd_draw_tree_partitions_full(tree, 0);

    int *neighbors = NULL;

    for (int i = 0; i < s->n; ++i) {
        float totalForceX = 0.0f;
        float totalForceY = 0.0f;

        neighbors = NULL; // ensure clean start
        kd_query_radius(tree, s->posx[i], s->posy[i], s->rmax, &neighbors);

        for (int k = 0; k < arrlen(neighbors); ++k) {
            int j = neighbors[k];
            if (j == i)
                continue;
            float rx, ry;
            min_image_disp(s->posx[i], s->posy[i], s->posx[j], s->posy[j], &rx, &ry);
            float r = hypotf(rx, ry);
            if (r > 0.0f && r < s->rmax) {
                float invr = 1.0f / r;
                float a = s->matrix[s->colors[i] * s->m + s->colors[j]];
                float f = force_law(r / s->rmax, a);
                totalForceX += rx * invr * f;
                totalForceY += ry * invr * f;
            }
        }

        totalForceX *= s->rmax * s->forceFactor;
        totalForceY *= s->rmax * s->forceFactor;

        s->velx[i] *= s->frictionFactor;
        s->vely[i] *= s->frictionFactor;

        s->velx[i] += totalForceX * s->dt;
        s->vely[i] += totalForceY * s->dt;

        arrfree(neighbors);
    }

    if (tree)
        kd_free(tree);
    arrfree(points);

    for (int i = 0; i < s->n; ++i) {
        s->posx[i] += s->velx[i] * s->dt;
        s->posy[i] += s->vely[i] * s->dt;
        // apply wrap/clamp
        if (USE_PERIODIC) {
            s->posx[i] = wrap_coord(s->posx[i]);
            s->posy[i] = wrap_coord(s->posy[i]);
        } else {
            // optional clamp to [0, WORLD_SIZE)
            if (s->posx[i] < 0.0f)
                s->posx[i] = 0.0f;
            if (s->posx[i] >= WORLD_SIZE)
                s->posx[i] = WORLD_SIZE - 1e-6f;
            if (s->posy[i] < 0.0f)
                s->posy[i] = 0.0f;
            if (s->posy[i] >= WORLD_SIZE)
                s->posy[i] = WORLD_SIZE - 1e-6f;
        }
    }

    CURRENT_SIM_INSTANCE = NULL;
}

void sim_draw_frame(sim* s)
{
    if (!s)
        return;
    for (int i = 0; i < s->n; ++i) {
        draw_particle_world(s->posx[i], s->posy[i], 1.5f, s->colors[i]);
    }
}

void sim_get_positions(sim* s, float* out_x, float* out_y, int* out_colors)
{
    if (!s)
        return;
    for (int i = 0; i < s->n; ++i) {
        if (out_x) out_x[i] = s->posx[i];
        if (out_y) out_y[i] = s->posy[i];
        if (out_colors) out_colors[i] = s->colors[i];
    }
}
