#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "stb_ds.h"

#include "maths.h"
#include "sim.h"
#if USE_GRID
#include "grid.h"
#else
#include "kd.h"
#endif

extern bool show_tree;
extern const int USE_PERIODIC;
extern const float WORLD_SIZE;
extern const int WIDTH;
extern const int HEIGHT;
extern const int DEFAULT_M;

#if !USE_GRID
extern int KD_CURRENT_POINT_COUNT;
extern sim *CURRENT_SIM_INSTANCE;
#endif

extern void draw_particle_world(float x, float y, float size, int color, float mass);
#if !USE_GRID
extern void draw_line_world(float x1, float y1, float x2, float y2);
#endif

static inline float force_law(float r, float a, float beta)
{
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
    if (v >= 0.0f && v < WORLD_SIZE)
        return v;
    float r = fmodf(v, WORLD_SIZE);
    if (r < 0.0f)
        r += WORLD_SIZE;
    return r;
}

static inline void min_image_disp(float ax, float ay, float bx, float by, float *out_dx, float *out_dy)
{
    float dx = bx - ax;
    float dy = by - ay;
    if (!USE_PERIODIC) {
        *out_dx = dx;
        *out_dy = dy;
        return;
    }
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
    s->masses = NULL;
    s->matrix = (float*)malloc(sizeof(float) * m * m);
    s->min_dist_matrix = (float*)malloc(sizeof(float) * m * m);
    s->radii_matrix = (float*)malloc(sizeof(float) * m * m);

    arrsetlen(s->colors, n);
    arrsetlen(s->posx, n);
    arrsetlen(s->posy, n);
    arrsetlen(s->velx, n);
    arrsetlen(s->vely, n);
    arrsetlen(s->masses, n);

    if (USE_PERIODIC) {
        s->neighbor_mark = (int*)calloc(n, sizeof(int));
    } else {
        s->neighbor_mark = NULL;
    }

#if USE_GRID
    s->grid = grid_create(rmax, WORLD_SIZE);
#endif

    sim_randomize_colors(s);
    sim_randomize_positions(s);
    sim_randomize_masses(s, 0.5f, 2.0f);
    sim_randomize_matrix(s);
    sim_randomize_min_dist_matrix(s);
    sim_randomize_radii_matrix(s);
    return s;
}

void sim_free(sim *s)
{
    if (!s)
        return;
    arrfree(s->colors);
    arrfree(s->posx);
    arrfree(s->posy);
    arrfree(s->velx);
    arrfree(s->vely);
    arrfree(s->masses);
    free(s->matrix);
    free(s->min_dist_matrix);
    free(s->radii_matrix);
    if (s->neighbor_mark)
        free(s->neighbor_mark);
#if USE_GRID
    if (s->grid)
        grid_free(s->grid);
#endif
    free(s);
}

static inline void compute_forces_for_particle(sim *s, int i, const int *neighbors, int n_neighbors,
                                                float *out_fx, float *out_fy)
{
    float totalForceX = 0.0f;
    float totalForceY = 0.0f;
    float rmax = s->rmax;
    int mi = s->m;

    for (int k = 0; k < n_neighbors; ++k) {
        int j = neighbors[k];
        if (j == i)
            continue;
        float rx, ry;
        min_image_disp(s->posx[i], s->posy[i], s->posx[j], s->posy[j], &rx, &ry);
        float r = hypotf(rx, ry);
        if (r > 0.0f && r < rmax) {
            float invr = 1.0f / r;
            int ci = s->colors[i];
            int cj = s->colors[j];
            float a = s->matrix[ci * mi + cj];
            float beta = s->min_dist_matrix[ci * mi + cj];
            float f = force_law(r / rmax, a, beta);
            totalForceX += rx * invr * f;
            totalForceY += ry * invr * f;
        }
    }

    *out_fx = totalForceX * rmax * s->forceFactor;
    *out_fy = totalForceY * rmax * s->forceFactor;
}

void sim_update(sim *s)
{
    if (!s || s->n <= 0)
        return;

#if USE_GRID
    // ---- grid path ----
    grid_clear(s->grid);
    for (int i = 0; i < s->n; ++i) {
        grid_insert(s->grid, i, s->posx[i], s->posy[i]);
    }

    if (show_tree)
        grid_draw_partitions(s->grid);

    int *neighbors = NULL;
    for (int i = 0; i < s->n; ++i) {
        neighbors = NULL;
        grid_query(s->grid, s->posx[i], s->posy[i], &neighbors);

        float totalForceX, totalForceY;
        compute_forces_for_particle(s, i, neighbors, arrlen(neighbors), &totalForceX, &totalForceY);

        s->velx[i] *= s->frictionFactor;
        s->vely[i] *= s->frictionFactor;

        float ax = totalForceX / s->masses[i];
        float ay = totalForceY / s->masses[i];
        s->velx[i] += ax * s->dt;
        s->vely[i] += ay * s->dt;

        arrfree(neighbors);
    }

#else
    // ---- k-d tree path ----
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

    if (show_tree)
        kd_draw_tree_partitions_full(tree, 0);

    int *neighbors = NULL;
    for (int i = 0; i < s->n; ++i) {
        neighbors = NULL;
        kd_query_radius(tree, s->posx[i], s->posy[i], s->rmax, &neighbors);

        float totalForceX, totalForceY;
        compute_forces_for_particle(s, i, neighbors, arrlen(neighbors), &totalForceX, &totalForceY);

        s->velx[i] *= s->frictionFactor;
        s->vely[i] *= s->frictionFactor;

        float ax = totalForceX / s->masses[i];
        float ay = totalForceY / s->masses[i];
        s->velx[i] += ax * s->dt;
        s->vely[i] += ay * s->dt;

        arrfree(neighbors);
    }

    if (tree)
        kd_free(tree);
    arrfree(points);

    CURRENT_SIM_INSTANCE = NULL;
#endif

    for (int i = 0; i < s->n; ++i) {
        s->posx[i] += s->velx[i] * s->dt;
        s->posy[i] += s->vely[i] * s->dt;
        if (USE_PERIODIC) {
            s->posx[i] = wrap_coord(s->posx[i]);
            s->posy[i] = wrap_coord(s->posy[i]);
        } else {
            if (s->posx[i] < 0.0f) s->posx[i] = 0.0f;
            if (s->posx[i] >= WORLD_SIZE) s->posx[i] = WORLD_SIZE - 1e-6f;
            if (s->posy[i] < 0.0f) s->posy[i] = 0.0f;
            if (s->posy[i] >= WORLD_SIZE) s->posy[i] = WORLD_SIZE - 1e-6f;
        }
    }
}

void sim_draw_frame(sim *s)
{
    if (!s) return;
    for (int i = 0; i < s->n; ++i) {
        draw_particle_world(s->posx[i], s->posy[i], 1.5f, s->colors[i], s->masses[i]);
    }
}

void sim_get_positions(sim *s, float* out_x, float* out_y, int* out_colors)
{
    if (!s) return;
    for (int i = 0; i < s->n; ++i) {
        if (out_x) out_x[i] = s->posx[i];
        if (out_y) out_y[i] = s->posy[i];
        if (out_colors) out_colors[i] = s->colors[i];
    }
}

void sim_randomize_all(sim *s, float min_mass, float max_mass)
{
    sim_randomize_matrix(s);
    sim_randomize_min_dist_matrix(s);
    sim_randomize_radii_matrix(s);
    sim_randomize_masses(s, min_mass, max_mass);
    sim_randomize_positions(s);
    sim_randomize_colors(s);
}

void sim_randomize_matrix(sim *s)
{
    for (int i = 0; i < s->m * s->m; ++i)
        s->matrix[i] = randf() * 2.0f - 1.0f;
}

void sim_randomize_min_dist_matrix(sim *s)
{
    for (int i = 0; i < s->m * s->m; ++i)
        s->min_dist_matrix[i] = randf() * 0.35f;
}

void sim_randomize_radii_matrix(sim *s)
{
    for (int i = 0; i < s->m * s->m; ++i)
        s->radii_matrix[i] = 0.3f + randf() * 0.7f;
}

void sim_randomize_masses(sim *s, float min_mass, float max_mass)
{
    for (int i = 0; i < s->n; i++) {
        s->masses[i] = min_mass + randf() * (max_mass - min_mass);
    }
}

void sim_randomize_positions(sim *s)
{
    for (int i = 0; i < s->n; i++) {
        s->posx[i] = randf() * WORLD_SIZE;
        s->posy[i] = randf() * WORLD_SIZE;
        s->velx[i] = 0.0f;
        s->vely[i] = 0.0f;
    }
}

void sim_randomize_colors(sim *s)
{
    for (int i = 0; i < s->n; i++) {
        s->colors[i] = (int)(randf() * s->m);
    }
}
