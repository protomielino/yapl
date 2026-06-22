#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "stb_ds.h"

#include "maths.h"
#include "sim.h"
#include "grid.h"
#include "kd.h"

extern bool show_tree;
extern const int USE_PERIODIC;
extern const float WORLD_SIZE;
extern const int WIDTH;
extern const int HEIGHT;

int KD_CURRENT_POINT_COUNT = 0;
sim *CURRENT_SIM_INSTANCE = NULL;

extern void draw_particle_world(float x, float y, float size, int color, float mass);

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

sim* sim_create(int n, int m, backend_type backend, float dt, float friction_half_life, float rmax, float forceFactor)
{
    sim* s = (sim*)malloc(sizeof(sim));
    s->n = n;
    s->m = m;
    s->backend = backend;
    s->dt = dt;
    s->friction_half_life = friction_half_life;
    s->rmax = rmax;
    s->forceFactor = forceFactor;
    s->frictionFactor = powf(0.5f, dt / friction_half_life);
    s->colors = NULL;
    s->positions = NULL;
    s->velocities = NULL;
    s->masses = NULL;
    s->matrix = (float*)malloc(sizeof(float) * m * m);
    s->min_dist_matrix = (float*)malloc(sizeof(float) * m * m);
    s->radii_matrix = (float*)malloc(sizeof(float) * m * m);

    arrsetlen(s->colors, n);
    arrsetlen(s->positions, n);
    arrsetlen(s->velocities, n);
    arrsetlen(s->masses, n);

    if (USE_PERIODIC) {
        s->neighbor_mark = (int*)calloc(n, sizeof(int));
    } else {
        s->neighbor_mark = NULL;
    }

    s->grid = NULL;
    s->kd_tree = NULL;
    s->kd_rebuild_counter = 0;
    if (backend == BACKEND_GRID)
        s->grid = grid_create(rmax, WORLD_SIZE);

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
    arrfree(s->positions);
    arrfree(s->velocities);
    arrfree(s->masses);
    free(s->matrix);
    free(s->min_dist_matrix);
    free(s->radii_matrix);
    if (s->neighbor_mark)
        free(s->neighbor_mark);
    if (s->grid)
        grid_free(s->grid);
    if (s->kd_tree)
        kd_free(s->kd_tree);
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
        min_image_disp(s->positions[i].x, s->positions[i].y,
                       s->positions[j].x, s->positions[j].y, &rx, &ry);
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

    if (s->backend == BACKEND_GRID) {
        // ---- grid path ----
        grid_clear(s->grid);
        for (int i = 0; i < s->n; ++i)
            grid_insert(s->grid, i, s->positions[i].x, s->positions[i].y);

        if (show_tree)
            grid_draw_partitions(s->grid);

        int *neighbors = NULL;
        for (int i = 0; i < s->n; ++i) {
            neighbors = NULL;
            grid_query(s->grid, s->positions[i].x, s->positions[i].y, &neighbors);

            float totalForceX, totalForceY;
            compute_forces_for_particle(s, i, neighbors, arrlen(neighbors), &totalForceX, &totalForceY);

            s->velocities[i].x *= s->frictionFactor;
            s->velocities[i].y *= s->frictionFactor;

            float ax = totalForceX / s->masses[i];
            float ay = totalForceY / s->masses[i];
            s->velocities[i].x += ax * s->dt;
            s->velocities[i].y += ay * s->dt;

            arrfree(neighbors);
        }
    } else {
        // ---- k-d tree path ----
        const int KD_REBUILD_INTERVAL = 10;

        CURRENT_SIM_INSTANCE = s;
        KD_CURRENT_POINT_COUNT = s->n;

        if (s->kd_rebuild_counter <= 0) {
            if (s->kd_tree) {
                kd_free(s->kd_tree);
                s->kd_tree = NULL;
            }
            kd_point *points = NULL;
            arrsetlen(points, s->n);
            for (int i = 0; i < s->n; ++i) {
                points[i].x = s->positions[i].x;
                points[i].y = s->positions[i].y;
                points[i].index = i;
            }
            s->kd_tree = kd_build(points, s->n);
            s->kd_rebuild_counter = KD_REBUILD_INTERVAL;
            arrfree(points);
        } else {
            s->kd_rebuild_counter--;
        }

        if (show_tree && s->kd_tree)
            kd_draw_tree_partitions_full(s->kd_tree, 0);

        int *neighbors = NULL;
        for (int i = 0; i < s->n; ++i) {
            neighbors = NULL;
            kd_query_radius(s->kd_tree, s->positions[i].x, s->positions[i].y, s->rmax, &neighbors);

            float totalForceX, totalForceY;
            compute_forces_for_particle(s, i, neighbors, arrlen(neighbors), &totalForceX, &totalForceY);

            s->velocities[i].x *= s->frictionFactor;
            s->velocities[i].y *= s->frictionFactor;

            float ax = totalForceX / s->masses[i];
            float ay = totalForceY / s->masses[i];
            s->velocities[i].x += ax * s->dt;
            s->velocities[i].y += ay * s->dt;

            arrfree(neighbors);
        }

        CURRENT_SIM_INSTANCE = NULL;
    }

    for (int i = 0; i < s->n; ++i) {
        s->positions[i].x += s->velocities[i].x * s->dt;
        s->positions[i].y += s->velocities[i].y * s->dt;
        if (USE_PERIODIC) {
            s->positions[i].x = wrap_coord(s->positions[i].x);
            s->positions[i].y = wrap_coord(s->positions[i].y);
        } else {
            if (s->positions[i].x < 0.0f) s->positions[i].x = 0.0f;
            if (s->positions[i].x >= WORLD_SIZE) s->positions[i].x = WORLD_SIZE - 1e-6f;
            if (s->positions[i].y < 0.0f) s->positions[i].y = 0.0f;
            if (s->positions[i].y >= WORLD_SIZE) s->positions[i].y = WORLD_SIZE - 1e-6f;
        }
    }
}

void sim_draw_frame(sim *s)
{
    if (!s) return;
    for (int i = 0; i < s->n; ++i) {
        draw_particle_world(s->positions[i].x, s->positions[i].y, 1.5f, s->colors[i], s->masses[i]);
    }
}

void sim_get_positions(sim *s, float* out_x, float* out_y, int* out_colors)
{
    if (!s) return;
    for (int i = 0; i < s->n; ++i) {
        if (out_x) out_x[i] = s->positions[i].x;
        if (out_y) out_y[i] = s->positions[i].y;
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
        s->matrix[i] = rand01() * 2.0f - 1.0f;
}

void sim_randomize_min_dist_matrix(sim *s)
{
    for (int i = 0; i < s->m * s->m; ++i)
        s->min_dist_matrix[i] = rand01() * 0.35f;
}

void sim_randomize_radii_matrix(sim *s)
{
    for (int i = 0; i < s->m * s->m; ++i)
        s->radii_matrix[i] = 0.3f + rand01() * 0.7f;
}

void sim_randomize_masses(sim *s, float min_mass, float max_mass)
{
    for (int i = 0; i < s->n; i++) {
        s->masses[i] = min_mass + rand01() * (max_mass - min_mass);
    }
}

void sim_randomize_positions(sim *s)
{
    for (int i = 0; i < s->n; i++) {
        s->positions[i].x = rand01() * WORLD_SIZE;
        s->positions[i].y = rand01() * WORLD_SIZE;
        s->velocities[i].x = 0.0f;
        s->velocities[i].y = 0.0f;
    }
    s->kd_rebuild_counter = 0;
}

void sim_randomize_colors(sim *s)
{
    for (int i = 0; i < s->n; i++) {
        s->colors[i] = (int)(rand01() * s->m);
    }
}
