#include <stdlib.h>
#include <math.h>

#include "grid.h"

extern void draw_line_world(float x1, float y1, float x2, float y2);
extern const float WORLD_SIZE;

static void int_array_push(int_array *a, int val)
{
    if (a->count >= a->capacity) {
        a->capacity = a->capacity ? a->capacity * 2 : 8;
        a->data = (int*)realloc(a->data, a->capacity * sizeof(int));
    }
    a->data[a->count++] = val;
}

grid* grid_create(float cell_size, float world_size)
{
    grid *g = (grid*)malloc(sizeof(grid));
    g->inv_cell_size = 1.0f / cell_size;
    g->ncells = (int)ceilf(world_size * g->inv_cell_size);
    if ((float)g->ncells * cell_size < world_size)
        g->ncells++;
    g->cells = (int_array*)calloc(g->ncells * g->ncells, sizeof(int_array));
    return g;
}

void grid_free(grid *g)
{
    if (!g) return;
    int total = g->ncells * g->ncells;
    for (int i = 0; i < total; i++)
        free(g->cells[i].data);
    free(g->cells);
    free(g);
}

void grid_clear(grid *g)
{
    int total = g->ncells * g->ncells;
    for (int i = 0; i < total; i++)
        g->cells[i].count = 0;
}

static inline int clip_cell(int c, int ncells)
{
    if (c < 0) return 0;
    if (c >= ncells) return ncells - 1;
    return c;
}

void grid_insert(grid *g, int idx, float x, float y)
{
    int cx = (int)(x * g->inv_cell_size);
    int cy = (int)(y * g->inv_cell_size);
    cx = clip_cell(cx, g->ncells);
    cy = clip_cell(cy, g->ncells);
    int linear = cx + cy * g->ncells;
    int_array_push(&g->cells[linear], idx);
}

int grid_query(grid *g, float x, float y, int *out, int max_out)
{
    if (!g || !out) return 0;

    int cx = (int)(x * g->inv_cell_size);
    int cy = (int)(y * g->inv_cell_size);
    int count = 0;

    for (int dy = -1; dy <= 1; dy++) {
        for (int dx = -1; dx <= 1; dx++) {
            int nx = cx + dx;
            int ny = cy + dy;
            if (nx < 0) nx += g->ncells;
            else if (nx >= g->ncells) nx -= g->ncells;
            if (ny < 0) ny += g->ncells;
            else if (ny >= g->ncells) ny -= g->ncells;
            int linear = nx + ny * g->ncells;
            int_array *cell = &g->cells[linear];
            for (int k = 0; k < cell->count && count < max_out; k++)
                out[count++] = cell->data[k];
        }
    }
    return count;
}

void grid_draw_partitions(grid *g)
{
    if (!g) return;
    float cell_w = WORLD_SIZE / (float)g->ncells;
    for (int i = 0; i <= g->ncells; i++) {
        float t = (float)i * cell_w;
        draw_line_world(t, 0.0f, t, WORLD_SIZE);
        draw_line_world(0.0f, t, WORLD_SIZE, t);
    }
}
