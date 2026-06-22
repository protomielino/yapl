#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "sim.h"
#include "qt.h"

extern const float WORLD_SIZE;
extern const int USE_PERIODIC;
extern sim *CURRENT_SIM_INSTANCE;

extern void draw_line_world(float x1, float y1, float x2, float y2);

#define QT_CAPACITY 8
#define QT_MAX_DEPTH 16

typedef struct {
    int index;
    float x, y;
} qt_point;

typedef struct qt_node_s {
    float x0, y0, x1, y1;
    struct qt_node_s *children[4];
    qt_point *points;
    int count;
    int capacity;
    int depth;
    int is_leaf;
} qt_node;

struct qt_s {
    qt_node *root;
};

static qt_node *qt_node_create(float x0, float y0, float x1, float y1, int depth)
{
    qt_node *n = (qt_node*)malloc(sizeof(qt_node));
    n->x0 = x0; n->y0 = y0;
    n->x1 = x1; n->y1 = y1;
    for (int i = 0; i < 4; i++)
        n->children[i] = NULL;
    n->points = NULL;
    n->count = 0;
    n->capacity = 0;
    n->depth = depth;
    n->is_leaf = 1;
    return n;
}

static void qt_node_free(qt_node *n)
{
    if (!n) return;
    for (int i = 0; i < 4; i++) {
        if (n->children[i]) {
            qt_node_free(n->children[i]);
            free(n->children[i]);
        }
    }
    free(n->points);
}

static int get_quadrant(qt_node *n, float x, float y)
{
    float mx = 0.5f * (n->x0 + n->x1);
    float my = 0.5f * (n->y0 + n->y1);
    int right = x >= mx;
    int top = y >= my;
    if (!right && !top) return 0;
    if (right && !top)  return 1;
    if (!right && top)  return 2;
    return 3;
}

static void subdivide(qt_node *n)
{
    float mx = 0.5f * (n->x0 + n->x1);
    float my = 0.5f * (n->y0 + n->y1);
    n->children[0] = qt_node_create(n->x0, n->y0, mx, my, n->depth + 1);
    n->children[1] = qt_node_create(mx, n->y0, n->x1, my, n->depth + 1);
    n->children[2] = qt_node_create(n->x0, my, mx, n->y1, n->depth + 1);
    n->children[3] = qt_node_create(mx, my, n->x1, n->y1, n->depth + 1);
    n->is_leaf = 0;
    for (int i = 0; i < n->count; i++) {
        int quad = get_quadrant(n, n->points[i].x, n->points[i].y);
        qt_node *child = n->children[quad];
        if (child->count >= child->capacity) {
            child->capacity = child->capacity ? child->capacity * 2 : 8;
            child->points = (qt_point*)realloc(child->points, child->capacity * sizeof(qt_point));
        }
        child->points[child->count++] = n->points[i];
    }
    free(n->points);
    n->points = NULL;
    n->count = 0;
    n->capacity = 0;
}

static void insert_rec(qt_node *n, int idx, float x, float y)
{
    if (!n->is_leaf) {
        int quad = get_quadrant(n, x, y);
        insert_rec(n->children[quad], idx, x, y);
        return;
    }
    if (n->count >= n->capacity) {
        n->capacity = n->capacity ? n->capacity * 2 : 8;
        n->points = (qt_point*)realloc(n->points, n->capacity * sizeof(qt_point));
    }
    n->points[n->count].index = idx;
    n->points[n->count].x = x;
    n->points[n->count].y = y;
    n->count++;
    if (n->count > QT_CAPACITY && n->depth < QT_MAX_DEPTH)
        subdivide(n);
}

static int box_circle_overlap(qt_node *n, float qx, float qy, float r)
{
    float closest_x = fmaxf(n->x0, fminf(qx, n->x1));
    float closest_y = fmaxf(n->y0, fminf(qy, n->y1));
    float dx = qx - closest_x;
    float dy = qy - closest_y;
    return (dx*dx + dy*dy) <= r*r;
}

static void query_rec(qt_node *n, float qx, float qy, float radius,
                       int *out, int *count, int max_out)
{
    if (!n || *count >= max_out) return;
    if (!box_circle_overlap(n, qx, qy, radius)) return;

    if (n->is_leaf) {
        float r2 = radius * radius;
        for (int i = 0; i < n->count; i++) {
            float dx = n->points[i].x - qx;
            float dy = n->points[i].y - qy;
            if (dx*dx + dy*dy <= r2)
                out[(*count)++] = n->points[i].index;
        }
        return;
    }
    for (int k = 0; k < 4; k++) {
        if (n->children[k])
            query_rec(n->children[k], qx, qy, radius, out, count, max_out);
    }
}

static void query_collect(qt_node *n, float qx, float qy, float radius,
                           int *out, int *count, int max_out, int *mark)
{
    if (!n || *count >= max_out) return;
    if (!box_circle_overlap(n, qx, qy, radius)) return;

    if (n->is_leaf) {
        float r2 = radius * radius;
        for (int i = 0; i < n->count; i++) {
            float dx = n->points[i].x - qx;
            float dy = n->points[i].y - qy;
            if (dx*dx + dy*dy <= r2) {
                int idx = n->points[i].index;
                if (!mark[idx]) {
                    out[(*count)++] = idx;
                    mark[idx] = 1;
                }
            }
        }
        return;
    }
    for (int k = 0; k < 4; k++) {
        if (n->children[k])
            query_collect(n->children[k], qx, qy, radius, out, count, max_out, mark);
    }
}

qt* qt_create(float world_size)
{
    (void)world_size;
    qt *q = (qt*)malloc(sizeof(qt));
    q->root = qt_node_create(0.0f, 0.0f, world_size, world_size, 0);
    return q;
}

void qt_free(qt *q)
{
    if (!q) return;
    if (q->root) {
        qt_node_free(q->root);
        free(q->root);
    }
    free(q);
}

void qt_insert(qt *q, int idx, float x, float y)
{
    if (!q || !q->root) return;
    insert_rec(q->root, idx, x, y);
}

int qt_query_radius(qt *q, float x, float y, float radius, int *out, int max_out)
{
    if (!q || !q->root || !out || max_out <= 0) return 0;
    int count = 0;

    if (USE_PERIODIC) {
        int *mark = CURRENT_SIM_INSTANCE && CURRENT_SIM_INSTANCE->neighbor_mark
                    ? CURRENT_SIM_INSTANCE->neighbor_mark
                    : NULL;
        int local_mark = 0;
        if (!mark) {
            mark = (int*)calloc(CURRENT_SIM_INSTANCE ? CURRENT_SIM_INSTANCE->n : 0, sizeof(int));
            local_mark = 1;
        } else {
            memset(mark, 0, sizeof(int) * CURRENT_SIM_INSTANCE->n);
        }
        for (int ox = -1; ox <= 1 && count < max_out; ox++)
            for (int oy = -1; oy <= 1 && count < max_out; oy++)
                query_collect(q->root, x + ox * WORLD_SIZE, y + oy * WORLD_SIZE,
                              radius, out, &count, max_out, mark);
        if (local_mark) free(mark);
    } else {
        query_rec(q->root, x, y, radius, out, &count, max_out);
    }
    return count;
}

static void draw_node(qt_node *n)
{
    if (!n) return;
    draw_line_world(n->x0, n->y0, n->x1, n->y0);
    draw_line_world(n->x1, n->y0, n->x1, n->y1);
    draw_line_world(n->x1, n->y1, n->x0, n->y1);
    draw_line_world(n->x0, n->y1, n->x0, n->y0);
    for (int k = 0; k < 4; k++) {
        if (n->children[k])
            draw_node(n->children[k]);
    }
}

void qt_draw_partitions(qt *q)
{
    if (!q) return;
    draw_node(q->root);
}
