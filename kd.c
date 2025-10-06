#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#include "stb_ds.h"

#include "sim.h"
#include "kd.h"

extern const float WORLD_SIZE;
extern const int USE_PERIODIC;
extern sim *CURRENT_SIM_INSTANCE;
extern int KD_CURRENT_POINT_COUNT;

extern void draw_line_world(float x1, float y1, float x2, float y2);

// squared distance point-rect for bbox pruning
static inline float sqr(float x)
{
    return x*x;
}
static float rect_dist_sq(kd_node *node, float x, float y)
{
    float dx = 0.0f, dy = 0.0f;
    if (x < node->minx)
        dx = node->minx - x;
    else if (x > node->maxx)
        dx = x - node->maxx;
    if (y < node->miny)
        dy = node->miny - y;
    else if (y > node->maxy)
        dy = y - node->maxy;
    return sqr(dx) + sqr(dy);
}

// helper swap
static void kd_swap(kd_point *a, kd_point *b)
{
    kd_point tmp = *a;
    *a = *b;
    *b = tmp;
}

// nth_element partition (quickselect) to choose median by axis
static int kd_partition(kd_point *arr, int lo, int hi, int pivot, int axis)
{
    float pivotVal = (axis == 0) ? arr[pivot].x : arr[pivot].y;
    kd_swap(&arr[pivot], &arr[hi]);
    int store = lo;
    for (int i = lo; i < hi; ++i) {
        float v = (axis == 0) ? arr[i].x : arr[i].y;
        if (v < pivotVal || (v == pivotVal && arr[i].index < arr[hi].index)) {
            kd_swap(&arr[i], &arr[store]);
            store++;
        }
    }
    kd_swap(&arr[store], &arr[hi]);
    return store;
}

static void kd_nth_element(kd_point *arr, int lo, int hi, int n, int axis)
{
    while (lo < hi) {
        int pivot = lo + (hi - lo) / 2;
        int pivotNew = kd_partition(arr, lo, hi, pivot, axis);
        if (n == pivotNew)
            return;
        else if (n < pivotNew)
            hi = pivotNew - 1;
        else
            lo = pivotNew + 1;
    }
}

// Build kd-tree recursively, returns root (allocates nodes)
static kd_node* kd_build_rec(kd_point *pts, int lo, int hi, int depth)
{
    if (lo > hi)
        return NULL;
    int axis = depth % 2;
    int mid = (lo + hi) / 2;
    kd_nth_element(pts, lo, hi, mid, axis);

    kd_node *node = (kd_node*)malloc(sizeof(kd_node));
    node->pt = pts[mid];
    node->axis = axis;
    node->left = node->right = NULL;
    node->minx = node->maxx = node->pt.x;
    node->miny = node->maxy = node->pt.y;

    node->left = kd_build_rec(pts, lo, mid - 1, depth + 1);
    node->right = kd_build_rec(pts, mid + 1, hi, depth + 1);

    // compute bounding box from children
    if (node->left) {
        node->minx = fminf(node->minx, node->left->minx);
        node->miny = fminf(node->miny, node->left->miny);
        node->maxx = fmaxf(node->maxx, node->left->maxx);
        node->maxy = fmaxf(node->maxy, node->left->maxy);
    }
    if (node->right) {
        node->minx = fminf(node->minx, node->right->minx);
        node->miny = fminf(node->miny, node->right->miny);
        node->maxx = fmaxf(node->maxx, node->right->maxx);
        node->maxy = fmaxf(node->maxy, node->right->maxy);
    }
    return node;
}

// helper to safely push unique indices; we rely on a mark array sized to number of points.
// caller must provide mark array of size n_points (initialized to 0).
static void kd_query_radius_collect(kd_node *node, float qx, float qy, float r, int **dst, int *mark, int n_points)
{
    if (!node)
        return;
    float r2 = r * r;
    if (rect_dist_sq(node, qx, qy) > r2)
        return;
    float dx = node->pt.x - qx;
    float dy = node->pt.y - qy;
    float dist2 = dx*dx + dy*dy;
    if (dist2 <= r2) {
        int idx = node->pt.index;
        if (!mark[idx]) {
            arrpush(*dst, idx);
            mark[idx] = 1;
        }
    }
    if (node->left)
        kd_query_radius_collect(node->left, qx, qy, r, dst, mark, n_points);
    if (node->right)
        kd_query_radius_collect(node->right, qx, qy, r, dst, mark, n_points);
}

// recursive radius search, fills dst (stb_ds int array)
static void kd_query_radius_rec(kd_node *node, float qx, float qy, float r, int **dst)
{
    if (!node)
        return;
    float r2 = r * r;
    // prune if bounding box is farther than r
    if (rect_dist_sq(node, qx, qy) > r2)
        return;

    float dx = node->pt.x - qx;
    float dy = node->pt.y - qy;
    float dist2 = dx*dx + dy*dy;
    if (dist2 <= r2) {
        arrpush(*dst, node->pt.index);
    }

    // search children - choose order heuristically
    if (node->left && node->right) {
        float dl = rect_dist_sq(node->left, qx, qy);
        float dr = rect_dist_sq(node->right, qx, qy);
        if (dl < dr) {
            if (dl <= r2)
                kd_query_radius_rec(node->left, qx, qy, r, dst);
            if (dr <= r2)
                kd_query_radius_rec(node->right, qx, qy, r, dst);
        } else {
            if (dr <= r2)
                kd_query_radius_rec(node->right, qx, qy, r, dst);
            if (dl <= r2)
                kd_query_radius_rec(node->left, qx, qy, r, dst);
        }
    } else {
        if (node->left && rect_dist_sq(node->left, qx, qy) <= r2)
            kd_query_radius_rec(node->left, qx, qy, r, dst);
        if (node->right && rect_dist_sq(node->right, qx, qy) <= r2)
            kd_query_radius_rec(node->right, qx, qy, r, dst);
    }
}

static int kd_query_radius_periodic_alloc(kd_node *tree, float px, float py, float radius, int **out_indices, int n_points)
{
    if (!out_indices)
        return 0;
    *out_indices = NULL;
    if (!tree)
        return 0;
    int *mark = (int*)calloc(n_points, sizeof(int));
    if (!mark)
        return 0;
    for (int ox = -1; ox <= 1; ++ox) {
        for (int oy = -1; oy <= 1; ++oy) {
            float qx = px + ox * WORLD_SIZE;
            float qy = py + oy * WORLD_SIZE;
            kd_query_radius_collect(tree, qx, qy, radius, out_indices, mark, n_points);
        }
    }
    int cnt = arrlen(*out_indices);
    free(mark); // free mark but keep out_indices for caller
    return cnt;
}

static int kd_query_radius_periodic_with_mark(kd_node *tree, float px, float py, float radius, int **out_indices, int *mark, int n_points)
{
    if (!out_indices || !mark)
        return 0;
    *out_indices = NULL; // important: ensure dst starts NULL for stb_ds
    if (!tree)
        return 0;
    memset(mark, 0, sizeof(int) * n_points);
    for (int ox = -1; ox <= 1; ++ox) {
        for (int oy = -1; oy <= 1; ++oy) {
            float qx = px + ox * WORLD_SIZE;
            float qy = py + oy * WORLD_SIZE;
            kd_query_radius_collect(tree, qx, qy, radius, out_indices, mark, n_points);
        }
    }
    return arrlen(*out_indices);
}

// public build: copies points into temp array (because build modifies order)
kd_node *kd_build(const kd_point *points, int n)
{
    if (n <= 0)
        return NULL;
    kd_point *copy = (kd_point*)malloc(sizeof(kd_point) * n);
    memcpy(copy, points, sizeof(kd_point) * n);
    kd_node* root = kd_build_rec(copy, 0, n - 1, 0);
    free(copy);
    return root;
}

// free tree
void kd_free(kd_node *node)
{
    if (!node)
        return;
    kd_free(node->left);
    kd_free(node->right);
    free(node);
}

// periodic-aware public query
int kd_query_radius_periodic(kd_node *tree, float px, float py, float radius, int **out_indices, int n_points)
{
    if (!out_indices)
        return 0;
    *out_indices = NULL;
    if (!tree)
        return 0;
    // allocate temp mark array on heap (size = number of points)
    int *mark = (int*)calloc(n_points, sizeof(int));
    if (!mark)
        return 0;

    // offsets to consider: -1, 0, +1 in each axis (world repeats)
    for (int ox = -1; ox <= 1; ++ox) {
        for (int oy = -1; oy <= 1; ++oy) {
            float qx = px + ox * WORLD_SIZE;
            float qy = py + oy * WORLD_SIZE;
            kd_query_radius_rec(tree, qx, qy, radius, out_indices);
            // kd_query_radius_rec doesn't check mark, so instead call collector:
            // (call collector instead of above if you replaced kd_query_radius_rec usage)
        }
    }

    // However kd_query_radius_rec pushes raw indices possibly duplicated; instead we should use collector:
    // Clear out_indices then use collector:
    arrsetlen(*out_indices, 0);
    for (int ox = -1; ox <= 1; ++ox) {
        for (int oy = -1; oy <= 1; ++oy) {
            float qx = px + ox * WORLD_SIZE;
            float qy = py + oy * WORLD_SIZE;
            kd_query_radius_collect(tree, qx, qy, radius, out_indices, mark, n_points);
        }
    }

    int cnt = arrlen(*out_indices);
    free(mark);
    return cnt;
}

// public query
int kd_query_radius(kd_node *tree, float px, float py, float radius, int **out_indices)
{
    if (!out_indices)
        return 0;
    *out_indices = NULL;
    if (!tree)
        return 0;
    if (USE_PERIODIC) {
        if (!CURRENT_SIM_INSTANCE || !CURRENT_SIM_INSTANCE->neighbor_mark) {
            // fallback to alloc-based (safety)
//            return kd_query_radius_periodic(tree, px, py, radius, out_indices, KD_CURRENT_POINT_COUNT);
            return kd_query_radius_periodic_alloc(tree, px, py, radius, out_indices, KD_CURRENT_POINT_COUNT);
        }
        return kd_query_radius_periodic_with_mark(tree, px, py, radius, out_indices, CURRENT_SIM_INSTANCE->neighbor_mark, KD_CURRENT_POINT_COUNT);
    } else {
        kd_query_radius_rec(tree, px, py, radius, out_indices);
        return arrlen(*out_indices);
    }
}

// draw kd tree partitions in world coords using mock draw functions.
// draw_line(x1,y1,x2,y2) must be implemented externally (mock).
void kd_draw_tree_partitions(kd_node *node)
{
    if (!node)
        return;
    // draw splitting line across node bbox
    if (node->axis == 0) {
        // vertical line at x = node->pt.x, spanning miny..maxy
        draw_line_world(node->pt.x, node->miny, node->pt.x, node->maxy);
    } else {
        // horizontal line at y = node->pt.y, spanning minx..maxx
        draw_line_world(node->minx, node->pt.y, node->maxx, node->pt.y);
    }
    kd_draw_tree_partitions(node->left);
    kd_draw_tree_partitions(node->right);
}

// draw kd-tree partitions and node points.
// node: current node
// depth: recursion depth (used for color/size if desired)
void kd_draw_tree_partitions_full(kd_node *node, int depth)
{
    if (!node)
        return;

    // draw splitting line across node bbox
    if (node->axis == 0) {
        // vertical at x = node->pt.x, from miny..maxy
        float x = node->pt.x;
        float y0 = node->miny;
        float y1 = node->maxy;
        // if periodic, draw also shifted copies by +-WORLD_SIZE
        if (USE_PERIODIC) {
            for (int sx = -1; sx <= 1; ++sx) {
                draw_line_world(x + sx, y0, x + sx, y1);
            }
        } else {
            draw_line_world(x, y0, x, y1);
        }
    } else {
        // horizontal at y = node->pt.y, from minx..maxx
        float y = node->pt.y;
        float x0 = node->minx;
        float x1 = node->maxx;
        if (USE_PERIODIC) {
            for (int sy = -1; sy <= 1; ++sy) {
                draw_line_world(x0, y + sy, x1, y + sy);
            }
        } else {
            draw_line_world(x0, y, x1, y);
        }
    }

    // draw bbox edges
    // draw rectangle (minx,miny) -> (maxx,miny) -> (maxx,maxy) -> (minx,maxy)
    if (USE_PERIODIC) {
        for (int ox = -1; ox <= 1; ++ox) {
            for (int oy = -1; oy <= 1; ++oy) {
                if (!(abs(ox) == 1 && abs(oy) == 1)) {
                    float mx = node->minx + ox;
                    float my = node->miny + oy;
                    float Mx = node->maxx + ox;
                    float My = node->maxy + oy;
                    draw_line_world(mx, my, Mx, my);
                    draw_line_world(Mx, my, Mx, My);
                    draw_line_world(Mx, My, mx, My);
                    draw_line_world(mx, My, mx, my);
                }
            }
        }
    } else {
        draw_line_world(node->minx, node->miny, node->maxx, node->miny);
        draw_line_world(node->maxx, node->miny, node->maxx, node->maxy);
        draw_line_world(node->maxx, node->maxy, node->minx, node->maxy);
        draw_line_world(node->minx, node->maxy, node->minx, node->miny);
    }

    // recurse
    kd_draw_tree_partitions_full(node->left, depth + 1);
    kd_draw_tree_partitions_full(node->right, depth + 1);
}

void kd_serialize(kd_node* tree, const char* filename)
{
    if (!tree || !filename)
        return;
    FILE *f = fopen(filename, "w");
    if (!f)
        return;
    // simple preorder write: index x y axis minx miny maxx maxy
    // recursion:
    void rec(kd_node* n) {
        if (!n)
            return;
        fprintf(f, "%d %.6f %.6f %d %.6f %.6f %.6f %.6f\n",
                n->pt.index,
                n->pt.x, n->pt.y,
                n->axis,
                n->minx, n->miny,
                n->maxx, n->maxy);
        rec(n->left);
        rec(n->right);
    }
    rec(tree);
    fclose(f);
}
