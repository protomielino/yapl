#ifndef KD_H_
#define KD_H_

typedef struct kd_point_s
{
    float x, y;
    int index; // original particle index
} kd_point;

typedef struct kd_node_s
{
    kd_point pt;
    int axis; // 0 = x, 1 = y
    struct kd_node_s *left;
    struct kd_node_s *right;
    // bounding box for possible pruning (not strictly needed but helpful)
    float minx, miny, maxx, maxy;
} kd_node;

kd_node *kd_build(const kd_point *points, int n);
void kd_free(kd_node *node);
int kd_query_radius_periodic(kd_node *tree, float px, float py, float radius, int **out_indices, int n_points);
int kd_query_radius(kd_node *tree, float px, float py, float radius, int **out_indices);
void kd_draw_tree_partitions(kd_node *node);
void kd_draw_tree_partitions_full(kd_node *node, int depth);
void kd_serialize(kd_node* tree, const char* filename);

#endif /* KD_H_ */
