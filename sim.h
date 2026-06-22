#ifndef SIM_H_
#define SIM_H_

typedef struct grid_s grid;
typedef struct kd_node_s kd_node;
typedef struct qt_s qt;

typedef struct { float x, y; } vec2;
typedef enum { BACKEND_GRID, BACKEND_KD, BACKEND_QUADTREE } backend_type;

typedef struct sim_s
{
    int n;
    int m;
    backend_type backend;
    float dt;
    float friction_half_life;
    float rmax;
    float frictionFactor;
    float forceFactor;
    int *colors;
    vec2 *positions;
    vec2 *velocities;
    float *matrix;
    float *masses;
    float *min_dist_matrix;
    float *radii_matrix;
    int *neighbor_mark;
    grid *grid;
    struct kd_node_s *kd_tree;
    int kd_rebuild_counter;
    qt *qt;
} sim;

sim* sim_create(int n, int m, backend_type backend, float dt, float friction_half_life, float rmax, float forceFactor);
void sim_free(sim* s);
void sim_update(sim* s);
void sim_draw_frame(sim* s);
void sim_get_positions(sim* s, float* out_x, float* out_y, int* out_colors);
void sim_randomize_matrix(sim *s);
void sim_randomize_masses(sim *s, float min_mass, float max_mass);
void sim_randomize_positions(sim *s);
void sim_randomize_colors(sim *s);
void sim_randomize_all(sim *s, float min_mass, float max_mass);
void sim_randomize_min_dist_matrix(sim *s);
void sim_randomize_radii_matrix(sim *s);

#endif /* SIM_H_ */
