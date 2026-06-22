#ifndef SIM_H_
#define SIM_H_

typedef struct grid_s grid;

typedef struct { float x, y; } vec2;

typedef struct sim_s
{
    int n;
    int m; // number of colors
    float dt;
    float friction_half_life;
    float rmax;
    float frictionFactor;
    float forceFactor;
    // arrays (SoA — Structure of Arrays)
    int *colors;        // size n
    vec2 *positions;    // size n
    vec2 *velocities;   // size n
    // interaction matrix (m x m)
    float *matrix;      // row-major, size m*m
    float *masses;      // massa per particella
    float *min_dist_matrix; // row-major, size m*m, beta per type pair
    float *radii_matrix;    // row-major, size m*m, display only
    int *neighbor_mark; // size n, reused across queries when USE_PERIODIC==1
#if USE_GRID
    grid *grid;
#else
    struct kd_node_s *kd_tree;
    int kd_rebuild_counter;
#endif
} sim;

sim* sim_create(int n, int m, float dt, float friction_half_life, float rmax, float forceFactor);
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
