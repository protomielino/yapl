#ifndef SIM_H_
#define SIM_H_

typedef struct sim_s
{
    int n;
    int m; // number of colors
    float dt;
    float friction_half_life;
    float rmax;
    float frictionFactor;
    float forceFactor;
    // arrays
    int *colors;        // size n
    float *posx;        // size n
    float *posy;        // size n
    float *velx;        // size n
    float *vely;        // size n
    // interaction matrix (m x m)
    float *matrix;      // row-major, size m*m
    float *masses;      // massa per particella
    int *neighbor_mark; // size n, reused across queries when USE_PERIODIC==1
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

#endif /* SIM_H_ */
