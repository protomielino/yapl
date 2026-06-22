#ifndef GRID_H_
#define GRID_H_

typedef struct grid_s {
    int ncells;
    float inv_cell_size;
    int **cells;
} grid;

grid* grid_create(float cell_size, float world_size);
void grid_free(grid *g);
void grid_clear(grid *g);
void grid_insert(grid *g, int idx, float x, float y);
int grid_query(grid *g, float x, float y, int **out_indices);
void grid_draw_partitions(grid *g);

#endif /* GRID_H_ */
