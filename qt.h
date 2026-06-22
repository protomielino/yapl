#ifndef QT_H_
#define QT_H_

typedef struct qt_s qt;

qt* qt_create(float world_size);
void qt_free(qt *q);
void qt_insert(qt *q, int idx, float x, float y);
int qt_query_radius(qt *q, float x, float y, float radius, int *out, int max_out);
void qt_draw_partitions(qt *q);

#endif
