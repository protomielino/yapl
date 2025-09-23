#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <raylib.h>

typedef struct
{
    Vector2 position;
    Vector2 velocity;
    float mass;
    int type;
    Vector2 newPosition;
} particle;

#endif /* PARTICLE_H_ */
