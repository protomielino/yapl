#include <stdlib.h>
#include "maths.h"

// pan-zoom variables
extern Vector2 offset;
extern Vector2 scale;

// Convert coordinates from World Space --> Screen Space
Vector2 WorldToScreen(Vector2 world)
{
    Vector2 screen = Vector2Multiply(Vector2Subtract(world, offset), scale);
    screen = (Vector2) { screen.x, screen.y };
    return screen;
}

// Convert coordinates from Screen Space --> World Space
Vector2 ScreenToWorld(Vector2 screen)
{
    Vector2 scr = { (float)screen.x, (float)screen.y};
    Vector2 world = Vector2Add(Vector2Divide(scr, scale), offset);
    return world;
}

float function (float x)
{
    return cos(x);
}

float map(float input, float input_start, float input_end, float output_start, float output_end)
{
    float slope = 1.0 * (output_end - output_start) / (input_end - input_start);
    float output = output_start + slope * (input - input_start);
    return output;
}

float random_range(float min, float max)
{
    return drand48() * (max - min) + min;
}

// helper random
extern unsigned int rnd_state;
float frand()
{
    rnd_state = rnd_state * 1664525u + 1013904223u;
    return (rnd_state & 0x00FFFFFF) / (float)0x01000000;
}

// clamp wrap into [0,1)
float wrap01(float v)
{
    if (v >= 1.0f) v -= floorf(v);
    if (v < 0.0f) v -= floorf(v);
    return v;
}

float randf()
{
    return (float)rand() / (float)RAND_MAX;
}

