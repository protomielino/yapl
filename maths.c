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

float map(float input, float input_start, float input_end, float output_start, float output_end)
{
    float slope = 1.0 * (output_end - output_start) / (input_end - input_start);
    float output = output_start + slope * (input - input_start);
    return output;
}

// clamp wrap into [0,1)
float wrap01(float v)
{
    if (v >= 1.0f) v -= floorf(v);
    if (v < 0.0f) v -= floorf(v);
    return v;
}

// --- xorshift64* PRNG ---
static uint64_t rng_state;

void rand_seed(unsigned int seed)
{
    rng_state = seed;
    if (rng_state == 0) rng_state = 1;
    // warmup to avoid seed correlation
    for (int i = 0; i < 4; i++) {
        uint64_t x = rng_state;
        x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
        rng_state = x;
    }
}

float rand01(void)
{
    uint64_t x = rng_state;
    x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
    rng_state = x;
    return (float)((x * 0x2545F4914F6CDD1DULL) >> 40) / 16777216.0f;
}

float rand_range(float min, float max)
{
    return min + rand01() * (max - min);
}

