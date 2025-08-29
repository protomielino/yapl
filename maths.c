#include <stdlib.h>
#include "maths.h"

// pan-zoom variables
extern Vector2 offset;
extern Vector2 scale;

// Convert coordinates from World Space --> Screen Space
Vector2 WorldToScreen(Vector2 world)
{
    Vector2 screen = Vector2Multiply(Vector2Subtract(world, offset), scale);
    screen = (Vector2) { (int)screen.x, (int)screen.y};
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
    double slope = 1.0 * (output_end - output_start) / (input_end - input_start);
    double output = output_start + round(slope * (input - input_start));
    return output;
}

float random_range(float min, float max)
{
    float ret = 0.0f;
    float r = drand48();
    ret = r * (max - min) + min;
    return ret;
}
