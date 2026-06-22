#ifndef MATHS_H_
#define MATHS_H_

#include <stdint.h>
#include <raylib.h>
#include <raymath.h>

// Convert coordinates from World to Screen Space
Vector2 WorldToScreen(Vector2 world);
// Convert coordinates from Screen to World Space
Vector2 ScreenToWorld(Vector2 screen);
float map(float input, float input_start, float input_end, float output_start, float output_end);
void rand_seed(unsigned int seed);
float rand01(void);

#endif /* MATHS_H_ */
