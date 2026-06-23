#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <raylib.h>
#define RAYMATH_IMPLEMENTATION
#include <raymath.h>

#include "maths.h"
#include "sim.h"

const int WIDTH  = 700;
const int HEIGHT = 700;

Vector2 mouse_screen;
Vector2 mouse_world;

Vector2 offset;
Vector2 scale;

Vector2 startPan;

int bx;
int by, ty, dy;
int textDim;
char text[1024] = {};

size_t frameCounter;

bool show_help;
bool show_HUD;
bool show_tree;

bool showForcesDBG;
bool showMinDistancesDBG;
bool showMassesDBG;
bool showRadiiDBG;

float cell_size;
int cell_row_hover;
int cell_col_hover;
int cell_row_selected;
int cell_col_selected;
bool mouse_over_matrix;

// -------------------- configuration --------------------
const int DEFAULT_N = 1000;
const float DEFAULT_DT = 0.01f;
const float DEFAULT_FRICTION_HALF_LIFE = 0.040f;
const float DEFAULT_RMAX = 0.1f;
const int DEFAULT_M = 6;
const float FORCE_FACTOR = 10.0f;

const int USE_PERIODIC = 1;
const float WORLD_SIZE = 1.0f;

void showHUD(sim *s)
{
    if(show_HUD) {
        sprintf(text, "counter: %ld", frameCounter);
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "offset: {%d, %d} [CLICK+DRAG]", (int)offset.x, (int)offset.y);
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "scale: {%.2f, %.2f} [WHEEL]", scale.x, scale.y);
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        const char *backend_name =
            s->backend == BACKEND_GRID ? "grid" :
            s->backend == BACKEND_KD ? "k-d tree" : "quadtree";
        sprintf(text, "backend: %s", backend_name);
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "force-law: %s", force_law_names[s->force_law]);
        int law_label_w = MeasureText(text, textDim);
        DrawText(text, 10, ty, textDim, WHITE);

        { // small force-law plot
            int gw = 64, gh = 36;
            int gx = 10 + law_label_w + 8;
            int gy = ty;
            float a_disp = 1.0f, beta_disp = 0.3f;
            int ns = 60;

            // sample to find y range
            float ymin = 0, ymax = 0;
            for (int i = 0; i <= ns; ++i) {
                float r = (float)i / ns;
                float f = sim_eval_force_law(s->force_law, r, a_disp, beta_disp);
                f = fminf(fmaxf(f, -3.0f), 3.0f);
                if (f < ymin) ymin = f;
                if (f > ymax) ymax = f;
            }
            float range = ymax - ymin;
            if (range < 0.01f) range = 1.0f;

            DrawRectangle(gx-1, gy-1, gw+2, gh+2, (Color){20,20,20,220});
            DrawRectangleLines(gx, gy, gw, gh, (Color){80,80,80,180});

            // zero line
            if (ymin < 0 && ymax > 0) {
                float zy = gy + gh - (0 - ymin) / range * gh;
                DrawLine(gx, (int)zy, gx + gw - 1, (int)zy, (Color){60,60,60,100});
            }

            // curve segments
            for (int i = 0; i < ns; ++i) {
                float r1 = (float)i / ns, r2 = (float)(i+1) / ns;
                float f1 = sim_eval_force_law(s->force_law, r1, a_disp, beta_disp);
                float f2 = sim_eval_force_law(s->force_law, r2, a_disp, beta_disp);
                f1 = fminf(fmaxf(f1, -3.0f), 3.0f);
                f2 = fminf(fmaxf(f2, -3.0f), 3.0f);
                int x1 = gx + (int)(r1 * gw);
                int y1 = gy + gh - (int)((f1 - ymin) / range * gh);
                int x2 = gx + (int)(r2 * gw);
                int y2 = gy + gh - (int)((f2 - ymin) / range * gh);
                y1 = (y1 < gy) ? gy : (y1 >= gy + gh ? gy + gh - 1 : y1);
                y2 = (y2 < gy) ? gy : (y2 >= gy + gh ? gy + gh - 1 : y2);
                DrawLine(x1, y1, x2, y2, WHITE);
            }
        }
        ty += dy;
        ty += dy;
    }
}

void showHelp()
{
    if(show_help) {
        sprintf(text, "[R] reset pan & zoom");
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "[MMB] reset pan");
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "[d] show/hide tree");
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "[p] pause / unpause");
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "[r] randomize positions");
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "[SPACE] randomize all");
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "[z,x,c,v] forces, minDist, masses, radii");
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "[f] cycle force law");
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "[scroll] edit cell under mouse");
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "[CTRL+click] select cell");
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        sprintf(text, "[CTRL+right-click] zero active panel");
        DrawText(text, 10, ty, textDim, WHITE);
        ty += dy;
        ty += dy;
    }
}

void cycle_force_law(sim *s)
{
    s->force_law = (force_law_type)((s->force_law + 1) % FORCE_LAW_COUNT);
}

void processInputs(sim *s)
{
    if (IsKeyPressed(KEY_F)) {
        cycle_force_law(s);
    }
    if(IsKeyDown(KEY_R) && IsKeyDown(KEY_LEFT_SHIFT)) {
        offset = Vector2Zero();
        offset = Vector2Divide(offset, scale);
        scale = Vector2One();
    }
    if (IsMouseButtonPressed(MOUSE_BUTTON_MIDDLE)) {
        offset = Vector2Zero();
        offset = Vector2Divide(offset, scale);
    }
    if(IsKeyPressed(KEY_D)) {
        show_tree = !show_tree;
    }
    if(IsKeyPressed(KEY_H) && IsKeyDown(KEY_LEFT_SHIFT)) {
        show_HUD = !show_HUD;
    }
    if(IsKeyPressed(KEY_H) && !IsKeyDown(KEY_LEFT_SHIFT)) {
        show_help = !show_help;
    }
    if(IsKeyPressed(KEY_Z)) {
        showForcesDBG = !showForcesDBG;
        showMinDistancesDBG = false;
        showMassesDBG = false;
        showRadiiDBG = false;
    }
    if(IsKeyPressed(KEY_X)) {
        showForcesDBG = false;
        showMinDistancesDBG = !showMinDistancesDBG;
        showMassesDBG = false;
        showRadiiDBG = false;
    }
    if(IsKeyPressed(KEY_C)) {
        showForcesDBG = false;
        showMinDistancesDBG = false;
        showMassesDBG = !showMassesDBG;
        showRadiiDBG = false;
    }
    if(IsKeyPressed(KEY_V)) {
        showForcesDBG = false;
        showMinDistancesDBG = false;
        showMassesDBG = false;
        showRadiiDBG = !showRadiiDBG;
    }
}

void updatePanZoom()
{
    mouse_screen = (Vector2){ (float)GetMouseX(), (float)GetMouseY() };
    mouse_world = ScreenToWorld(mouse_screen);

    if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
        startPan = mouse_screen;
    }
    if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
        offset = Vector2Subtract(
                offset,
                Vector2Divide(
                        Vector2Subtract(mouse_screen, startPan),
                                scale));
        startPan = mouse_screen;
    }
}

static void zoom_at_cursor(float wheel)
{
    if (wheel == 0.0f) return;
    Vector2 mouseWorld_BeforeZoom = ScreenToWorld(mouse_screen);

    if (wheel > 0) {
        scale.x *= 1.1f;
        scale.y *= 1.1f;
    } else {
        scale.x *= 0.9f;
        scale.y *= 0.9f;
    }

    Vector2 mouseWorld_AfterZoom = ScreenToWorld(mouse_screen);
    offset = Vector2Add(
            offset,
            Vector2Subtract(
                    mouseWorld_BeforeZoom,
                    mouseWorld_AfterZoom));
}

void initPanZoom(Vector2 pan, Vector2 sca)
{
    startPan = pan;
    offset = Vector2Zero();
    scale = sca;
}

void initHUD()
{
    textDim = 20;
    by = 50;
    ty = by;
    dy = textDim;
    bx = 10;
    cell_size = 50.0f;
    cell_row_hover = 0;
    cell_col_hover = 0;
    cell_row_selected = -1;
    cell_col_selected = -1;
    mouse_over_matrix = false;
}

// -------------------- drawing API --------------------
void draw_line_world(float x1, float y1, float x2, float y2)
{
    Vector2 v1 = { x1 * WIDTH, y1 * HEIGHT };
    Vector2 v2 = { x2 * WIDTH, y2 * HEIGHT };
    v1 = WorldToScreen(v1);
    v2 = WorldToScreen(v2);
    Color c = WHITE; c.a = 32;
    DrawLine(v1.x, v1.y, v2.x, v2.y, c);
}

void draw_particle_world(float x, float y, int color, float mass)
{
    Color c = ColorFromHSV(360.0f / DEFAULT_M * color, 1.0f, 1.0f);
    Vector2 v = { x * WIDTH, y * HEIGHT };
    v = WorldToScreen(v);
    float radius = fmaxf(0.5, mass * 1.5);
    radius *= scale.x;
    radius = radius >= 1.0f ? radius * 0.75f : 0.75f;
    DrawCircleLines(v.x, v.y, radius, c);
}

// -------------------- matrix panel drawing --------------------
static void draw_color_header(int type, float x, float y, float size)
{
    Color hc = ColorFromHSV(360.0f / DEFAULT_M * type, 0.8f, 1.0f);
    DrawRectangleV((Vector2){x, y}, (Vector2){size, size}, hc);
}

static void draw_matrix_cell(float x, float y, float size, float val, float range_min, float range_max)
{
    float t = (val - range_min) / (range_max - range_min); // 0..1
    Color c;
    if (val <= 0.0f) {
        c = ColorFromHSV(360.0f, 1.0f, 1.0f - t);
    } else {
        c = ColorFromHSV(360.0f / 3.0f, 1.0f, t);
    }
    if (val == 0.0f) c = BLACK;
    c.a = 128;
    DrawRectangleV((Vector2){x, y}, (Vector2){size, size}, c);
}

static void draw_cell_border(float x, float y, float size, int hover, int selected)
{
    if (hover) {
        DrawRectangleLinesEx((Rectangle){x-1, y-1, size+2, size+2}, 1, WHITE);
    }
    if (selected) {
        DrawRectangleLinesEx((Rectangle){x-2, y-2, size+4, size+4}, 2, YELLOW);
    }
}

static void draw_cell_text(float x, float y, float size, const char *label)
{
    DrawText(label,
        x + size / 2.0f - textDim / 2.0f,
        y + size / 2.0f - textDim / 3.0f,
        textDim, WHITE);
}

static void draw_matrix_panel(const char *title, float *matrix, int m,
                               float range_min, float range_max, int is_mass_panel,
                               sim *s)
{
    ty += dy;
    sprintf(text, "%s:", title);
    DrawText(text, 10, ty, textDim, WHITE);
    ty += dy;
    sprintf(text, "min: %.2f", range_min);
    DrawText(text, 10, ty, textDim, WHITE);
    ty += dy;
    sprintf(text, "max: %.2f", range_max);
    DrawText(text, 10, ty, textDim, WHITE);
    ty += dy;
    ty += dy;

    // column headers
    for (int col = 0; col < m; ++col) {
        draw_color_header(col, bx + cell_size + (cell_size + 2) * col + 5, ty, cell_size);
    }

    if (!is_mass_panel) {
        // row headers + grid
        for (int row = 0; row < m; ++row) {
            float ry = ty + cell_size + (cell_size + 2) * row + 5;
            draw_color_header(row, bx, ry, cell_size);
        }
    }

    // cells
    int nrows = is_mass_panel ? 1 : m;
    for (int row = 0; row < nrows; ++row) {
        for (int col = 0; col < m; ++col) {
            float cx = bx + cell_size + (cell_size + 2) * col + 5;
            float cy = ty + cell_size + (cell_size + 2) * row + 5;

            float val;
            if (is_mass_panel) {
                // compute average mass for this type from particle data
                float sum = 0;
                int cnt = 0;
                for (int i = 0; i < s->n; ++i) {
                    if (s->colors[i] == col) {
                        sum += s->masses[i];
                        cnt++;
                    }
                }
                val = cnt > 0 ? sum / cnt : 0.0f;
            } else {
                val = matrix[row * m + col];
            }

            draw_matrix_cell(cx, cy, cell_size, val, range_min, range_max);

            int hover = 0, sel = 0;
            if (!is_mass_panel) {
                hover = (cell_row_hover == row && cell_col_hover == col);
                sel = (cell_row_selected == row && cell_col_selected == col);
            } else {
                hover = (cell_col_hover == col);
                sel = (cell_col_selected == col);
            }
            draw_cell_border(cx, cy, cell_size, hover, sel);

            sprintf(text, "%.2f", val);
            draw_cell_text(cx, cy, cell_size, text);
        }
    }
}

static void handle_mouse_hover(int ty_now)
{
    Vector2 mp = GetMousePosition();
    int header_h = 5 * dy;
    int nrows = showMassesDBG ? 1 : DEFAULT_M;
    int raw_col = (int)((mp.x - bx - cell_size - 5) / (cell_size + 2));
    int raw_row = (int)((mp.y - ty_now - header_h - cell_size - 5) / (cell_size + 2));

    mouse_over_matrix = (raw_col >= 0 && raw_col < DEFAULT_M && raw_row >= 0 && raw_row < nrows);

    cell_col_hover = (raw_col < 0) ? 0 : (raw_col >= DEFAULT_M ? DEFAULT_M - 1 : raw_col);
    cell_row_hover = (raw_row < 0) ? 0 : (raw_row >= nrows ? nrows - 1 : raw_row);
}

static void handle_cell_edit(sim *s, float delta)
{
    if (showForcesDBG && cell_row_selected >= 0 && cell_col_selected >= 0) {
        int idx = cell_row_selected * s->m + cell_col_selected;
        s->matrix[idx] += delta;
        if (s->matrix[idx] < -1.0f) s->matrix[idx] = -1.0f;
        if (s->matrix[idx] > 1.0f) s->matrix[idx] = 1.0f;
    }
    if (showMinDistancesDBG && cell_row_selected >= 0 && cell_col_selected >= 0) {
        int idx = cell_row_selected * s->m + cell_col_selected;
        s->min_dist_matrix[idx] += delta;
        if (s->min_dist_matrix[idx] < 0.0f) s->min_dist_matrix[idx] = 0.0f;
        if (s->min_dist_matrix[idx] > 1.0f) s->min_dist_matrix[idx] = 1.0f;
    }
    if (showMassesDBG && cell_col_selected >= 0) {
        for (int i = 0; i < s->n; ++i) {
            if (s->colors[i] == cell_col_selected) {
                s->masses[i] += delta;
                if (s->masses[i] < 0.01f) s->masses[i] = 0.01f;
            }
        }
    }
    if (showRadiiDBG && cell_row_selected >= 0 && cell_col_selected >= 0) {
        int idx = cell_row_selected * s->m + cell_col_selected;
        s->radii_matrix[idx] += delta;
        if (s->radii_matrix[idx] < 0.0f) s->radii_matrix[idx] = 0.0f;
        if (s->radii_matrix[idx] > 1.0f) s->radii_matrix[idx] = 1.0f;
    }
}

static void zero_active_panel(sim *s)
{
    if (showForcesDBG) {
        memset(s->matrix, 0, sizeof(float) * s->m * s->m);
    }
    if (showMinDistancesDBG) {
        memset(s->min_dist_matrix, 0, sizeof(float) * s->m * s->m);
    }
    if (showMassesDBG) {
        for (int i = 0; i < s->n; ++i) s->masses[i] = 1.0f;
    }
    if (showRadiiDBG) {
        memset(s->radii_matrix, 0, sizeof(float) * s->m * s->m);
    }
}

int main(int argc, char *argv[])
{
    int n = DEFAULT_N;
    int m_colors = DEFAULT_M;
    backend_type backend = BACKEND_GRID;
    force_law_type force_law = FORCE_LAW_STANDARD;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-n") == 0 && i + 1 < argc) {
            n = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-m") == 0 && i + 1 < argc) {
            m_colors = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--backend") == 0 && i + 1 < argc) {
            i++;
            if (strcmp(argv[i], "kd") == 0) backend = BACKEND_KD;
            else if (strcmp(argv[i], "grid") == 0) backend = BACKEND_GRID;
            else if (strcmp(argv[i], "qt") == 0) backend = BACKEND_QUADTREE;
            else { fprintf(stderr, "Unknown backend: %s (use grid, kd, or qt)\n", argv[i]); return 1; }
        } else if (strcmp(argv[i], "--force-law") == 0 && i + 1 < argc) {
            i++;
            int found = 0;
            if (strcmp(argv[i], "lj") == 0) {
                force_law = FORCE_LAW_LENNARD_JONES;
                found = 1;
            }
            if (!found) {
                for (int fl = 0; fl < FORCE_LAW_COUNT; ++fl) {
                    if (strcmp(argv[i], force_law_names[fl]) == 0) {
                        force_law = (force_law_type)fl;
                        found = 1;
                        break;
                    }
                }
            }
            if (!found) {
                fprintf(stderr, "Unknown force law: %s (use", argv[i]);
                for (int fl = 0; fl < FORCE_LAW_COUNT; ++fl)
                    fprintf(stderr, "%s %s", fl == 0 ? "" : ", ", force_law_names[fl]);
                fprintf(stderr, ")\n");
                return 1;
            }
        } else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            printf("Usage: yapl [-n N] [-m M] [--backend grid|kd|qt] [--force-law standard|linear|lj|lennard-jones|smooth|damped-wave]\n");
            return 0;
        } else {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            return 1;
        }
    }

    rand_seed((unsigned)clock());
    frameCounter = 0;

    show_help = false;
    show_HUD = false;
    show_tree = false;
    showForcesDBG = false;
    showMinDistancesDBG = false;
    showMassesDBG = false;
    showRadiiDBG = false;

    sim *s = sim_create(n, m_colors, backend, force_law, DEFAULT_DT, DEFAULT_FRICTION_HALF_LIFE, DEFAULT_RMAX, FORCE_FACTOR);

    initPanZoom(Vector2Zero(), Vector2One());
    initHUD();

    InitWindow(WIDTH, HEIGHT, "Particle Life");
    SetTargetFPS(30);

    while (!WindowShouldClose())
    {
        updatePanZoom();
        processInputs(s);

        if (IsKeyPressed(KEY_SPACE)) {
            sim_randomize_all(s, 0.5f, 2.0f);
        }
        if (IsKeyPressed(KEY_R) && !IsKeyDown(KEY_LEFT_SHIFT)) {
            sim_randomize_positions(s);
        }

        if (IsKeyPressed(KEY_P)) {
            s->dt = s->dt == DEFAULT_DT ? 0.0f : DEFAULT_DT;
        }

        // Ctrl+LeftClick to select cell
        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT) && IsKeyDown(KEY_LEFT_CONTROL)) {
            cell_row_selected = cell_row_hover;
            cell_col_selected = cell_col_hover;
        }

        // Ctrl+RightClick to zero active panel
        if (IsMouseButtonPressed(MOUSE_BUTTON_RIGHT) && IsKeyDown(KEY_LEFT_CONTROL)) {
            zero_active_panel(s);
        }

        frameCounter++;

        BeginDrawing(); {
            ClearBackground(BLACK);

            sim_update(s);
            sim_draw_frame(s);

            Vector2 cpos = Vector2Zero();
            DrawLineV(
                    WorldToScreen((Vector2){cpos.x, 0.0f}),
                    WorldToScreen((Vector2){cpos.x, HEIGHT / 2}), GREEN);
            DrawLineV(
                    WorldToScreen((Vector2){0.0f, cpos.y}),
                    WorldToScreen((Vector2){WIDTH / 2, cpos.y}), RED);
            DrawCircleLinesV(WorldToScreen(cpos), 2.5f * Vector2Length(scale), WHITE);

            sprintf(text, "[0, 0]");
            cpos = WorldToScreen(cpos);
            DrawText(text, cpos.x + 5, cpos.y + 5, textDim / 3 * Vector2Length(scale), WHITE);

            ty = by;
            showHUD(s);
            showHelp();
            sprintf(text, "help, HUD: [h,H]");
            DrawText(text, 10, ty, textDim, WHITE);
            ty += dy;

            // --- debug panels ---
            bool any_panel = showForcesDBG || showMinDistancesDBG || showMassesDBG || showRadiiDBG;
            if (any_panel) {
                handle_mouse_hover(ty);
            }

            // Wheel: edit hovered cell or zoom at cursor
            float wheel = GetMouseWheelMove();
            if (any_panel && mouse_over_matrix) {
                float step = 0.05f;
                if (IsKeyDown(KEY_LEFT_SHIFT)) step = 0.25f;
                if (IsKeyDown(KEY_LEFT_CONTROL)) step = 0.01f;
                if (wheel > 0) {
                    cell_row_selected = cell_row_hover;
                    cell_col_selected = cell_col_hover;
                    handle_cell_edit(s, step);
                } else if (wheel < 0) {
                    cell_row_selected = cell_row_hover;
                    cell_col_selected = cell_col_hover;
                    handle_cell_edit(s, -step);
                }
            } else if (wheel != 0.0f) {
                zoom_at_cursor(wheel);
            }

            if (showForcesDBG) {
                draw_matrix_panel("forces", s->matrix, s->m, -1.0f, 1.0f, 0, s);
            }
            if (showMinDistancesDBG) {
                draw_matrix_panel("minDist", s->min_dist_matrix, s->m, 0.0f, 0.35f, 0, s);
            }
            if (showMassesDBG) {
                draw_matrix_panel("masses (avg)", NULL, s->m, 0.0f, 2.5f, 1, s);
            }
            if (showRadiiDBG) {
                draw_matrix_panel("radii", s->radii_matrix, s->m, 0.3f, 1.0f, 0, s);
            }

            DrawFPS(20, 20);
        } EndDrawing();
    }

    sim_free(s);
    CloseWindow();
    return 0;
}
