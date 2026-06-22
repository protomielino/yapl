# yapl - Yet Another Particle Life

Particle Life simulation with spatial partitioning for performance.

## Build

```
make          # build both backends
make grid     # grid backend (default, recommended)
make kd       # k-d tree backend
make run      # build grid + run
make run-kd   # build k-d tree + run
make clean
```

## Controls

| Key | Action |
|---|---|
| `[z]` | Show forces matrix (edit with arrows) |
| `[x]` | Show minDist matrix (beta per type pair) |
| `[c]` | Show masses per type |
| `[v]` | Show radii matrix (display only) |
| `[arrows]` | Edit selected cell ±0.1 |
| `[CTRL+click]` | Select cell (yellow border) |
| `[CTRL+right-click]` | Zero active panel |
| `[SPACE]` | Randomize all |
| `[r]` | Randomize positions |
| `[p]` | Pause / unpause |
| `[d]` | Show tree partitions overlay |
| `[h]` | Toggle help |
| `[H]` | Toggle HUD |
| `[R]` | Reset pan & zoom |
| `[MMB]` | Reset pan |
| `[CLICK+DRAG]` | Pan |
| `[WHEEL]` | Zoom |
