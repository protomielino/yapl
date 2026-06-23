# yapl - Yet Another Particle Life

Particle Life simulation with spatial partitioning for performance.

## Build

```
make          # build
make run      # build + run (grid backend, default)
make run-kd   # build + run with k-d tree backend
make clean
```

## Usage

```
yapl [-n <count>] [-m <colors>] [--backend grid|kd|qt] [--force-law standard|linear|lj|lennard-jones|smooth|damped-wave]
```

| Argument | Default | Description |
|---|---|---|
| `-n` | 1000 | Number of particles |
| `-m` | 6 | Number of colors / types |
| `--backend` | `grid` | Spatial partition backend (`grid`, `kd`, or `qt`) |
| `--force-law` | `standard` | Force law for pair interactions |

## Controls

| Key | Action |
|---|---|
| `[z]` | Show forces matrix (scroll to edit) |
| `[x]` | Show minDist matrix (beta per type pair) |
| `[c]` | Show masses per type |
| `[v]` | Show radii matrix (display only) |
| `[scroll]` | Hover + scroll to edit cell (±0.05, Shift=±0.25, Ctrl=±0.01) |
| `[CTRL+click]` | Select cell (yellow border) |
| `[CTRL+right-click]` | Zero active panel |
| `[SPACE]` | Randomize all |
| `[r]` | Randomize positions |
| `[p]` | Pause / unpause |
| `[d]` | Show tree partitions overlay |
| `[f]` | Cycle force law at runtime |
| `[h]` | Toggle help |
| `[H]` | Toggle HUD |
| `[R]` | Reset pan & zoom |
| `[MMB]` | Reset pan |
| `[CLICK+DRAG]` | Pan |
| `[WHEEL]` | Zoom (or edit cell when over matrix panel) |
