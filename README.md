# Acrobot
Simulation of the *Acrobot* with plain-vanilla `C++` running in the browser.

Run it directly online <https://raw.githack.com/olofer/acrobot/main/payload/index.html> or clone and build/run locally (see below).

## Usage

- `R` reset to random state
- `Z` reset to (arbitrary) "zero" (potential energy) state
- `U` reset to upright position (maximum pot. energy)
- `D` reset to hanging position (minimum pot. energy)
- `up/down` apply torque at common joint in either direection
- `left/right` increase/decrease friction at common joint

More options: TBD

## Build WebAssembly & Run locally
Requires (emscripten) `emcc` compiler and `emrun` utility. Scipts may or may not assume a WSL2 environment.

```
./build.sh
```

This should recreate the `payload` folder.

```
./run-html.sh
```

This should start a local browser and load the *Acrobot* application.

## Octave-based test
From the `test-octave` folder run:

```
octave --no-window-system --eval "acrobot_rebuild_test"
```

Several other `octave` programs and scripts are available. These may be useful as documentation of the governing equations.
