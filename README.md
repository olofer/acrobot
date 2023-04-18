# Acrobot
Simulation of the *Acrobot* with plain-vanilla `C++` running in the browser.

Run it directly online <https://raw.githack.com/olofer/acrobot/main/payload/index.html> or clone and build/run locally (see below).

## Usage
- `R` reset to random state
- `Z` reset to (arbitrary) "zero" (potential energy) state
- `U` reset to upright position (maximum pot. energy)
- `D` reset to hanging position (minimum pot. energy)
- `up/down` apply torque at common joint (or change lock/hold angle)
- `left/right` increase/decrease friction at common joint
- `F` toggle lock of the arm relative angle
- `H` toggle hold of the second arm absolute angle
- `B` apply brake friction at main joint
- `P` reset and enter pump-mode (or exit pump-mode into lock mode) 

## Build WebAssembly & run locally
Requires (emscripten) `emcc` compiler and `emrun` utility. Scripts may or may not assume a WSL2 environment.

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

### References
Tutorial on combining C++/WASM found here: 
- https://github.com/olafurw/talk-accu-webassembly
