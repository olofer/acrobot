#!/bin/bash
rm -rf payload
rm -f browser_sim.wasm;
emcc -std=c++17 -Wall -O2 --no-entry -s STANDALONE_WASM browser_sim.cpp -o browser_sim.wasm;

# The complete browser-runnable program files are put in folder "payload"
mkdir payload
mv browser_sim.wasm payload/.
cp acrobot.js payload/.
cp index.html payload/.
