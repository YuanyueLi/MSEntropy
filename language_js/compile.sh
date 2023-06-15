#!/bin/bash
docker run \
    --rm \
    -v $(pwd):/src \
    -u $(id -u):$(id -g) \
    emscripten/emsdk \
    emcc -O3 -s EXPORTED_FUNCTIONS="['_malloc','_free']" -s WASM=1 -s EXTRA_EXPORTED_RUNTIME_METHODS="['cwrap']" \
    -o /src/entropy.js /src/SpectralEntropy.c
