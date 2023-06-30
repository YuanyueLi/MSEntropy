#!/bin/bash
docker run \
    --rm \
    -v $(pwd)/..:/src \
    -u $(id -u):$(id -g) \
    emscripten/emsdk \
    emcc -O3 -s EXPORTED_FUNCTIONS="['_malloc','_free']" -s WASM=1 -s EXPORTED_RUNTIME_METHODS="['cwrap']" \
    -o /src/javascript/msentropy.js \
    /src/javascript/MSEntropy.c /src/javascript/SpectralEntropy.c /src/javascript/CleanSpectrum.c
