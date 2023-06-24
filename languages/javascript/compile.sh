#!/bin/bash
docker run \
    --rm \
    -v $(pwd)/..:/src \
    -u $(id -u):$(id -g) \
    emscripten/emsdk \
    emcc -O3 -s EXPORTED_FUNCTIONS="['_malloc','_free']" -s WASM=1 -s EXTRA_EXPORTED_RUNTIME_METHODS="['cwrap']" \
    -o /src/language_js/msentropy.js /src/language_js/MSEntropy.c /src/language_js/SpectralEntropy.c /src/language_js/CleanSpectrum.c
