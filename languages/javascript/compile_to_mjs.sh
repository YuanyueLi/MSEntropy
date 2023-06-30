#!/bin/bash
docker run \
    --rm \
    -v $(pwd)/..:/src \
    -u $(id -u):$(id -g) \
    emscripten/emsdk \
    emcc --no-entry \
    -s ENVIRONMENT='web' \
    -s SINGLE_FILE=1 \
    -s EXPORT_NAME='createModule' \
    -s USE_ES6_IMPORT_META=0 \
    -s EXPORTED_FUNCTIONS='["_malloc", "_free"]' \
    -s EXPORTED_RUNTIME_METHODS='["ccall", "cwrap"]' \
    -O3 \
    -o /src/javascript/msentropy.mjs \
    /src/javascript/MSEntropy.c /src/javascript/SpectralEntropy.c /src/javascript/CleanSpectrum.c
