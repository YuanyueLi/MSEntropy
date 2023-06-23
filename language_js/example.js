var Module = require('./msentropy.js');

arrayMalloc = (arr) => {
    // Get data byte size, allocate memory on Emscripten heap, and get pointer
    var nDataBytes = arr.length * arr.BYTES_PER_ELEMENT;
    var dataPtr = Module._malloc(nDataBytes);

    // Copy data to Emscripten heap (directly accessed from Module.HEAPU8)
    var dataHeap = new Uint8Array(
        Module.HEAPU8.buffer,
        dataPtr,
        nDataBytes
    );
    dataHeap.set(new Uint8Array(arr.buffer));

    return dataHeap;
};
arrayFreeAndExtract = (dataHeap, arrayType, bytes) => {
    // Free memory
    let result = new arrayType(
        dataHeap.buffer,
        dataHeap.byteOffset,
        dataHeap.length / bytes
    );
    Module._free(dataHeap.byteOffset);
    return result;
};
arrayFree = (dataHeap) => Module._free(dataHeap.byteOffset);
Module.onRuntimeInitialized = async () => {
    const api = {
        clean_spectrum: Module.cwrap("wasm_clean_spectrum", "number",
            ["number", "number", "number", "number", "number", "number", "number", "number", "number",]),
        calculate_entropy_similarity: Module.cwrap("wasm_calculate_entropy_similarity", "number",
            ["number", "number", "number", "number", "number", "number", "number", "number", "number", "number", "number", "number",]),
        calculate_unweighted_entropy_similarity: Module.cwrap("wasm_calculate_unweighted_entropy_similarity", "number",
            ["number", "number", "number", "number", "number", "number", "number", "number", "number", "number", "number", "number",]),
    };

    let arraryRaw = new Float32Array([
        41.04, 0.3716, 0, 0.3716, 69.07, 7.917962, 69.07, -7.917962, 69.071,
        100, 86.0969, 66.83, 86.01, 10,
    ]);

    let arrayHeap = arrayMalloc(arraryRaw);
    let arr_len = api.clean_spectrum(arrayHeap.byteOffset, arraryRaw.length / 2, 0, -1, 0.01, 0.05, -1, 5, 1);
    let result = arrayFreeAndExtract(arrayHeap, Float32Array, 4).slice(0, arr_len * 2);
    console.log(result);


    let arrayHeapA = arrayMalloc(new Float32Array([69.071, 7.917962, 86.066, 1.021589, 86.0969, 100.0]));
    let arrayHeapB = arrayMalloc(new Float32Array([41.04, 37.16, 69.07, 66.83, 86.1, 999.0]));
    let similarity = api.calculate_entropy_similarity(
        arrayHeapA.byteOffset, 3,
        arrayHeapB.byteOffset, 3,
        0.02, -1, 1, -1, -1, 0.01, -1);
    arrayFree(arrayHeapA);
    arrayFree(arrayHeapB);
    console.log("Entropy similarity:", similarity);

    arrayHeapA = arrayMalloc(new Float32Array([69.071, 7.917962, 86.066, 1.021589, 86.0969, 100.0]));
    arrayHeapB = arrayMalloc(new Float32Array([41.04, 37.16, 69.07, 66.83, 86.1, 999.0]));
    similarity = api.calculate_unweighted_entropy_similarity(
        arrayHeapA.byteOffset, 3,
        arrayHeapB.byteOffset, 3,
        0.02, -1, 1, -1, -1, 0.01, -1);
    arrayFree(arrayHeapA);
    arrayFree(arrayHeapB);
    console.log("Unweighted entropy similarity:", similarity);
}
