# simd-minimizers-c

C/C++ bindings for the [simd-minimizers](https://crates.io/crates/simd-minimizers) library.

## Building

This crate is designed to be built as a static or dynamic library.

```bash
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

This will generate:
- `target/release/libsimd_minimizers_c.so` (Dynamic library)
- `target/release/libsimd_minimizers_c.a` (Static library)
- `simd_minimizers.h` (C header)

## Usage in C++

Include `simd_minimizers.h` and link against the library.

### Canonical Minimizers

Use `canonical_minimizer_positions` for just positions, or `canonical_minimizers` for full results (positions, values, super-kmers).

```cpp
#include "simd_minimizers.h"
#include <iostream>
#include <vector>
#include <cstring>

int main() {
    const char* seq = "ACGTACGTACGTACGTACGTACGT";
    size_t len = strlen(seq);
    uint8_t k = 10;
    uint8_t w = 5;

    // Create a sketcher
    SimdSketcher* sketcher = simd_sketcher_new(k, w);

    // 1. Get just positions
    MinimizerList list = canonical_minimizer_positions(sketcher, seq, len);
    
    std::cout << "Positions (" << list.len << "):" << std::endl;
    for (size_t i = 0; i < list.len; ++i) {
        std::cout << list.data[i] << " ";
    }
    std::cout << std::endl;

    free_minimizer_list(list);

    // 2. Get full details (positions, values, super-kmers)
    MinimizerResult result = canonical_minimizers(sketcher, seq, len);
    
    std::cout << "Full Result:" << std::endl;
    std::cout << "Positions: " << result.positions.len << std::endl;
    std::cout << "Values: " << result.values.len << std::endl;
    std::cout << "Super-kmers: " << result.super_kmers.len << std::endl;

    free_minimizer_result(result);

    // Free the sketcher
    simd_sketcher_free(sketcher);

    return 0;
}
```

### Forward-only Minimizers

If you only need forward-strand minimizers (not canonical), use `minimizer_positions` and `minimizers`.

```cpp
#include "simd_minimizers.h"
#include <iostream>
#include <vector>
#include <cstring>

int main() {
    const char* seq = "ACGTACGTACGTACGTACGTACGT";
    size_t len = strlen(seq);
    uint8_t k = 10;
    uint8_t w = 5;

    // Create a sketcher
    SimdSketcher* sketcher = simd_sketcher_new(k, w);

    // 1. Get forward-only positions
    MinimizerList fwd_list = minimizer_positions(sketcher, seq, len);
    
    std::cout << "Forward Positions (" << fwd_list.len << "):" << std::endl;
    for (size_t i = 0; i < fwd_list.len; ++i) {
        std::cout << fwd_list.data[i] << " ";
    }
    std::cout << std::endl;

    free_minimizer_list(fwd_list);

    // 2. Get forward-only full details
    MinimizerResult fwd_result = minimizers(sketcher, seq, len);
    
    std::cout << "Forward Full Result:" << std::endl;
    std::cout << "Positions: " << fwd_result.positions.len << std::endl;
    std::cout << "Values: " << fwd_result.values.len << std::endl;
    // Super-kmers are also available if needed

    free_minimizer_result(fwd_result);

    // Free the sketcher
    simd_sketcher_free(sketcher);

    return 0;
}
```

## API Reference

### Types

- `MinimizerList`: List of `u32` positions.
- `MinimizerValues`: List of `u64` minimizer values.
- `SuperKmerList`: List of `u64` super-kmer values.
- `MinimizerResult`: Struct containing all three lists.
- `SimdSketcher`: Opaque struct for stateful sketching.

### Functions

- `simd_sketcher_new(k, w)`: Creates a new sketcher.
- `simd_sketcher_free(sketcher)`: Frees the sketcher.
- `canonical_minimizer_positions(sketcher, seq, len)`: Computes positions using the sketcher.
- `canonical_minimizers(sketcher, seq, len)`: Computes full result using the sketcher.
- `minimizer_positions(sketcher, seq, len)`: Computes forward-only positions.
- `minimizers(sketcher, seq, len)`: Computes forward-only full result.

### Memory Management
- `free_minimizer_list(list)`: Frees `MinimizerList`.
- `free_minimizer_result(result)`: Frees `MinimizerResult`.
