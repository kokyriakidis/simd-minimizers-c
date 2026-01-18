# simd-minimizers-c

C/C++ bindings for the [simd-minimizers](https://github.com/rust-seq/simd-minimizers) library - a SIMD-accelerated library to compute random minimizers and syncmers.

> **Note**: This library uses the latest `master` branch from GitHub which includes syncmers support (not yet released on crates.io).

## Requirements

- **Rust 1.87+** (for the `use<...>` precise capturing syntax)
- AVX2 (x86-64) or NEON (ARM) instruction sets

## Building

```bash
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

This will generate:
- `target/release/libsimd_minimizers_c.dylib` (macOS) or `.so` (Linux)
- `target/release/libsimd_minimizers_c.a` (Static library)
- `simd_minimizers.h` (C++ header)

## Quick Start

### Basic Usage - Positions Only

```cpp
#include "simd_minimizers.h"
#include <cstring>
#include <iostream>

int main() {
    const char* seq = "ACGTGCTCAGAGACTCAGAGGA";
    size_t len = strlen(seq);
    
    // Create sketcher with k=5, w=7
    // Note: k + w - 1 = 11 (odd) - required for canonical syncmers
    SimdSketcher* sketcher = simd_sketcher_new(5, 7);
    
    // Get minimizer positions
    MinimizerList min_pos = canonical_minimizer_positions(sketcher, seq, len);
    std::cout << "Found " << min_pos.len << " minimizers" << std::endl;
    free_minimizer_list(min_pos);
    
    // Get closed syncmer positions
    SyncmerList closed_pos = canonical_syncmer_positions(sketcher, seq, len);
    std::cout << "Found " << closed_pos.len << " closed syncmers" << std::endl;
    free_syncmer_list(closed_pos);
    
    // Get open syncmer positions (w must be odd for open syncmers)
    SyncmerList open_pos = canonical_open_syncmer_positions(sketcher, seq, len);
    std::cout << "Found " << open_pos.len << " open syncmers" << std::endl;
    free_syncmer_list(open_pos);
    
    simd_sketcher_free(sketcher);
    return 0;
}
```

### Getting Positions and Values

```cpp
#include "simd_minimizers.h"
#include <cstring>
#include <iostream>

int main() {
    const char* seq = "ACGTGCTCAGAGACTCAGAGGA";
    size_t len = strlen(seq);
    
    SimdSketcher* sketcher = simd_sketcher_new(5, 7);
    
    // Get full minimizer result (positions + values + super-kmers)
    MinimizerResult min_result = canonical_minimizers(sketcher, seq, len);
    
    std::cout << "Minimizers:" << std::endl;
    for (size_t i = 0; i < min_result.positions.len; ++i) {
        std::cout << "  Position: " << min_result.positions.data[i] 
                  << ", Value: " << min_result.values.data[i] << std::endl;
    }
    std::cout << "  Super-kmers count: " << min_result.super_kmers.len << std::endl;
    free_minimizer_result(min_result);
    
    // Get closed syncmer result (positions + values, no super-kmers)
    SyncmerResult sync_result = canonical_syncmers(sketcher, seq, len);
    
    std::cout << "Closed Syncmers:" << std::endl;
    for (size_t i = 0; i < sync_result.positions.len; ++i) {
        std::cout << "  Position: " << sync_result.positions.data[i] 
                  << ", Value: " << sync_result.values.data[i] << std::endl;
    }
    free_syncmer_result(sync_result);
    
    // Get open syncmer result
    SyncmerResult open_result = canonical_open_syncmers(sketcher, seq, len);
    
    std::cout << "Open Syncmers:" << std::endl;
    for (size_t i = 0; i < open_result.positions.len; ++i) {
        std::cout << "  Position: " << open_result.positions.data[i] 
                  << ", Value: " << open_result.values.data[i] << std::endl;
    }
    free_syncmer_result(open_result);
    
    simd_sketcher_free(sketcher);
    return 0;
}
```

## Compiling and Running Tests

```bash
# Build the library
RUSTFLAGS="-C target-cpu=native" cargo build --release

# Compile the test
g++ -std=c++17 -o test_sketcher test_sketcher.cpp -L target/release -lsimd_minimizers_c -Wl,-rpath,target/release

# Run
./test_sketcher
```

## Closed vs Open Syncmers

Syncmers are k-mers selected based on the position of the minimum s-mer (sub-kmer) within them. The library supports two types:

### Closed Syncmers

A k-mer is a **closed syncmer** if its minimum s-mer occurs at the **first or last** position within the k-mer.

```
k-mer:     [ACGTG CTCAG AGACT]  (k+w-1 = 11)
s-mers:     ^^^^^              position 0 (first)
                  ^^^^^        position w-1 (last)
```

- Selected if minimum s-mer is at position 0 **or** position `w-1`
- Higher density (more k-mers selected)
- Good for applications needing denser coverage

### Open Syncmers

A k-mer is an **open syncmer** if its minimum s-mer occurs at any position **except** the first or last (i.e., in the middle).

```
k-mer:     [ACGTG CTCAG AGACT]  (k+w-1 = 11)
s-mers:           ^^^^^        positions 1 to w-2 (middle)
```

- Selected if minimum s-mer is at positions `1` to `w-2`
- Lower density (fewer k-mers selected)
- Requires `w` to be **odd** so there's a unique middle position
- Better conservation under mutations (selected k-mers are more stable)

### Comparison

| Property | Closed Syncmers | Open Syncmers |
|----------|-----------------|---------------|
| Selection | min s-mer at first/last | min s-mer in middle |
| Density | Higher | Lower |
| Conservation | Lower | Higher |
| `w` requirement | Any | Must be odd |

## API Overview

### Minimizers

| Function | Description |
|----------|-------------|
| `canonical_minimizer_positions` | Canonical minimizer positions only |
| `canonical_minimizers` | Full result: positions, values, super-kmers |
| `minimizer_positions` | Forward-only minimizer positions |
| `minimizers` | Forward-only full result |

### Closed Syncmers

| Function | Description |
|----------|-------------|
| `canonical_syncmer_positions` | Canonical closed syncmer positions |
| `canonical_syncmers` | Positions and values (no super-kmers) |
| `syncmer_positions` | Forward-only closed syncmer positions |
| `syncmers` | Forward-only positions and values |

### Open Syncmers

| Function | Description |
|----------|-------------|
| `canonical_open_syncmer_positions` | Canonical open syncmer positions |
| `canonical_open_syncmers` | Positions and values |
| `open_syncmer_positions` | Forward-only open syncmer positions |
| `open_syncmers` | Forward-only positions and values |

### Memory Management

| Function | Description |
|----------|-------------|
| `free_minimizer_list` | Free a MinimizerList/SyncmerList |
| `free_minimizer_result` | Free a MinimizerResult/SyncmerResult |
| `free_syncmer_list` | Alias for free_minimizer_list |
| `free_syncmer_result` | Alias for free_minimizer_result |

## Data Structures

### MinimizerList / SyncmerList

```cpp
struct MinimizerList {
    uint32_t* data;    // Array of positions
    size_t len;        // Number of elements
    size_t capacity;   // Internal capacity
};
```

### MinimizerResult / SyncmerResult

```cpp
struct MinimizerResult {
    MinimizerList positions;   // Position indices
    MinimizerValues values;    // Hash values (uint64_t*)
    SuperKmerList super_kmers; // Super-kmer indices (minimizers only)
};
```

## Parameters

- **k**: k-mer length (the s-mer length for syncmers, typically 5-31)
- **w**: Window size for minimizers, or determines the syncmer length as `k + w - 1`

### Important Constraints

| Constraint | Reason |
|------------|--------|
| Canonical syncmers: `k + w - 1` must be **odd** | Guarantees canonicality |
| Open syncmers: `w` must be **odd** | Ensures unique middle position |
| Super-kmers: only for minimizers | Not supported for syncmers |

## Performance

The library uses SIMD (AVX2/NEON) for high performance:
- Processes sequences in 8 parallel chunks
- Uses rolling hash (ntHash) for efficient k-mer hashing
- Sliding window minimum via "two stacks" algorithm

Always compile with `-C target-cpu=native` for best results.

## License

MIT - see the [simd-minimizers](https://github.com/rust-seq/simd-minimizers) repository for details.
