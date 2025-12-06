#include <cstdarg>
#include <cstdint>
#include <cstdlib>
#include <ostream>
#include <new>

/// Wrapper struct for the Sketcher to be used via C pointer.
struct SimdSketcher;

/// A list of minimizer positions (indices).
struct MinimizerList {
  /// Pointer to the array of positions.
  uint32_t *data;
  /// Number of elements in the array.
  size_t len;
  /// Capacity of the allocated array (internal use).
  size_t capacity;
};

/// A list of minimizer values (hashes).
struct MinimizerValues {
  /// Pointer to the array of values.
  uint64_t *data;
  /// Number of elements in the array.
  size_t len;
  /// Capacity of the allocated array (internal use).
  size_t capacity;
};

/// A list of super-kmer positions.
struct SuperKmerList {
  /// Pointer to the array of super-kmer positions.
  uint32_t *data;
  /// Number of elements in the array.
  size_t len;
  /// Capacity of the allocated array (internal use).
  size_t capacity;
};

/// Result structure containing positions, values, and super-kmers.
struct MinimizerResult {
  MinimizerList positions;
  MinimizerValues values;
  SuperKmerList super_kmers;
};

extern "C" {

/// Creates a new SimdSketcher.
SimdSketcher *simd_sketcher_new(uint8_t k, uint8_t w);

/// Frees the SimdSketcher.
void simd_sketcher_free(SimdSketcher *sketcher);

/// Computes minimizer positions using the sketcher.
MinimizerList canonical_minimizer_positions(SimdSketcher *sketcher, const char *seq, size_t len);

/// Computes minimizers (positions, values, super-kmers) using the sketcher.
MinimizerResult canonical_minimizers(SimdSketcher *sketcher, const char *seq, size_t len);

/// Computes forward-only minimizer positions using the sketcher.
MinimizerList minimizer_positions(SimdSketcher *sketcher, const char *seq, size_t len);

/// Computes forward-only minimizers (positions, values, super-kmers) using the sketcher.
MinimizerResult minimizers(SimdSketcher *sketcher, const char *seq, size_t len);

/// Frees the memory associated with a MinimizerList.
void free_minimizer_list(MinimizerList list);

/// Frees the memory associated with a MinimizerResult.
void free_minimizer_result(MinimizerResult result);

}  // extern "C"
