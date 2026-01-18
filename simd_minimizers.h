#include <cstdarg>
#include <cstdint>
#include <cstdlib>
#include <ostream>
#include <new>

/// Opaque handle to the SIMD sketcher for computing minimizers and syncmers.
///
/// The sketcher stores the k-mer size (k), window size (w), and an internal cache
/// for SIMD computations. Create with `simd_sketcher_new` and free with `simd_sketcher_free`.
struct SimdSketcher;

/// A list of minimizer/syncmer positions (0-based indices into the sequence).
struct MinimizerList {
  /// Pointer to the array of positions.
  uint32_t *data;
  /// Number of elements in the array.
  size_t len;
  /// Capacity of the allocated array (for internal memory management).
  size_t capacity;
};

/// A list of minimizer/syncmer values (2-bit packed k-mer encodings as u64).
///
/// Each value represents the canonical or forward k-mer at the corresponding position,
/// encoded with 2 bits per nucleotide (A=0, C=1, G=2, T=3).
struct MinimizerValues {
  /// Pointer to the array of k-mer values.
  uint64_t *data;
  /// Number of elements in the array.
  size_t len;
  /// Capacity of the allocated array (for internal memory management).
  size_t capacity;
};

/// A list of super-kmer indices.
///
/// Each index corresponds to a minimizer position and indicates which super-kmer
/// that minimizer belongs to. Consecutive minimizers with the same super-kmer index
/// are part of the same super-kmer (a maximal sequence where all windows share
/// the same minimizer).
///
/// Note: Super-kmers are only available for minimizers, not for syncmers.
struct SuperKmerList {
  /// Pointer to the array of super-kmer indices.
  uint32_t *data;
  /// Number of elements in the array.
  size_t len;
  /// Capacity of the allocated array (for internal memory management).
  size_t capacity;
};

/// Result structure containing positions, values, and super-kmer indices.
struct MinimizerResult {
  /// The positions (0-based indices) of each minimizer/syncmer.
  MinimizerList positions;
  /// The k-mer values (2-bit packed encodings) at each position.
  MinimizerValues values;
  /// Super-kmer indices for each minimizer (empty for syncmers).
  SuperKmerList super_kmers;
};

/// Type alias for syncmer position lists (same structure as MinimizerList).
using SyncmerList = MinimizerList;

/// Type alias for syncmer results (same structure as MinimizerResult).
/// Note: The super_kmers field will always be empty for syncmers.
using SyncmerResult = MinimizerResult;

extern "C" {

/// Creates a new SimdSketcher with the given parameters.
///
/// # Parameters
/// - `k`: k-mer size (length of the k-mers to extract, typically 15-31).
/// - `w`: window size for minimizers (number of consecutive k-mers per window),
///        or s-mer size for syncmers (size of sub-k-mer used for selection).
///
/// # Returns
/// A pointer to the new SimdSketcher, or null on allocation failure.
/// Must be freed with `simd_sketcher_free` when no longer needed.
SimdSketcher *simd_sketcher_new(uint8_t k, uint8_t w);

/// Frees the SimdSketcher and releases all associated memory.
///
/// # Safety
/// - `sketcher` must be a valid pointer returned by `simd_sketcher_new`, or null.
/// - The pointer must not be used after this call.
/// - Must not be called twice on the same pointer.
void simd_sketcher_free(SimdSketcher *sketcher);

/// Computes canonical minimizer positions using SIMD acceleration.
///
/// Canonical minimizers consider both the forward k-mer and its reverse complement,
/// selecting the lexicographically smaller one. This ensures that a sequence and
/// its reverse complement produce the same minimizers.
///
/// # Parameters
/// - `sketcher`: A valid SimdSketcher pointer (created with `simd_sketcher_new`).
/// - `seq`: Pointer to the ASCII DNA sequence (A, C, G, T characters).
/// - `len`: Length of the sequence in bytes.
///
/// # Returns
/// A MinimizerList containing 0-based positions of all minimizers.
/// Must be freed with `free_minimizer_list` when no longer needed.
///
/// # Safety
/// - `sketcher` must be a valid pointer or null (returns empty list if null).
/// - `seq` must point to a valid memory region of at least `len` bytes.
MinimizerList canonical_minimizer_positions(SimdSketcher *sketcher, const char *seq, size_t len);

/// Computes canonical minimizers with positions, values, and super-kmer indices.
///
/// Canonical minimizers consider both the forward k-mer and its reverse complement,
/// selecting the lexicographically smaller one. This ensures that a sequence and
/// its reverse complement produce the same minimizers.
///
/// # Parameters
/// - `sketcher`: A valid SimdSketcher pointer (created with `simd_sketcher_new`).
/// - `seq`: Pointer to the ASCII DNA sequence (A, C, G, T characters).
/// - `len`: Length of the sequence in bytes.
///
/// # Returns
/// A MinimizerResult containing:
/// - `positions`: 0-based indices of minimizers in the sequence.
/// - `values`: 2-bit packed canonical k-mer encodings at each position.
/// - `super_kmers`: Indices grouping minimizers by their super-kmer.
///
/// Must be freed with `free_minimizer_result` when no longer needed.
///
/// # Safety
/// - `sketcher` must be a valid pointer or null (returns empty result if null).
/// - `seq` must point to a valid memory region of at least `len` bytes.
MinimizerResult canonical_minimizers(SimdSketcher *sketcher, const char *seq, size_t len);

/// Computes forward-only minimizer positions using SIMD acceleration.
///
/// Forward-only minimizers only consider the forward strand k-mers, without
/// comparing to reverse complements. Use this when strand orientation matters
/// or for slightly faster computation when canonicality is not required.
///
/// # Parameters
/// - `sketcher`: A valid SimdSketcher pointer (created with `simd_sketcher_new`).
/// - `seq`: Pointer to the ASCII DNA sequence (A, C, G, T characters).
/// - `len`: Length of the sequence in bytes.
///
/// # Returns
/// A MinimizerList containing 0-based positions of all minimizers.
/// Must be freed with `free_minimizer_list` when no longer needed.
///
/// # Safety
/// - `sketcher` must be a valid pointer or null (returns empty list if null).
/// - `seq` must point to a valid memory region of at least `len` bytes.
MinimizerList minimizer_positions(SimdSketcher *sketcher, const char *seq, size_t len);

/// Computes forward-only minimizers with positions, values, and super-kmer indices.
///
/// Forward-only minimizers only consider the forward strand k-mers, without
/// comparing to reverse complements. Use this when strand orientation matters
/// or for slightly faster computation when canonicality is not required.
///
/// # Parameters
/// - `sketcher`: A valid SimdSketcher pointer (created with `simd_sketcher_new`).
/// - `seq`: Pointer to the ASCII DNA sequence (A, C, G, T characters).
/// - `len`: Length of the sequence in bytes.
///
/// # Returns
/// A MinimizerResult containing:
/// - `positions`: 0-based indices of minimizers in the sequence.
/// - `values`: 2-bit packed forward k-mer encodings at each position.
/// - `super_kmers`: Indices grouping minimizers by their super-kmer.
///
/// Must be freed with `free_minimizer_result` when no longer needed.
///
/// # Safety
/// - `sketcher` must be a valid pointer or null (returns empty result if null).
/// - `seq` must point to a valid memory region of at least `len` bytes.
MinimizerResult minimizers(SimdSketcher *sketcher, const char *seq, size_t len);

/// Computes canonical closed syncmer positions using SIMD acceleration.
///
/// Closed syncmers select k-mers where the minimum s-mer appears at the start
/// or end of the k-mer. Canonical mode considers both strands and selects the
/// lexicographically smaller orientation.
///
/// # Parameters
/// - `sketcher`: A valid SimdSketcher pointer. Here, `k` is the k-mer size and
///   `w` is the s-mer size (size of sub-k-mers used for selection).
/// - `seq`: Pointer to the ASCII DNA sequence (A, C, G, T characters).
/// - `len`: Length of the sequence in bytes.
///
/// # Returns
/// A SyncmerList containing 0-based positions of all closed syncmers.
/// Must be freed with `free_syncmer_list` when no longer needed.
///
/// # Constraints
/// For canonical closed syncmers, `k + w - 1` must be odd.
///
/// # Safety
/// - `sketcher` must be a valid pointer or null (returns empty list if null).
/// - `seq` must point to a valid memory region of at least `len` bytes.
SyncmerList canonical_syncmer_positions(SimdSketcher *sketcher, const char *seq, size_t len);

/// Computes canonical closed syncmers with positions and values.
///
/// Closed syncmers select k-mers where the minimum s-mer appears at the start
/// or end of the k-mer. Canonical mode considers both strands and selects the
/// lexicographically smaller orientation.
///
/// # Parameters
/// - `sketcher`: A valid SimdSketcher pointer. Here, `k` is the k-mer size and
///   `w` is the s-mer size (size of sub-k-mers used for selection).
/// - `seq`: Pointer to the ASCII DNA sequence (A, C, G, T characters).
/// - `len`: Length of the sequence in bytes.
///
/// # Returns
/// A SyncmerResult containing:
/// - `positions`: 0-based indices of syncmers in the sequence.
/// - `values`: 2-bit packed canonical k-mer encodings at each position.
/// - `super_kmers`: Always empty (super-kmers are not supported for syncmers).
///
/// Must be freed with `free_syncmer_result` when no longer needed.
///
/// # Constraints
/// For canonical closed syncmers, `k + w - 1` must be odd.
///
/// # Safety
/// - `sketcher` must be a valid pointer or null (returns empty result if null).
/// - `seq` must point to a valid memory region of at least `len` bytes.
SyncmerResult canonical_syncmers(SimdSketcher *sketcher, const char *seq, size_t len);

/// Computes forward-only closed syncmer positions using SIMD acceleration.
///
/// Closed syncmers select k-mers where the minimum s-mer appears at the start
/// or end of the k-mer. Forward-only mode considers only the forward strand.
///
/// # Parameters
/// - `sketcher`: A valid SimdSketcher pointer. Here, `k` is the k-mer size and
///   `w` is the s-mer size (size of sub-k-mers used for selection).
/// - `seq`: Pointer to the ASCII DNA sequence (A, C, G, T characters).
/// - `len`: Length of the sequence in bytes.
///
/// # Returns
/// A SyncmerList containing 0-based positions of all closed syncmers.
/// Must be freed with `free_syncmer_list` when no longer needed.
///
/// # Safety
/// - `sketcher` must be a valid pointer or null (returns empty list if null).
/// - `seq` must point to a valid memory region of at least `len` bytes.
SyncmerList syncmer_positions(SimdSketcher *sketcher, const char *seq, size_t len);

/// Computes forward-only closed syncmers with positions and values.
///
/// Closed syncmers select k-mers where the minimum s-mer appears at the start
/// or end of the k-mer. Forward-only mode considers only the forward strand.
///
/// # Parameters
/// - `sketcher`: A valid SimdSketcher pointer. Here, `k` is the k-mer size and
///   `w` is the s-mer size (size of sub-k-mers used for selection).
/// - `seq`: Pointer to the ASCII DNA sequence (A, C, G, T characters).
/// - `len`: Length of the sequence in bytes.
///
/// # Returns
/// A SyncmerResult containing:
/// - `positions`: 0-based indices of syncmers in the sequence.
/// - `values`: 2-bit packed forward k-mer encodings at each position.
/// - `super_kmers`: Always empty (super-kmers are not supported for syncmers).
///
/// Must be freed with `free_syncmer_result` when no longer needed.
///
/// # Safety
/// - `sketcher` must be a valid pointer or null (returns empty result if null).
/// - `seq` must point to a valid memory region of at least `len` bytes.
SyncmerResult syncmers(SimdSketcher *sketcher, const char *seq, size_t len);

/// Computes canonical open syncmer positions using SIMD acceleration.
///
/// Open syncmers select k-mers where the minimum s-mer appears at a specific
/// offset (typically the middle position). Canonical mode considers both strands
/// and selects the lexicographically smaller orientation.
///
/// # Parameters
/// - `sketcher`: A valid SimdSketcher pointer. Here, `k` is the k-mer size and
///   `w` is the s-mer size (size of sub-k-mers used for selection).
/// - `seq`: Pointer to the ASCII DNA sequence (A, C, G, T characters).
/// - `len`: Length of the sequence in bytes.
///
/// # Returns
/// A SyncmerList containing 0-based positions of all open syncmers.
/// Must be freed with `free_syncmer_list` when no longer needed.
///
/// # Constraints
/// For open syncmers, `w` (the s-mer size) must be odd.
///
/// # Safety
/// - `sketcher` must be a valid pointer or null (returns empty list if null).
/// - `seq` must point to a valid memory region of at least `len` bytes.
SyncmerList canonical_open_syncmer_positions(SimdSketcher *sketcher, const char *seq, size_t len);

/// Computes canonical open syncmers with positions and values.
///
/// Open syncmers select k-mers where the minimum s-mer appears at a specific
/// offset (typically the middle position). Canonical mode considers both strands
/// and selects the lexicographically smaller orientation.
///
/// # Parameters
/// - `sketcher`: A valid SimdSketcher pointer. Here, `k` is the k-mer size and
///   `w` is the s-mer size (size of sub-k-mers used for selection).
/// - `seq`: Pointer to the ASCII DNA sequence (A, C, G, T characters).
/// - `len`: Length of the sequence in bytes.
///
/// # Returns
/// A SyncmerResult containing:
/// - `positions`: 0-based indices of syncmers in the sequence.
/// - `values`: 2-bit packed canonical k-mer encodings at each position.
/// - `super_kmers`: Always empty (super-kmers are not supported for syncmers).
///
/// Must be freed with `free_syncmer_result` when no longer needed.
///
/// # Constraints
/// For open syncmers, `w` (the s-mer size) must be odd.
///
/// # Safety
/// - `sketcher` must be a valid pointer or null (returns empty result if null).
/// - `seq` must point to a valid memory region of at least `len` bytes.
SyncmerResult canonical_open_syncmers(SimdSketcher *sketcher, const char *seq, size_t len);

/// Computes forward-only open syncmer positions using SIMD acceleration.
///
/// Open syncmers select k-mers where the minimum s-mer appears at a specific
/// offset (typically the middle position). Forward-only mode considers only
/// the forward strand.
///
/// # Parameters
/// - `sketcher`: A valid SimdSketcher pointer. Here, `k` is the k-mer size and
///   `w` is the s-mer size (size of sub-k-mers used for selection).
/// - `seq`: Pointer to the ASCII DNA sequence (A, C, G, T characters).
/// - `len`: Length of the sequence in bytes.
///
/// # Returns
/// A SyncmerList containing 0-based positions of all open syncmers.
/// Must be freed with `free_syncmer_list` when no longer needed.
///
/// # Constraints
/// For open syncmers, `w` (the s-mer size) must be odd.
///
/// # Safety
/// - `sketcher` must be a valid pointer or null (returns empty list if null).
/// - `seq` must point to a valid memory region of at least `len` bytes.
SyncmerList open_syncmer_positions(SimdSketcher *sketcher, const char *seq, size_t len);

/// Computes forward-only open syncmers with positions and values.
///
/// Open syncmers select k-mers where the minimum s-mer appears at a specific
/// offset (typically the middle position). Forward-only mode considers only
/// the forward strand.
///
/// # Parameters
/// - `sketcher`: A valid SimdSketcher pointer. Here, `k` is the k-mer size and
///   `w` is the s-mer size (size of sub-k-mers used for selection).
/// - `seq`: Pointer to the ASCII DNA sequence (A, C, G, T characters).
/// - `len`: Length of the sequence in bytes.
///
/// # Returns
/// A SyncmerResult containing:
/// - `positions`: 0-based indices of syncmers in the sequence.
/// - `values`: 2-bit packed forward k-mer encodings at each position.
/// - `super_kmers`: Always empty (super-kmers are not supported for syncmers).
///
/// Must be freed with `free_syncmer_result` when no longer needed.
///
/// # Constraints
/// For open syncmers, `w` (the s-mer size) must be odd.
///
/// # Safety
/// - `sketcher` must be a valid pointer or null (returns empty result if null).
/// - `seq` must point to a valid memory region of at least `len` bytes.
SyncmerResult open_syncmers(SimdSketcher *sketcher, const char *seq, size_t len);

/// Frees the memory associated with a MinimizerList.
///
/// # Safety
/// - `list` must be a valid MinimizerList returned by a minimizer positions function.
/// - Must not be called twice on the same list.
/// - After calling, the list's data pointer is invalid and must not be accessed.
void free_minimizer_list(MinimizerList list);

/// Frees the memory associated with a MinimizerResult.
///
/// Frees all three arrays (positions, values, super_kmers) contained in the result.
///
/// # Safety
/// - `result` must be a valid MinimizerResult returned by a minimizer function.
/// - Must not be called twice on the same result.
/// - After calling, all data pointers in the result are invalid and must not be accessed.
void free_minimizer_result(MinimizerResult result);

/// Frees the memory associated with a SyncmerList.
///
/// This is an alias for `free_minimizer_list` since SyncmerList and MinimizerList
/// have the same structure.
///
/// # Safety
/// - `list` must be a valid SyncmerList returned by a syncmer positions function.
/// - Must not be called twice on the same list.
/// - After calling, the list's data pointer is invalid and must not be accessed.
void free_syncmer_list(SyncmerList list);

/// Frees the memory associated with a SyncmerResult.
///
/// This is an alias for `free_minimizer_result` since SyncmerResult and MinimizerResult
/// have the same structure.
///
/// # Safety
/// - `result` must be a valid SyncmerResult returned by a syncmer function.
/// - Must not be called twice on the same result.
/// - After calling, all data pointers in the result are invalid and must not be accessed.
void free_syncmer_result(SyncmerResult result);

}  // extern "C"
