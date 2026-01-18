//! C/C++ FFI bindings for the simd-minimizers library.
//!
//! This module provides C-compatible functions for computing minimizers and syncmers
//! on DNA sequences using SIMD acceleration.

use libc::{c_char, size_t};
use simd_minimizers::seq_hash::NtHasher;
use simd_minimizers::Cache;
use simd_minimizers::collect::CollectAndDedup;
use simd_minimizers::private::minimizers::canonical_minimizers_seq_simd;
use packed_seq::{PackedSeqVec, SeqVec, Seq};
use std::slice;

/// A list of minimizer/syncmer positions (0-based indices into the sequence).
#[repr(C)]
pub struct MinimizerList {
    /// Pointer to the array of positions.
    pub data: *mut u32,
    /// Number of elements in the array.
    pub len: size_t,
    /// Capacity of the allocated array (for internal memory management).
    pub capacity: size_t,
}

/// A list of minimizer/syncmer values (2-bit packed k-mer encodings as u64).
///
/// Each value represents the canonical or forward k-mer at the corresponding position,
/// encoded with 2 bits per nucleotide (A=0, C=1, G=2, T=3).
#[repr(C)]
pub struct MinimizerValues {
    /// Pointer to the array of k-mer values.
    pub data: *mut u64,
    /// Number of elements in the array.
    pub len: size_t,
    /// Capacity of the allocated array (for internal memory management).
    pub capacity: size_t,
}

/// A list of super-kmer indices.
///
/// Each index corresponds to a minimizer position and indicates which super-kmer
/// that minimizer belongs to. Consecutive minimizers with the same super-kmer index
/// are part of the same super-kmer (a maximal sequence where all windows share
/// the same minimizer).
///
/// Note: Super-kmers are only available for minimizers, not for syncmers.
#[repr(C)]
pub struct SuperKmerList {
    /// Pointer to the array of super-kmer indices.
    pub data: *mut u32,
    /// Number of elements in the array.
    pub len: size_t,
    /// Capacity of the allocated array (for internal memory management).
    pub capacity: size_t,
}

/// Result structure containing positions, values, and super-kmer indices.
#[repr(C)]
pub struct MinimizerResult {
    /// The positions (0-based indices) of each minimizer/syncmer.
    pub positions: MinimizerList,
    /// The k-mer values (2-bit packed encodings) at each position.
    pub values: MinimizerValues,
    /// Super-kmer indices for each minimizer (empty for syncmers).
    pub super_kmers: SuperKmerList,
}

/// Type alias for syncmer position lists (same structure as MinimizerList).
pub type SyncmerList = MinimizerList;

/// Type alias for syncmer value lists (same structure as MinimizerValues).
pub type SyncmerValues = MinimizerValues;

/// Type alias for syncmer results (same structure as MinimizerResult).
/// Note: The super_kmers field will always be empty for syncmers.
pub type SyncmerResult = MinimizerResult;

/// Opaque handle to the SIMD sketcher for computing minimizers and syncmers.
///
/// The sketcher stores the k-mer size (k), window size (w), and an internal cache
/// for SIMD computations. Create with `simd_sketcher_new` and free with `simd_sketcher_free`.
pub struct SimdSketcher {
    /// k-mer size: the length of the k-mers to extract.
    k: usize,
    /// Window size: the number of consecutive k-mers to consider for minimizer selection.
    /// For syncmers, this is the s-mer size (size of the sub-k-mer used for selection).
    w: usize,
    /// Internal cache for SIMD operations.
    cache: Cache,
}

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
#[unsafe(no_mangle)]
pub extern "C" fn simd_sketcher_new(k: u8, w: u8) -> *mut SimdSketcher {
    let sketcher = SimdSketcher {
        k: k as usize,
        w: w as usize,
        cache: Cache::default(),
    };
    Box::into_raw(Box::new(sketcher))
}

/// Frees the SimdSketcher and releases all associated memory.
///
/// # Safety
/// - `sketcher` must be a valid pointer returned by `simd_sketcher_new`, or null.
/// - The pointer must not be used after this call.
/// - Must not be called twice on the same pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn simd_sketcher_free(sketcher: *mut SimdSketcher) {
    if !sketcher.is_null() {
        unsafe {
            let _ = Box::from_raw(sketcher);
        }
    }
}

// =============================================================================
// MINIMIZERS
// =============================================================================

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
#[unsafe(no_mangle)]
pub unsafe extern "C" fn canonical_minimizer_positions(
    sketcher: *mut SimdSketcher,
    seq: *const c_char,
    len: size_t,
) -> MinimizerList {
    if sketcher.is_null() {
        return MinimizerList { data: std::ptr::null_mut(), len: 0, capacity: 0 };
    }

    let sketcher = unsafe { &mut *sketcher };
    let seq_slice = unsafe { slice::from_raw_parts(seq as *const u8, len) };
    let packed_seq = PackedSeqVec::from_ascii(seq_slice);

    let hasher = NtHasher::<true>::new(sketcher.k);
    let mut positions = Vec::new();

    canonical_minimizers_seq_simd(packed_seq.as_slice(), &hasher, sketcher.w, &mut sketcher.cache)
        .collect_and_dedup_into::<false>(&mut positions);

    let ptr = positions.as_mut_ptr();
    let len = positions.len();
    let cap = positions.capacity();
    std::mem::forget(positions);

    MinimizerList { data: ptr, len, capacity: cap }
}

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
#[unsafe(no_mangle)]
pub unsafe extern "C" fn canonical_minimizers(
    sketcher: *mut SimdSketcher,
    seq: *const c_char,
    len: size_t,
) -> MinimizerResult {
    if sketcher.is_null() {
        return MinimizerResult {
            positions: MinimizerList { data: std::ptr::null_mut(), len: 0, capacity: 0 },
            values: MinimizerValues { data: std::ptr::null_mut(), len: 0, capacity: 0 },
            super_kmers: SuperKmerList { data: std::ptr::null_mut(), len: 0, capacity: 0 },
        };
    }

    let sketcher = unsafe { &mut *sketcher };
    let seq_slice = unsafe { slice::from_raw_parts(seq as *const u8, len) };
    let packed_seq = PackedSeqVec::from_ascii(seq_slice);

    let hasher = NtHasher::<true>::new(sketcher.k);
    let mut positions = Vec::new();
    let mut super_kmers = Vec::new();

    canonical_minimizers_seq_simd(packed_seq.as_slice(), &hasher, sketcher.w, &mut sketcher.cache)
        .collect_and_dedup_with_index_into(&mut positions, &mut super_kmers);

    let mut values = Vec::with_capacity(positions.len());
    for &pos in &positions {
        let val = {
            let a = packed_seq.as_slice().read_kmer(sketcher.k, pos as usize);
            let b = packed_seq.as_slice().read_revcomp_kmer(sketcher.k, pos as usize);
            core::cmp::min(a, b)
        };
        values.push(val);
    }

    let pos_ptr = positions.as_mut_ptr();
    let pos_len = positions.len();
    let pos_cap = positions.capacity();
    std::mem::forget(positions);

    let val_ptr = values.as_mut_ptr();
    let val_len = values.len();
    let val_cap = values.capacity();
    std::mem::forget(values);

    let sk_ptr = super_kmers.as_mut_ptr();
    let sk_len = super_kmers.len();
    let sk_cap = super_kmers.capacity();
    std::mem::forget(super_kmers);

    MinimizerResult {
        positions: MinimizerList { data: pos_ptr, len: pos_len, capacity: pos_cap },
        values: MinimizerValues { data: val_ptr, len: val_len, capacity: val_cap },
        super_kmers: SuperKmerList { data: sk_ptr, len: sk_len, capacity: sk_cap },
    }
}

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
#[unsafe(no_mangle)]
pub unsafe extern "C" fn minimizer_positions(
    sketcher: *mut SimdSketcher,
    seq: *const c_char,
    len: size_t,
) -> MinimizerList {
    if sketcher.is_null() {
        return MinimizerList { data: std::ptr::null_mut(), len: 0, capacity: 0 };
    }

    let sketcher = unsafe { &mut *sketcher };
    let seq_slice = unsafe { slice::from_raw_parts(seq as *const u8, len) };
    let packed_seq = PackedSeqVec::from_ascii(seq_slice);

    let hasher = NtHasher::<false>::new(sketcher.k);
    let mut positions = Vec::new();

    use simd_minimizers::private::minimizers::minimizers_seq_simd;

    minimizers_seq_simd(packed_seq.as_slice(), &hasher, sketcher.w, &mut sketcher.cache)
        .collect_and_dedup_into::<false>(&mut positions);

    let ptr = positions.as_mut_ptr();
    let len = positions.len();
    let cap = positions.capacity();
    std::mem::forget(positions);

    MinimizerList { data: ptr, len, capacity: cap }
}

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
#[unsafe(no_mangle)]
pub unsafe extern "C" fn minimizers(
    sketcher: *mut SimdSketcher,
    seq: *const c_char,
    len: size_t,
) -> MinimizerResult {
    if sketcher.is_null() {
        return MinimizerResult {
            positions: MinimizerList { data: std::ptr::null_mut(), len: 0, capacity: 0 },
            values: MinimizerValues { data: std::ptr::null_mut(), len: 0, capacity: 0 },
            super_kmers: SuperKmerList { data: std::ptr::null_mut(), len: 0, capacity: 0 },
        };
    }

    let sketcher = unsafe { &mut *sketcher };
    let seq_slice = unsafe { slice::from_raw_parts(seq as *const u8, len) };
    let packed_seq = PackedSeqVec::from_ascii(seq_slice);

    let hasher = NtHasher::<false>::new(sketcher.k);
    let mut positions = Vec::new();
    let mut super_kmers = Vec::new();

    use simd_minimizers::private::minimizers::minimizers_seq_simd;

    minimizers_seq_simd(packed_seq.as_slice(), &hasher, sketcher.w, &mut sketcher.cache)
        .collect_and_dedup_with_index_into(&mut positions, &mut super_kmers);

    let mut values = Vec::with_capacity(positions.len());
    for &pos in &positions {
        let val = packed_seq.as_slice().read_kmer(sketcher.k, pos as usize);
        values.push(val);
    }

    let pos_ptr = positions.as_mut_ptr();
    let pos_len = positions.len();
    let pos_cap = positions.capacity();
    std::mem::forget(positions);

    let val_ptr = values.as_mut_ptr();
    let val_len = values.len();
    let val_cap = values.capacity();
    std::mem::forget(values);

    let sk_ptr = super_kmers.as_mut_ptr();
    let sk_len = super_kmers.len();
    let sk_cap = super_kmers.capacity();
    std::mem::forget(super_kmers);

    MinimizerResult {
        positions: MinimizerList { data: pos_ptr, len: pos_len, capacity: pos_cap },
        values: MinimizerValues { data: val_ptr, len: val_len, capacity: val_cap },
        super_kmers: SuperKmerList { data: sk_ptr, len: sk_len, capacity: sk_cap },
    }
}

// =============================================================================
// CLOSED SYNCMERS
// =============================================================================
//
// Closed syncmers are k-mers where the minimum s-mer (sub-k-mer) occurs at
// either the first or last position within the k-mer. This provides a
// deterministic and context-independent selection of k-mers.

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
#[unsafe(no_mangle)]
pub unsafe extern "C" fn canonical_syncmer_positions(
    sketcher: *mut SimdSketcher,
    seq: *const c_char,
    len: size_t,
) -> SyncmerList {
    if sketcher.is_null() {
        return SyncmerList { data: std::ptr::null_mut(), len: 0, capacity: 0 };
    }

    let sketcher = unsafe { &mut *sketcher };
    let seq_slice = unsafe { slice::from_raw_parts(seq as *const u8, len) };
    let packed_seq = PackedSeqVec::from_ascii(seq_slice);

    let mut positions = Vec::new();

    use simd_minimizers::canonical_closed_syncmers;
    canonical_closed_syncmers(sketcher.k, sketcher.w)
        .run(packed_seq.as_slice(), &mut positions);

    let ptr = positions.as_mut_ptr();
    let len = positions.len();
    let cap = positions.capacity();
    std::mem::forget(positions);

    SyncmerList { data: ptr, len, capacity: cap }
}

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
#[unsafe(no_mangle)]
pub unsafe extern "C" fn canonical_syncmers(
    sketcher: *mut SimdSketcher,
    seq: *const c_char,
    len: size_t,
) -> SyncmerResult {
    if sketcher.is_null() {
        return SyncmerResult {
            positions: SyncmerList { data: std::ptr::null_mut(), len: 0, capacity: 0 },
            values: SyncmerValues { data: std::ptr::null_mut(), len: 0, capacity: 0 },
            super_kmers: SuperKmerList { data: std::ptr::null_mut(), len: 0, capacity: 0 },
        };
    }

    let sketcher = unsafe { &mut *sketcher };
    let seq_slice = unsafe { slice::from_raw_parts(seq as *const u8, len) };
    let packed_seq = PackedSeqVec::from_ascii(seq_slice);

    let mut positions = Vec::new();

    use simd_minimizers::canonical_closed_syncmers;
    let values: Vec<u64> = canonical_closed_syncmers(sketcher.k, sketcher.w)
        .run(packed_seq.as_slice(), &mut positions)
        .values_u64()
        .collect();

    let pos_ptr = positions.as_mut_ptr();
    let pos_len = positions.len();
    let pos_cap = positions.capacity();
    std::mem::forget(positions);

    let mut values = values;
    let val_ptr = values.as_mut_ptr();
    let val_len = values.len();
    let val_cap = values.capacity();
    std::mem::forget(values);

    SyncmerResult {
        positions: SyncmerList { data: pos_ptr, len: pos_len, capacity: pos_cap },
        values: SyncmerValues { data: val_ptr, len: val_len, capacity: val_cap },
        super_kmers: SuperKmerList { data: std::ptr::null_mut(), len: 0, capacity: 0 },
    }
}

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
#[unsafe(no_mangle)]
pub unsafe extern "C" fn syncmer_positions(
    sketcher: *mut SimdSketcher,
    seq: *const c_char,
    len: size_t,
) -> SyncmerList {
    if sketcher.is_null() {
        return SyncmerList { data: std::ptr::null_mut(), len: 0, capacity: 0 };
    }

    let sketcher = unsafe { &mut *sketcher };
    let seq_slice = unsafe { slice::from_raw_parts(seq as *const u8, len) };
    let packed_seq = PackedSeqVec::from_ascii(seq_slice);

    let mut positions = Vec::new();

    use simd_minimizers::closed_syncmers;
    closed_syncmers(sketcher.k, sketcher.w)
        .run(packed_seq.as_slice(), &mut positions);

    let ptr = positions.as_mut_ptr();
    let len = positions.len();
    let cap = positions.capacity();
    std::mem::forget(positions);

    SyncmerList { data: ptr, len, capacity: cap }
}

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
#[unsafe(no_mangle)]
pub unsafe extern "C" fn syncmers(
    sketcher: *mut SimdSketcher,
    seq: *const c_char,
    len: size_t,
) -> SyncmerResult {
    if sketcher.is_null() {
        return SyncmerResult {
            positions: SyncmerList { data: std::ptr::null_mut(), len: 0, capacity: 0 },
            values: SyncmerValues { data: std::ptr::null_mut(), len: 0, capacity: 0 },
            super_kmers: SuperKmerList { data: std::ptr::null_mut(), len: 0, capacity: 0 },
        };
    }

    let sketcher = unsafe { &mut *sketcher };
    let seq_slice = unsafe { slice::from_raw_parts(seq as *const u8, len) };
    let packed_seq = PackedSeqVec::from_ascii(seq_slice);

    let mut positions = Vec::new();

    use simd_minimizers::closed_syncmers;
    let values: Vec<u64> = closed_syncmers(sketcher.k, sketcher.w)
        .run(packed_seq.as_slice(), &mut positions)
        .values_u64()
        .collect();

    let pos_ptr = positions.as_mut_ptr();
    let pos_len = positions.len();
    let pos_cap = positions.capacity();
    std::mem::forget(positions);

    let mut values = values;
    let val_ptr = values.as_mut_ptr();
    let val_len = values.len();
    let val_cap = values.capacity();
    std::mem::forget(values);

    SyncmerResult {
        positions: SyncmerList { data: pos_ptr, len: pos_len, capacity: pos_cap },
        values: SyncmerValues { data: val_ptr, len: val_len, capacity: val_cap },
        super_kmers: SuperKmerList { data: std::ptr::null_mut(), len: 0, capacity: 0 },
    }
}

// =============================================================================
// OPEN SYNCMERS
// =============================================================================
//
// Open syncmers are k-mers where the minimum s-mer (sub-k-mer) occurs at a
// specific offset position within the k-mer (typically the middle). This
// provides different density characteristics compared to closed syncmers.

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
#[unsafe(no_mangle)]
pub unsafe extern "C" fn canonical_open_syncmer_positions(
    sketcher: *mut SimdSketcher,
    seq: *const c_char,
    len: size_t,
) -> SyncmerList {
    if sketcher.is_null() {
        return SyncmerList { data: std::ptr::null_mut(), len: 0, capacity: 0 };
    }

    let sketcher = unsafe { &mut *sketcher };
    let seq_slice = unsafe { slice::from_raw_parts(seq as *const u8, len) };
    let packed_seq = PackedSeqVec::from_ascii(seq_slice);

    let mut positions = Vec::new();

    use simd_minimizers::canonical_open_syncmers;
    canonical_open_syncmers(sketcher.k, sketcher.w)
        .run(packed_seq.as_slice(), &mut positions);

    let ptr = positions.as_mut_ptr();
    let len = positions.len();
    let cap = positions.capacity();
    std::mem::forget(positions);

    SyncmerList { data: ptr, len, capacity: cap }
}

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
#[unsafe(no_mangle)]
pub unsafe extern "C" fn canonical_open_syncmers(
    sketcher: *mut SimdSketcher,
    seq: *const c_char,
    len: size_t,
) -> SyncmerResult {
    if sketcher.is_null() {
        return SyncmerResult {
            positions: SyncmerList { data: std::ptr::null_mut(), len: 0, capacity: 0 },
            values: SyncmerValues { data: std::ptr::null_mut(), len: 0, capacity: 0 },
            super_kmers: SuperKmerList { data: std::ptr::null_mut(), len: 0, capacity: 0 },
        };
    }

    let sketcher = unsafe { &mut *sketcher };
    let seq_slice = unsafe { slice::from_raw_parts(seq as *const u8, len) };
    let packed_seq = PackedSeqVec::from_ascii(seq_slice);

    let mut positions = Vec::new();

    use simd_minimizers::canonical_open_syncmers;
    let values: Vec<u64> = canonical_open_syncmers(sketcher.k, sketcher.w)
        .run(packed_seq.as_slice(), &mut positions)
        .values_u64()
        .collect();

    let pos_ptr = positions.as_mut_ptr();
    let pos_len = positions.len();
    let pos_cap = positions.capacity();
    std::mem::forget(positions);

    let mut values = values;
    let val_ptr = values.as_mut_ptr();
    let val_len = values.len();
    let val_cap = values.capacity();
    std::mem::forget(values);

    SyncmerResult {
        positions: SyncmerList { data: pos_ptr, len: pos_len, capacity: pos_cap },
        values: SyncmerValues { data: val_ptr, len: val_len, capacity: val_cap },
        super_kmers: SuperKmerList { data: std::ptr::null_mut(), len: 0, capacity: 0 },
    }
}

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
#[unsafe(no_mangle)]
pub unsafe extern "C" fn open_syncmer_positions(
    sketcher: *mut SimdSketcher,
    seq: *const c_char,
    len: size_t,
) -> SyncmerList {
    if sketcher.is_null() {
        return SyncmerList { data: std::ptr::null_mut(), len: 0, capacity: 0 };
    }

    let sketcher = unsafe { &mut *sketcher };
    let seq_slice = unsafe { slice::from_raw_parts(seq as *const u8, len) };
    let packed_seq = PackedSeqVec::from_ascii(seq_slice);

    let mut positions = Vec::new();

    use simd_minimizers::open_syncmers;
    open_syncmers(sketcher.k, sketcher.w)
        .run(packed_seq.as_slice(), &mut positions);

    let ptr = positions.as_mut_ptr();
    let len = positions.len();
    let cap = positions.capacity();
    std::mem::forget(positions);

    SyncmerList { data: ptr, len, capacity: cap }
}

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
#[unsafe(no_mangle)]
pub unsafe extern "C" fn open_syncmers(
    sketcher: *mut SimdSketcher,
    seq: *const c_char,
    len: size_t,
) -> SyncmerResult {
    if sketcher.is_null() {
        return SyncmerResult {
            positions: SyncmerList { data: std::ptr::null_mut(), len: 0, capacity: 0 },
            values: SyncmerValues { data: std::ptr::null_mut(), len: 0, capacity: 0 },
            super_kmers: SuperKmerList { data: std::ptr::null_mut(), len: 0, capacity: 0 },
        };
    }

    let sketcher = unsafe { &mut *sketcher };
    let seq_slice = unsafe { slice::from_raw_parts(seq as *const u8, len) };
    let packed_seq = PackedSeqVec::from_ascii(seq_slice);

    let mut positions = Vec::new();

    use simd_minimizers::open_syncmers;
    let values: Vec<u64> = open_syncmers(sketcher.k, sketcher.w)
        .run(packed_seq.as_slice(), &mut positions)
        .values_u64()
        .collect();

    let pos_ptr = positions.as_mut_ptr();
    let pos_len = positions.len();
    let pos_cap = positions.capacity();
    std::mem::forget(positions);

    let mut values = values;
    let val_ptr = values.as_mut_ptr();
    let val_len = values.len();
    let val_cap = values.capacity();
    std::mem::forget(values);

    SyncmerResult {
        positions: SyncmerList { data: pos_ptr, len: pos_len, capacity: pos_cap },
        values: SyncmerValues { data: val_ptr, len: val_len, capacity: val_cap },
        super_kmers: SuperKmerList { data: std::ptr::null_mut(), len: 0, capacity: 0 },
    }
}

// =============================================================================
// MEMORY MANAGEMENT
// =============================================================================
//
// All results returned by minimizer/syncmer functions allocate memory that must
// be explicitly freed using these functions. Failure to free results in memory leaks.

/// Frees the memory associated with a MinimizerList.
///
/// # Safety
/// - `list` must be a valid MinimizerList returned by a minimizer positions function.
/// - Must not be called twice on the same list.
/// - After calling, the list's data pointer is invalid and must not be accessed.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn free_minimizer_list(list: MinimizerList) {
    if !list.data.is_null() {
        unsafe {
            let _ = Vec::from_raw_parts(list.data, list.len, list.capacity);
        }
    }
}

/// Frees the memory associated with a MinimizerResult.
///
/// Frees all three arrays (positions, values, super_kmers) contained in the result.
///
/// # Safety
/// - `result` must be a valid MinimizerResult returned by a minimizer function.
/// - Must not be called twice on the same result.
/// - After calling, all data pointers in the result are invalid and must not be accessed.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn free_minimizer_result(result: MinimizerResult) {
    if !result.positions.data.is_null() {
        unsafe {
            let _ = Vec::from_raw_parts(result.positions.data, result.positions.len, result.positions.capacity);
        }
    }
    if !result.values.data.is_null() {
        unsafe {
            let _ = Vec::from_raw_parts(result.values.data, result.values.len, result.values.capacity);
        }
    }
    if !result.super_kmers.data.is_null() {
        unsafe {
            let _ = Vec::from_raw_parts(result.super_kmers.data, result.super_kmers.len, result.super_kmers.capacity);
        }
    }
}

/// Frees the memory associated with a SyncmerList.
///
/// This is an alias for `free_minimizer_list` since SyncmerList and MinimizerList
/// have the same structure.
///
/// # Safety
/// - `list` must be a valid SyncmerList returned by a syncmer positions function.
/// - Must not be called twice on the same list.
/// - After calling, the list's data pointer is invalid and must not be accessed.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn free_syncmer_list(list: SyncmerList) {
    unsafe { free_minimizer_list(list) }
}

/// Frees the memory associated with a SyncmerResult.
///
/// This is an alias for `free_minimizer_result` since SyncmerResult and MinimizerResult
/// have the same structure.
///
/// # Safety
/// - `result` must be a valid SyncmerResult returned by a syncmer function.
/// - Must not be called twice on the same result.
/// - After calling, all data pointers in the result are invalid and must not be accessed.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn free_syncmer_result(result: SyncmerResult) {
    unsafe { free_minimizer_result(result) }
}
