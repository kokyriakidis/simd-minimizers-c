use libc::{c_char, size_t};
use simd_minimizers::seq_hash::NtHasher;
use simd_minimizers::Cache;
use simd_minimizers::collect::CollectAndDedup;
use simd_minimizers::private::minimizers::canonical_minimizers_seq_simd;
use packed_seq::{PackedSeqVec, SeqVec, Seq};
use std::slice;

/// A list of minimizer positions (indices).
#[repr(C)]
pub struct MinimizerList {
    /// Pointer to the array of positions.
    pub data: *mut u32,
    /// Number of elements in the array.
    pub len: size_t,
    /// Capacity of the allocated array (internal use).
    pub capacity: size_t,
}

/// A list of minimizer values (hashes).
#[repr(C)]
pub struct MinimizerValues {
    /// Pointer to the array of values.
    pub data: *mut u64,
    /// Number of elements in the array.
    pub len: size_t,
    /// Capacity of the allocated array (internal use).
    pub capacity: size_t,
}

/// A list of super-kmer positions.
#[repr(C)]
pub struct SuperKmerList {
    /// Pointer to the array of super-kmer positions.
    pub data: *mut u32,
    /// Number of elements in the array.
    pub len: size_t,
    /// Capacity of the allocated array (internal use).
    pub capacity: size_t,
}

/// Result structure containing positions, values, and super-kmers.
#[repr(C)]
pub struct MinimizerResult {
    pub positions: MinimizerList,
    pub values: MinimizerValues,
    pub super_kmers: SuperKmerList,
}

/// Wrapper struct for the Sketcher to be used via C pointer.
pub struct SimdSketcher {
    k: usize,
    w: usize,
    cache: Cache,
}

/// Creates a new SimdSketcher.
#[unsafe(no_mangle)]
pub extern "C" fn simd_sketcher_new(k: u8, w: u8) -> *mut SimdSketcher {
    let sketcher = SimdSketcher {
        k: k as usize,
        w: w as usize,
        cache: Cache::default(),
    };
    Box::into_raw(Box::new(sketcher))
}

/// Frees the SimdSketcher.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn simd_sketcher_free(sketcher: *mut SimdSketcher) {
    if !sketcher.is_null() {
        unsafe {
            let _ = Box::from_raw(sketcher);
        }
    }
}

/// Computes minimizer positions using the sketcher.
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

    // Create canonical hasher
    let hasher = NtHasher::<true>::new(sketcher.k);

    let mut positions = Vec::new();
    
    // Call the low-level function using the cached buffers
    canonical_minimizers_seq_simd(packed_seq.as_slice(), &hasher, sketcher.w, &mut sketcher.cache)
        .collect_and_dedup_into::<false>(&mut positions);

    let ptr = positions.as_mut_ptr();
    let len = positions.len();
    let cap = positions.capacity();
    std::mem::forget(positions);

    MinimizerList {
        data: ptr,
        len,
        capacity: cap,
    }
}

/// Computes minimizers (positions, values, super-kmers) using the sketcher.
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
    
    // Collect positions and super-kmer indices
    canonical_minimizers_seq_simd(packed_seq.as_slice(), &hasher, sketcher.w, &mut sketcher.cache)
        .collect_and_dedup_with_index_into(&mut positions, &mut super_kmers);

    // Compute values from positions
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

/// Computes forward-only minimizer positions using the sketcher.
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

    // Use forward hasher
    let hasher = NtHasher::<false>::new(sketcher.k);

    let mut positions = Vec::new();
    
    // Call the low-level function using the cached buffers
    use simd_minimizers::private::minimizers::minimizers_seq_simd;
    
    minimizers_seq_simd(packed_seq.as_slice(), &hasher, sketcher.w, &mut sketcher.cache)
        .collect_and_dedup_into::<false>(&mut positions);

    let ptr = positions.as_mut_ptr();
    let len = positions.len();
    let cap = positions.capacity();
    std::mem::forget(positions);

    MinimizerList {
        data: ptr,
        len,
        capacity: cap,
    }
}

/// Computes forward-only minimizers (positions, values, super-kmers) using the sketcher.
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

/// Frees the memory associated with a MinimizerList.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn free_minimizer_list(list: MinimizerList) {
    if !list.data.is_null() {
        unsafe {
            let _ = Vec::from_raw_parts(list.data, list.len, list.capacity);
        }
    }
}

/// Frees the memory associated with a MinimizerResult.
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

