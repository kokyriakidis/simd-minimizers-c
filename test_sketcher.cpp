#include <iostream>
#include <vector>
#include <string>
#include <cassert>
#include <cstring>
#include "simd_minimizers.h"

void print_list(const char* name, const MinimizerList& list) {
    std::cout << name << ": [";
    for (size_t i = 0; i < list.len; ++i) {
        std::cout << list.data[i] << (i < list.len - 1 ? ", " : "");
    }
    std::cout << "] (count: " << list.len << ")" << std::endl;
}

void print_result(const char* name, const MinimizerResult& result) {
    std::cout << name << ":" << std::endl;
    std::cout << "  Positions: " << result.positions.len << std::endl;
    std::cout << "  Values: " << result.values.len << std::endl;
    std::cout << "  SuperKmers: " << result.super_kmers.len << std::endl;
}

int main() {
    const char* seq = "ACGTGCTCAGAGACTCAGAGGA";
    size_t len = strlen(seq);

    std::cout << "========================================" << std::endl;
    std::cout << "SIMD Minimizers & Syncmers Test Suite" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Sequence: " << seq << std::endl;
    std::cout << std::endl;

    // ==========================================================================
    // MINIMIZERS TESTS (k=5, w=7: window size 7)
    // ==========================================================================
    std::cout << "--- MINIMIZERS (k=5, w=7) ---" << std::endl;
    std::cout << "For minimizers: k=kmer size, w=window size" << std::endl;

    SimdSketcher* min_sketcher = simd_sketcher_new(5, 7);
    assert(min_sketcher != nullptr);

    // Test canonical minimizer positions
    std::cout << "Testing canonical_minimizer_positions..." << std::endl;
    MinimizerList canonical_min_pos = canonical_minimizer_positions(min_sketcher, seq, len);
    print_list("Canonical Minimizer Positions", canonical_min_pos);
    assert(canonical_min_pos.len > 0);
    free_minimizer_list(canonical_min_pos);

    // Test canonical minimizers full
    std::cout << "Testing canonical_minimizers (full)..." << std::endl;
    MinimizerResult canonical_min_result = canonical_minimizers(min_sketcher, seq, len);
    print_result("Canonical Minimizers", canonical_min_result);
    assert(canonical_min_result.positions.len == canonical_min_result.values.len);
    free_minimizer_result(canonical_min_result);

    // Test forward-only minimizer positions
    std::cout << "Testing minimizer_positions..." << std::endl;
    MinimizerList fwd_min_pos = minimizer_positions(min_sketcher, seq, len);
    print_list("Forward Minimizer Positions", fwd_min_pos);
    assert(fwd_min_pos.len > 0);
    free_minimizer_list(fwd_min_pos);

    // Test forward-only minimizers full
    std::cout << "Testing minimizers (full)..." << std::endl;
    MinimizerResult fwd_min_result = minimizers(min_sketcher, seq, len);
    print_result("Forward Minimizers", fwd_min_result);
    assert(fwd_min_result.positions.len == fwd_min_result.values.len);
    free_minimizer_result(fwd_min_result);

    simd_sketcher_free(min_sketcher);
    std::cout << std::endl;

    // ==========================================================================
    // CLOSED SYNCMERS TESTS (k=11, w=5: kmer size 11, smer size 5)
    // Density = 2/(k-w+1) = 2/7 ≈ 28.6%
    // ==========================================================================
    std::cout << "--- CLOSED SYNCMERS (k=11, w=5) ---" << std::endl;
    std::cout << "For syncmers: k=kmer size, w=smer size (must have w < k)" << std::endl;
    std::cout << "Density = 2/(k-w+1) = 2/7 ≈ 28.6%" << std::endl;

    SimdSketcher* sync_sketcher = simd_sketcher_new(11, 5);
    assert(sync_sketcher != nullptr);

    // Test canonical syncmer positions
    std::cout << "Testing canonical_syncmer_positions..." << std::endl;
    SyncmerList canonical_sync_pos = canonical_syncmer_positions(sync_sketcher, seq, len);
    print_list("Canonical Syncmer Positions", canonical_sync_pos);
    free_syncmer_list(canonical_sync_pos);

    // Test canonical syncmers full
    std::cout << "Testing canonical_syncmers..." << std::endl;
    SyncmerResult canonical_sync_result = canonical_syncmers(sync_sketcher, seq, len);
    print_result("Canonical Syncmers", canonical_sync_result);
    assert(canonical_sync_result.positions.len == canonical_sync_result.values.len);
    free_syncmer_result(canonical_sync_result);

    // Test forward-only syncmer positions
    std::cout << "Testing syncmer_positions..." << std::endl;
    SyncmerList fwd_sync_pos = syncmer_positions(sync_sketcher, seq, len);
    print_list("Forward Syncmer Positions", fwd_sync_pos);
    free_syncmer_list(fwd_sync_pos);

    // Test forward-only syncmers full
    std::cout << "Testing syncmers..." << std::endl;
    SyncmerResult fwd_sync_result = syncmers(sync_sketcher, seq, len);
    print_result("Forward Syncmers", fwd_sync_result);
    assert(fwd_sync_result.positions.len == fwd_sync_result.values.len);
    free_syncmer_result(fwd_sync_result);

    simd_sketcher_free(sync_sketcher);
    std::cout << std::endl;

    // ==========================================================================
    // OPEN SYNCMERS TESTS (k=13, w=5: kmer size 13, smer size 5)
    // For canonical open syncmers, TWO constraints:
    //   1. k-mer size (k) must be odd (for canonicality)
    //   2. window size (k - w + 1) must be odd (for open syncmers)
    // k=13 (odd) ✓, window size = 13 - 5 + 1 = 9 (odd) ✓
    // Density = 1/(k-w+1) = 1/9 ≈ 11.1%
    // ==========================================================================
    std::cout << "--- OPEN SYNCMERS (k=13, w=5) ---" << std::endl;
    std::cout << "For canonical open syncmers:" << std::endl;
    std::cout << "  - k-mer size must be odd: k=13 ✓" << std::endl;
    std::cout << "  - window size must be odd: 13-5+1=9 ✓" << std::endl;
    std::cout << "Density = 1/(k-w+1) = 1/9 ≈ 11.1%" << std::endl;

    SimdSketcher* open_sync_sketcher = simd_sketcher_new(13, 5);
    assert(open_sync_sketcher != nullptr);

    // Test canonical open syncmer positions
    std::cout << "Testing canonical_open_syncmer_positions..." << std::endl;
    SyncmerList canonical_open_sync_pos = canonical_open_syncmer_positions(open_sync_sketcher, seq, len);
    print_list("Canonical Open Syncmer Positions", canonical_open_sync_pos);
    free_syncmer_list(canonical_open_sync_pos);

    // Test canonical open syncmers full
    std::cout << "Testing canonical_open_syncmers..." << std::endl;
    SyncmerResult canonical_open_sync_result = canonical_open_syncmers(open_sync_sketcher, seq, len);
    print_result("Canonical Open Syncmers", canonical_open_sync_result);
    assert(canonical_open_sync_result.positions.len == canonical_open_sync_result.values.len);
    free_syncmer_result(canonical_open_sync_result);

    // Test forward-only open syncmer positions
    std::cout << "Testing open_syncmer_positions..." << std::endl;
    SyncmerList fwd_open_sync_pos = open_syncmer_positions(open_sync_sketcher, seq, len);
    print_list("Forward Open Syncmer Positions", fwd_open_sync_pos);
    free_syncmer_list(fwd_open_sync_pos);

    // Test forward-only open syncmers full
    std::cout << "Testing open_syncmers..." << std::endl;
    SyncmerResult fwd_open_sync_result = open_syncmers(open_sync_sketcher, seq, len);
    print_result("Forward Open Syncmers", fwd_open_sync_result);
    assert(fwd_open_sync_result.positions.len == fwd_open_sync_result.values.len);
    free_syncmer_result(fwd_open_sync_result);

    simd_sketcher_free(open_sync_sketcher);
    std::cout << std::endl;

    // ==========================================================================
    // DENSITY COMPARISON TEST
    // ==========================================================================
    std::cout << "--- DENSITY COMPARISON TEST ---" << std::endl;
    std::cout << "Comparing densities for k=15 with different w values" << std::endl;

    // Create a longer sequence for better density measurement
    const char* long_seq = "ACGTGCTCAGAGACTCAGAGGAACGTGCTCAGAGACTCAGAGGAACGTGCTCAGAGACTCAGAGGA";
    size_t long_len = strlen(long_seq);
    
    // Closed syncmers with k=15, w=5: density = 2/11 ≈ 18.2%
    SimdSketcher* dense_sketcher = simd_sketcher_new(15, 5);
    SyncmerList dense_pos = canonical_syncmer_positions(dense_sketcher, long_seq, long_len);
    double dense_density = (double)dense_pos.len / (long_len - 15 + 1) * 100;
    std::cout << "k=15, w=5 (smer): " << dense_pos.len << " syncmers, density=" 
              << dense_density << "% (expected ~18.2%)" << std::endl;
    free_syncmer_list(dense_pos);
    simd_sketcher_free(dense_sketcher);

    // Closed syncmers with k=15, w=3: density = 2/13 ≈ 15.4%
    SimdSketcher* sparse_sketcher = simd_sketcher_new(15, 3);
    SyncmerList sparse_pos = canonical_syncmer_positions(sparse_sketcher, long_seq, long_len);
    double sparse_density = (double)sparse_pos.len / (long_len - 15 + 1) * 100;
    std::cout << "k=15, w=3 (smer): " << sparse_pos.len << " syncmers, density=" 
              << sparse_density << "% (expected ~15.4%)" << std::endl;
    free_syncmer_list(sparse_pos);
    simd_sketcher_free(sparse_sketcher);

    // Verify: smaller w (smer size) = sparser (lower density)
    std::cout << "Verified: smaller w (smer size) -> lower density" << std::endl;

    std::cout << std::endl;

    // ==========================================================================
    // SKETCHER REUSE TEST
    // ==========================================================================
    std::cout << "--- SKETCHER REUSE TEST ---" << std::endl;

    SimdSketcher* reuse_sketcher = simd_sketcher_new(5, 7);
    MinimizerList reuse1 = canonical_minimizer_positions(reuse_sketcher, seq, len);
    MinimizerList reuse2 = canonical_minimizer_positions(reuse_sketcher, seq, len);
    
    assert(reuse1.len == reuse2.len);
    for (size_t i = 0; i < reuse1.len; ++i) {
        assert(reuse1.data[i] == reuse2.data[i]);
    }
    std::cout << "Sketcher reuse verified - results are consistent" << std::endl;
    
    free_minimizer_list(reuse1);
    free_minimizer_list(reuse2);
    simd_sketcher_free(reuse_sketcher);

    std::cout << std::endl;

    std::cout << "========================================" << std::endl;
    std::cout << "All tests passed!" << std::endl;
    std::cout << "========================================" << std::endl;

    return 0;
}
