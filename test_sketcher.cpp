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
    uint8_t k = 5;
    uint8_t w = 7;

    std::cout << "========================================" << std::endl;
    std::cout << "SIMD Minimizers & Syncmers Test Suite" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Sequence: " << seq << std::endl;
    std::cout << "k=" << (int)k << ", w=" << (int)w << std::endl;
    std::cout << std::endl;

    // Create sketcher
    SimdSketcher* sketcher = simd_sketcher_new(k, w);
    assert(sketcher != nullptr);

    // ==========================================================================
    // MINIMIZERS TESTS
    // ==========================================================================
    std::cout << "--- MINIMIZERS ---" << std::endl;

    // Test canonical minimizer positions
    std::cout << "Testing canonical_minimizer_positions..." << std::endl;
    MinimizerList canonical_min_pos = canonical_minimizer_positions(sketcher, seq, len);
    print_list("Canonical Minimizer Positions", canonical_min_pos);
    assert(canonical_min_pos.len > 0);
    free_minimizer_list(canonical_min_pos);

    // Test canonical minimizers full
    std::cout << "Testing canonical_minimizers (full)..." << std::endl;
    MinimizerResult canonical_min_result = canonical_minimizers(sketcher, seq, len);
    print_result("Canonical Minimizers", canonical_min_result);
    assert(canonical_min_result.positions.len == canonical_min_result.values.len);
    free_minimizer_result(canonical_min_result);

    // Test forward-only minimizer positions
    std::cout << "Testing minimizer_positions..." << std::endl;
    MinimizerList fwd_min_pos = minimizer_positions(sketcher, seq, len);
    print_list("Forward Minimizer Positions", fwd_min_pos);
    assert(fwd_min_pos.len > 0);
    free_minimizer_list(fwd_min_pos);

    // Test forward-only minimizers full
    std::cout << "Testing minimizers (full)..." << std::endl;
    MinimizerResult fwd_min_result = minimizers(sketcher, seq, len);
    print_result("Forward Minimizers", fwd_min_result);
    assert(fwd_min_result.positions.len == fwd_min_result.values.len);
    free_minimizer_result(fwd_min_result);

    std::cout << std::endl;

    // ==========================================================================
    // SYNCMERS TESTS
    // ==========================================================================
    std::cout << "--- SYNCMERS ---" << std::endl;

    // Test canonical syncmer positions
    std::cout << "Testing canonical_syncmer_positions..." << std::endl;
    SyncmerList canonical_sync_pos = canonical_syncmer_positions(sketcher, seq, len);
    print_list("Canonical Syncmer Positions", canonical_sync_pos);
    free_syncmer_list(canonical_sync_pos);

    // Test canonical syncmers full
    std::cout << "Testing canonical_syncmers..." << std::endl;
    SyncmerResult canonical_sync_result = canonical_syncmers(sketcher, seq, len);
    print_result("Canonical Syncmers", canonical_sync_result);
    assert(canonical_sync_result.positions.len == canonical_sync_result.values.len);
    free_syncmer_result(canonical_sync_result);

    // Test forward-only syncmer positions
    std::cout << "Testing syncmer_positions..." << std::endl;
    SyncmerList fwd_sync_pos = syncmer_positions(sketcher, seq, len);
    print_list("Forward Syncmer Positions", fwd_sync_pos);
    free_syncmer_list(fwd_sync_pos);

    // Test forward-only syncmers full
    std::cout << "Testing syncmers..." << std::endl;
    SyncmerResult fwd_sync_result = syncmers(sketcher, seq, len);
    print_result("Forward Syncmers", fwd_sync_result);
    assert(fwd_sync_result.positions.len == fwd_sync_result.values.len);
    free_syncmer_result(fwd_sync_result);

    std::cout << std::endl;

    // ==========================================================================
    // OPEN SYNCMERS TESTS
    // ==========================================================================
    std::cout << "--- OPEN SYNCMERS ---" << std::endl;

    // Test canonical open syncmer positions
    std::cout << "Testing canonical_open_syncmer_positions..." << std::endl;
    SyncmerList canonical_open_sync_pos = canonical_open_syncmer_positions(sketcher, seq, len);
    print_list("Canonical Open Syncmer Positions", canonical_open_sync_pos);
    free_syncmer_list(canonical_open_sync_pos);

    // Test canonical open syncmers full
    std::cout << "Testing canonical_open_syncmers..." << std::endl;
    SyncmerResult canonical_open_sync_result = canonical_open_syncmers(sketcher, seq, len);
    print_result("Canonical Open Syncmers", canonical_open_sync_result);
    assert(canonical_open_sync_result.positions.len == canonical_open_sync_result.values.len);
    free_syncmer_result(canonical_open_sync_result);

    // Test forward-only open syncmer positions
    std::cout << "Testing open_syncmer_positions..." << std::endl;
    SyncmerList fwd_open_sync_pos = open_syncmer_positions(sketcher, seq, len);
    print_list("Forward Open Syncmer Positions", fwd_open_sync_pos);
    free_syncmer_list(fwd_open_sync_pos);

    // Test forward-only open syncmers full
    std::cout << "Testing open_syncmers..." << std::endl;
    SyncmerResult fwd_open_sync_result = open_syncmers(sketcher, seq, len);
    print_result("Forward Open Syncmers", fwd_open_sync_result);
    assert(fwd_open_sync_result.positions.len == fwd_open_sync_result.values.len);
    free_syncmer_result(fwd_open_sync_result);

    std::cout << std::endl;

    // ==========================================================================
    // SKETCHER REUSE TEST
    // ==========================================================================
    std::cout << "--- SKETCHER REUSE TEST ---" << std::endl;

    // Test reuse of sketcher
    std::cout << "Testing sketcher reuse..." << std::endl;
    MinimizerList reuse1 = canonical_minimizer_positions(sketcher, seq, len);
    MinimizerList reuse2 = canonical_minimizer_positions(sketcher, seq, len);
    
    assert(reuse1.len == reuse2.len);
    for (size_t i = 0; i < reuse1.len; ++i) {
        assert(reuse1.data[i] == reuse2.data[i]);
    }
    std::cout << "Sketcher reuse verified - results are consistent" << std::endl;
    
    free_minimizer_list(reuse1);
    free_minimizer_list(reuse2);

    std::cout << std::endl;

    // ==========================================================================
    // DIFFERENT PARAMETERS TEST
    // ==========================================================================
    std::cout << "--- DIFFERENT PARAMETERS TEST ---" << std::endl;

    // Create a new sketcher with different parameters
    // Note: For canonical syncmers, window length (k+w-1) must be odd
    SimdSketcher* sketcher2 = simd_sketcher_new(3, 5);  // k=3, w=5 -> window=7 (odd)
    assert(sketcher2 != nullptr);

    std::cout << "Testing with k=3, w=5..." << std::endl;
    MinimizerList small_k_min = canonical_minimizer_positions(sketcher2, seq, len);
    print_list("Minimizers (k=3, w=5)", small_k_min);
    assert(small_k_min.len > 0);
    free_minimizer_list(small_k_min);

    SyncmerList small_k_sync = canonical_syncmer_positions(sketcher2, seq, len);
    print_list("Syncmers (k=3, w=5)", small_k_sync);
    free_syncmer_list(small_k_sync);

    simd_sketcher_free(sketcher2);

    // ==========================================================================
    // CLEANUP
    // ==========================================================================
    simd_sketcher_free(sketcher);

    std::cout << "========================================" << std::endl;
    std::cout << "All tests passed!" << std::endl;
    std::cout << "========================================" << std::endl;

    return 0;
}
