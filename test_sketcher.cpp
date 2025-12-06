#include <iostream>
#include <vector>
#include <string>
#include <cassert>
#include <cstring>
#include "simd_minimizers.h"

void print_minimizers(const MinimizerList& list) {
    std::cout << "Minimizers: [";
    for (size_t i = 0; i < list.len; ++i) {
        std::cout << list.data[i] << (i < list.len - 1 ? ", " : "");
    }
    std::cout << "]" << std::endl;
}

int main() {
    const char* seq = "ACGTGCTCAGAGACTCAGAGGA";
    size_t len = strlen(seq);
    uint8_t k = 5;
    uint8_t w = 7;

    std::cout << "Sequence: " << seq << std::endl;

    // Test stateful API
    std::cout << "Testing stateful API..." << std::endl;
    SimdSketcher* sketcher = simd_sketcher_new(k, w);
    assert(sketcher != nullptr);

    MinimizerList list2 = canonical_minimizer_positions(sketcher, seq, len);
    print_minimizers(list2);

    // Test reuse
    std::cout << "Testing reuse..." << std::endl;
    MinimizerList list3 = canonical_minimizer_positions(sketcher, seq, len);
    print_minimizers(list3);
    assert(list3.len > 0);
    
    // Verify results match
    assert(list2.len == list3.len);
    for (size_t i = 0; i < list2.len; ++i) {
        assert(list2.data[i] == list3.data[i]);
    }

    free_minimizer_list(list2);
    free_minimizer_list(list3);

    // Test full result
    std::cout << "Testing full result..." << std::endl;
    MinimizerResult result = canonical_minimizers(sketcher, seq, len);
    std::cout << "Positions: " << result.positions.len << std::endl;
    std::cout << "Values: " << result.values.len << std::endl;
    std::cout << "SuperKmers: " << result.super_kmers.len << std::endl;
    
    assert(result.positions.len == result.values.len);
    
    free_minimizer_result(result);

    // Test forward-only API
    std::cout << "Testing forward-only API..." << std::endl;
    MinimizerList fwd_list = minimizer_positions(sketcher, seq, len);
    print_minimizers(fwd_list);
    assert(fwd_list.len > 0);
    free_minimizer_list(fwd_list);

    MinimizerResult fwd_result = minimizers(sketcher, seq, len);
    std::cout << "Forward Positions: " << fwd_result.positions.len << std::endl;
    std::cout << "Forward Values: " << fwd_result.values.len << std::endl;
    assert(fwd_result.positions.len == fwd_result.values.len);
    free_minimizer_result(fwd_result);

    simd_sketcher_free(sketcher);

    std::cout << "All tests passed!" << std::endl;
    return 0;
}
