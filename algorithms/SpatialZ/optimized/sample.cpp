#include <cstdint>

extern "C" void sample_expression(
    const float* reference, const std::int64_t* candidates,
    const double* cumulative, const std::int64_t* lengths,
    const std::int64_t* draw_rows, const double* uniforms, float* output,
    std::int64_t rows, std::int64_t genes, std::int64_t stride, int threads) {
    #pragma omp parallel for num_threads(threads) schedule(static)
    for (std::int64_t row = 0; row < rows; ++row) {
        const auto count = lengths[row];
        if (count == 0)
            continue;
        const auto* choices = candidates + row * stride;
        const auto* bounds = cumulative + row * stride;
        const auto* draws = uniforms + draw_rows[row] * genes;
        auto* values = output + row * genes;
        for (std::int64_t gene = 0; gene < genes; ++gene) {
            std::int64_t selected = 0;
            while (selected + 1 < count && draws[gene] >= bounds[selected])
                ++selected;
            values[gene] = reference[choices[selected] * genes + gene];
        }
    }
}
