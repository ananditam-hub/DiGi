#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cstring>
#include <cstdint>

//-----------------------------------------------------------------------------
// Platform-specific functions and macros
#if defined(_MSC_VER)
#define FORCE_INLINE __forceinline
#define ROTL32(x,y) _rotl(x,y)
#else
#define FORCE_INLINE inline __attribute__((always_inline))
inline uint32_t rotl32(uint32_t x, int8_t r) {
    return (x << r) | (x >> (32 - r));
}
#define ROTL32(x,y) rotl32(x,y)
#endif

//-----------------------------------------------------------------------------
// Finalization mix
FORCE_INLINE uint32_t fmix32(uint32_t h) {
    h ^= h >> 16;
    h *= 0x85ebca6b;
    h ^= h >> 13;
    h *= 0xc2b2ae35;
    h ^= h >> 16;
    return h;
}

// Memory-safe MurmurHash3 implementation (using memcpy for alignment safety)
uint32_t murmur3_32(const uint8_t* key, size_t len, uint32_t seed) {
    const uint32_t c1 = 0xcc9e2d51;
    const uint32_t c2 = 0x1b873593;
    const uint32_t r1 = 15;
    const uint32_t r2 = 13;
    const uint32_t m = 5;
    const uint32_t n = 0xe6546b64;

    uint32_t hash = seed;
    const int nblocks = len / 4;

    // Process 4-byte chunks safely
    for (int i = 0; i < nblocks; i++) {
        uint32_t k;
        std::memcpy(&k, key + (i * 4), sizeof(uint32_t)); // Prevents unaligned access crash

        k *= c1;
        k = ROTL32(k, r1);
        k *= c2;

        hash ^= k;
        hash = ROTL32(hash, r2);
        hash = hash * m + n;
    }

    // Handle remaining bytes (tail)
    const uint8_t* tail = (const uint8_t*)(key + nblocks * 4);
    uint32_t k1 = 0;
    switch (len & 3) {
        case 3: k1 ^= tail[2] << 16;
        case 2: k1 ^= tail[1] << 8;
        case 1: k1 ^= tail[0];
                k1 *= c1; k1 = ROTL32(k1, r1); k1 *= c2; hash ^= k1;
    };

    hash ^= len;
    hash = fmix32(hash);
    return hash;
}

// DNA Helper: Get the reverse complement of a k-mer
std::string get_reverse_complement(const std::string& kmer) {
    std::string rev = kmer;
    std::reverse(rev.begin(), rev.end());
    for (char& c : rev) {
        switch (c) {
            case 'A': c = 'T'; break;
            case 'C': c = 'G'; break;
            case 'G': c = 'C'; break;
            case 'T': c = 'A'; break;
            case 'a': c = 't'; break;
            case 'c': c = 'g'; break;
            case 'g': c = 'c'; break;
            case 't': c = 'a'; break;
        }
    }
    return rev;
}

// FracMinHash Selection: Returns true if the k-mer should be sampled
bool should_sample(const std::string& kmer, uint32_t seed, double threshold, uint32_t& out_hash) {
    out_hash = murmur3_32(reinterpret_cast<const uint8_t*>(kmer.c_str()), kmer.length(), seed);
    return (static_cast<double>(out_hash) < threshold);
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " [fraction] [kmer size] [reads fastq file] [optional output file]\n";
        return 1;
    }

    double fraction = std::stod(argv[1]);
    int kmer_size = std::stoi(argv[2]);
    std::ifstream fasta_file(argv[3]);
    
    if (!fasta_file.is_open()) {
        std::cerr << "Error: cannot open input file " << argv[3] << "\n";
        return 1;
    }

    // Pre-calculate threshold once
    double threshold = fraction * (static_cast<double>(4294967295UL) + 1.0);

    std::unordered_map<std::string, uint64_t> kmer_counts;
    std::string line;
    int line_count = 0;

    while (std::getline(fasta_file, line)) {
        line_count++;
        // FASTQ format: sequence is every 4th line, starting at line 2
        if (line_count % 4 == 2) {
            if (line.size() < static_cast<size_t>(kmer_size)) continue;

            for (size_t i = 0; i <= line.size() - kmer_size; ++i) {
                std::string kmer = line.substr(i, kmer_size);
                
                // --- 1. Canonical Logic ---
                std::string rev_comp = get_reverse_complement(kmer);
                std::string canonical = (kmer < rev_comp) ? kmer : rev_comp;

                // --- 2. Sampling Logic ---
                uint32_t hash_val;
                if (should_sample(canonical, 0, threshold, hash_val)) {
                    ++kmer_counts[canonical];
                }
            }
        }
    }

    // --- 3. Output ---
    std::ostream* out_stream = &std::cout;
    std::ofstream file_out;
    if (argc >= 5) {
        file_out.open(argv[4]);
        if (file_out.is_open()) out_stream = &file_out;
    }

    for (const auto& p : kmer_counts) {
        *out_stream << p.first << '\t' << p.second << '\n';
    }

    return 0;
}