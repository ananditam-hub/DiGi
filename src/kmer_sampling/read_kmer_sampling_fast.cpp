// Faster, OpenMP-parallel version of read_kmer_sampling.cpp.
//
// Same sampling scheme (MurmurHash3_32, seed 0, over the canonical k-mer string),
// so output matches the original for reads made of A/C/G/T, but:
//   - k-mers are rolled as 2-bit forward/reverse-complement codes (O(1) per base,
//     no substr/reverse-complement string allocations)
//   - counts are keyed by uint64_t instead of std::string
//   - k-mers containing non-ACGT bases (e.g. N) are skipped; lowercase is folded to uppercase
// Requires k <= 32.
//
// Parallelism:
//   - a reader thread parses FASTQ into batches of reads (double-buffered, overlaps I/O)
//   - phase 1, parallel over reads: sampled canonical k-mers go into per-thread,
//     per-shard buckets (shard = hash % NUM_SHARDS)
//   - phase 2, parallel over shards: each shard's count table is owned by one thread,
//     so no locks and no k-mer is stored in more than one table
// Thread count comes from OMP_NUM_THREADS (default: all cores).
//
// Build: g++ -O3 -fopenmp -pthread read_kmer_sampling_fast.cpp -o read_kmer_sampling_fast

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cstring>
#include <cstdint>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <deque>
#include <omp.h>

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

// Memory-safe MurmurHash3 implementation (identical to read_kmer_sampling.cpp)
uint32_t murmur3_32(const uint8_t* key, size_t len, uint32_t seed) {
    const uint32_t c1 = 0xcc9e2d51;
    const uint32_t c2 = 0x1b873593;
    const uint32_t r1 = 15;
    const uint32_t r2 = 13;
    const uint32_t m = 5;
    const uint32_t n = 0xe6546b64;

    uint32_t hash = seed;
    const int nblocks = len / 4;

    for (int i = 0; i < nblocks; i++) {
        uint32_t k;
        std::memcpy(&k, key + (i * 4), sizeof(uint32_t));

        k *= c1;
        k = ROTL32(k, r1);
        k *= c2;

        hash ^= k;
        hash = ROTL32(hash, r2);
        hash = hash * m + n;
    }

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

// 2-bit encoding: A=0, C=1, G=2, T=3 (so integer order == lexicographic string order)
static int8_t ENC[256];
static const char DEC[4] = {'A', 'C', 'G', 'T'};

static void init_encoding() {
    std::memset(ENC, -1, sizeof(ENC));
    ENC['A'] = ENC['a'] = 0;
    ENC['C'] = ENC['c'] = 1;
    ENC['G'] = ENC['g'] = 2;
    ENC['T'] = ENC['t'] = 3;
}

// Write the k-mer encoded in `code` as uppercase ASCII into buf[0..k)
FORCE_INLINE void decode_kmer(uint64_t code, int k, char* buf) {
    for (int i = k - 1; i >= 0; --i) {
        buf[i] = DEC[code & 3];
        code >>= 2;
    }
}

//-----------------------------------------------------------------------------
// Batches of reads passed from the reader thread to the workers

static const size_t BATCH_READS = 1 << 18;
static const int NUM_SHARDS = 1024;

struct Batch {
    std::string seq;               // all read sequences, concatenated
    std::vector<size_t> starts;    // read i is seq[starts[i], starts[i+1])
    size_t num_reads() const { return starts.size() - 1; }
};

// Two-slot handoff: the reader fills one batch while the workers process the other
class BatchQueue {
public:
    BatchQueue() {
        for (auto& b : storage_) free_.push_back(&b);
    }
    Batch* get_free() {
        std::unique_lock<std::mutex> lock(mu_);
        cv_.wait(lock, [&] { return !free_.empty(); });
        Batch* b = free_.front(); free_.pop_front();
        return b;
    }
    void put_free(Batch* b) {
        { std::lock_guard<std::mutex> lock(mu_); free_.push_back(b); }
        cv_.notify_all();
    }
    void put_full(Batch* b) {
        { std::lock_guard<std::mutex> lock(mu_); full_.push_back(b); }
        cv_.notify_all();
    }
    void close() {
        { std::lock_guard<std::mutex> lock(mu_); done_ = true; }
        cv_.notify_all();
    }
    // Returns nullptr once the reader has finished and every batch has been taken
    Batch* get_full() {
        std::unique_lock<std::mutex> lock(mu_);
        cv_.wait(lock, [&] { return !full_.empty() || done_; });
        if (full_.empty()) return nullptr;
        Batch* b = full_.front(); full_.pop_front();
        return b;
    }
private:
    Batch storage_[2];
    std::deque<Batch*> free_, full_;
    bool done_ = false;
    std::mutex mu_;
    std::condition_variable cv_;
};

int main(int argc, char *argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " [fraction] [kmer size] [reads fastq file] [optional output file]\n";
        return 1;
    }

    double fraction = std::stod(argv[1]);
    int k = std::stoi(argv[2]);
    if (k < 1 || k > 32) {
        std::cerr << "Error: kmer size must be between 1 and 32\n";
        return 1;
    }
    std::ifstream fastq_file(argv[3]);

    if (!fastq_file.is_open()) {
        std::cerr << "Error: cannot open input file " << argv[3] << "\n";
        return 1;
    }

    init_encoding();

    double threshold = fraction * (static_cast<double>(4294967295UL) + 1.0);
    const uint64_t mask = (k == 32) ? ~0ULL : ((1ULL << (2 * k)) - 1);
    const int rc_shift = 2 * (k - 1);
    const int num_threads = omp_get_max_threads();

    std::vector<std::unordered_map<uint64_t, uint64_t>> shard_counts(NUM_SHARDS);
    // buckets[t][s]: sampled k-mers found by thread t that belong to shard s
    std::vector<std::vector<std::vector<uint64_t>>> buckets(
        num_threads, std::vector<std::vector<uint64_t>>(NUM_SHARDS));

    // --- Reader thread ---
    BatchQueue queue;
    std::thread reader([&] {
        std::string line;
        uint64_t line_count = 0;
        Batch* b = queue.get_free();
        b->seq.clear(); b->starts.assign(1, 0);
        while (std::getline(fastq_file, line)) {
            // FASTQ format: sequence is every 4th line, starting at line 2
            if (++line_count % 4 != 2) continue;
            b->seq += line;
            b->starts.push_back(b->seq.size());
            if (b->num_reads() == BATCH_READS) {
                queue.put_full(b);
                b = queue.get_free();
                b->seq.clear(); b->starts.assign(1, 0);
            }
        }
        if (b->num_reads() > 0) queue.put_full(b);
        else queue.put_free(b);
        queue.close();
    });

    // --- Workers ---
    while (Batch* b = queue.get_full()) {
        const char* seq = b->seq.data();
        const long n = static_cast<long>(b->num_reads());

        #pragma omp parallel
        {
            // Phase 1: sample k-mers from reads into this thread's shard buckets
            auto& my_buckets = buckets[omp_get_thread_num()];
            char buf[32];

            #pragma omp for schedule(dynamic, 1024)
            for (long r = 0; r < n; ++r) {
                uint64_t fwd = 0, rev = 0;
                int valid = 0;  // number of consecutive ACGT bases ending at the current position
                for (size_t i = b->starts[r]; i < b->starts[r + 1]; ++i) {
                    int base = ENC[static_cast<uint8_t>(seq[i])];
                    if (base < 0) { valid = 0; fwd = rev = 0; continue; }
                    fwd = ((fwd << 2) | base) & mask;
                    rev = (rev >> 2) | (static_cast<uint64_t>(3 - base) << rc_shift);
                    if (++valid < k) continue;

                    // --- 1. Canonical Logic ---
                    uint64_t canonical = std::min(fwd, rev);

                    // --- 2. Sampling Logic (hash the ASCII k-mer, as in the original) ---
                    decode_kmer(canonical, k, buf);
                    uint32_t hash_val = murmur3_32(reinterpret_cast<const uint8_t*>(buf), k, 0);
                    if (static_cast<double>(hash_val) < threshold) {
                        my_buckets[hash_val % NUM_SHARDS].push_back(canonical);
                    }
                }
            }
            // implicit barrier: all buckets are filled

            // Phase 2: each shard is counted by exactly one thread
            #pragma omp for schedule(dynamic, 4)
            for (int s = 0; s < NUM_SHARDS; ++s) {
                auto& counts = shard_counts[s];
                for (int t = 0; t < num_threads; ++t) {
                    for (uint64_t km : buckets[t][s]) ++counts[km];
                    buckets[t][s].clear();
                }
            }
        }
        queue.put_free(b);
    }
    reader.join();

    // --- 3. Output (shards formatted in parallel, written in order) ---
    std::ostream* out_stream = &std::cout;
    std::ofstream file_out;
    if (argc >= 5) {
        file_out.open(argv[4]);
        if (file_out.is_open()) out_stream = &file_out;
    }

    const int OUT_BLOCK = 64;
    std::vector<std::string> text(OUT_BLOCK);
    for (int s0 = 0; s0 < NUM_SHARDS; s0 += OUT_BLOCK) {
        #pragma omp parallel for schedule(dynamic, 1)
        for (int j = 0; j < OUT_BLOCK; ++j) {
            auto& counts = shard_counts[s0 + j];
            std::string& out = text[j];
            out.clear();
            out.reserve(counts.size() * (k + 8));
            char buf[32];
            for (const auto& p : counts) {
                decode_kmer(p.first, k, buf);
                out.append(buf, k);
                out += '\t';
                out += std::to_string(p.second);
                out += '\n';
            }
            std::unordered_map<uint64_t, uint64_t>().swap(counts);  // free as we go
        }
        for (int j = 0; j < OUT_BLOCK; ++j) out_stream->write(text[j].data(), text[j].size());
    }

    return 0;
}
