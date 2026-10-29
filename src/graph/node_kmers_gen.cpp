#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cstring>
#include <cstdint>
#include <unordered_map>

// --- MurmurHash3 Core ---

#if defined(_MSC_VER)
#define ROTL32(x,y) _rotl(x,y)
#else
inline uint32_t rotl32(uint32_t x, int8_t r) {
    return (x << r) | (x >> (32 - r));
}
#define ROTL32(x,y) rotl32(x,y)
#endif

inline uint32_t fmix32(uint32_t h) {
    h ^= h >> 16;
    h *= 0x85ebca6b;
    h ^= h >> 13;
    h *= 0xc2b2ae35;
    h ^= h >> 16;
    return h;
}

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

        std::memcpy(
            &k,
            key + (i * 4),
            sizeof(uint32_t)
        );

        k *= c1;
        k = ROTL32(k, r1);
        k *= c2;

        hash ^= k;
        hash = ROTL32(hash, r2);
        hash = hash * m + n;
    }

    const uint8_t* tail =
        (const uint8_t*)(key + nblocks * 4);

    uint32_t k1 = 0;

    switch (len & 3) {

        case 3:
            k1 ^= tail[2] << 16;

        case 2:
            k1 ^= tail[1] << 8;

        case 1:
            k1 ^= tail[0];
            k1 *= c1;
            k1 = ROTL32(k1, r1);
            k1 *= c2;
            hash ^= k1;
    };

    hash ^= len;

    return fmix32(hash);
}

// --- DNA Helper Functions ---

std::string get_canonical(const std::string& kmer) {

    std::string rev = kmer;

    std::reverse(rev.begin(), rev.end());

    for (char& c : rev) {

        switch (c) {

            case 'A':
            case 'a':
                c = 'T';
                break;

            case 'C':
            case 'c':
                c = 'G';
                break;

            case 'G':
            case 'g':
                c = 'C';
                break;

            case 'T':
            case 't':
                c = 'A';
                break;
        }
    }

    return (kmer < rev) ? kmer : rev;
}

// --- Main Processing Logic ---

int main(int argc, char* argv[]) {

    std::unordered_map<std::string, std::string> args;

    for (int i = 1; i < argc - 1; i += 2) {
        args[argv[i]] = argv[i + 1];
    }

    if (
        args.find("--fraction") == args.end() ||
        args.find("--k") == args.end() ||
        args.find("--input") == args.end()
    ) {

        std::cerr << "\nUsage:\n";

        std::cerr
            << argv[0]
            << " --fraction <float>"
            << " --k <int>"
            << " --input <graph.txt>"
            << " [--output output.txt]\n\n";

        return 1;
    }

    double fraction =
        std::stod(args["--fraction"]);

    int k =
        std::stoi(args["--k"]);

    std::string input_file =
        args["--input"];


    std::ifstream fin(input_file);

    if (!fin.is_open()) {

        std::cerr
            << "Error: Cannot open input file "
            << input_file
            << "\n";

        return 1;
    }
    double threshold =
        fraction * 4294967296.0;

    std::ostream* out = &std::cout;

    std::ofstream fout;

    if (args.find("--output") != args.end()) {

        fout.open(args["--output"]);

        if (fout.is_open()) {
            out = &fout;
        }
    }

   std::string line;

    while (std::getline(fin, line)) {

        if (line.empty() || line[0] == '#')
            continue;

        // input format:
        // <node_id>\t<sequence>

        std::stringstream ss(line);

        std::string node_id;
        std::string seq;

        if (!(ss >> node_id >> seq))
            continue;

        // skip null / tiny sequences

        if (
            seq == "_" ||
            seq.length() < (size_t)k
        ) {
            continue;
        }

        std::vector<std::string> sampled_kmers;
        // Generate k-mers
        for (size_t i = 0;
             i <= seq.length() - k;
             ++i) {

            std::string kmer =
                seq.substr(i, k);
            // skip ambiguous bases
            if (
                kmer.find_first_of("Nn")
                != std::string::npos
            ) {
                continue;
            }
            // canonical k-mer
            std::string canonical =
                get_canonical(kmer);
            // MurmurHash3
            uint32_t h =
                murmur3_32(
                    reinterpret_cast<const uint8_t*>(
                        canonical.c_str()
                    ),
                    k,
                    0
                );
            // FracMinHash sampling
            if (
                static_cast<double>(h)
                < threshold
            ) {
                sampled_kmers.push_back(canonical);
            }
        }

        // Write sampled k-mers

        if (!sampled_kmers.empty()) {

            *out << node_id << "\t";

            for (size_t i = 0;
                 i < sampled_kmers.size();
                 ++i) {

                *out
                    << sampled_kmers[i]
                    << (
                        i == sampled_kmers.size() - 1
                        ? ""
                        : ","
                    );
            }

            *out << "\n";
        }
    }

    fin.close();

    if (fout.is_open()) {
        fout.close();
    }

    return 0;
}