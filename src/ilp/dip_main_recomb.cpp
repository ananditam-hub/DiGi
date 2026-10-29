#include <fstream>
#include <map>
#include <unordered_set>
#include <set>
#include <vector>
#include <iostream>
#include <chrono>
#include <cstdint>
#include <algorithm>
#include <string>
#include <filesystem>
#include <thread>

#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>

#include "gurobi_c++.h"

// --- Graph Structure ---
class Graph{
    public:
        int n;
        std::vector<std::vector<int>> adj;
        std::vector<std::string> label;

        Graph(int n);
        void loadLabels(const std::string& filename);
        void addEdge(int u, int v);
        void print() const;
};

// --- Graph Structure Methods ---
Graph::Graph(int n): n(n), adj(n), label(n) {}
void Graph::addEdge(int u, int v){ adj[u].push_back(v); }

void Graph::loadLabels(const std::string& filename) {
    std::ifstream fin(filename);
    if (!fin.is_open()) return;
    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        size_t space = line.find(' ');
        if (space == std::string::npos) space = line.find('\t');
        if (space == std::string::npos) continue;

        int node_id = std::stoi(line.substr(0, space));
        std::string node_label = line.substr(space + 1);
        if (node_label == "_") node_label = "";
        if (node_id >= 0 && node_id < n) label[node_id] = node_label;
    }
}

void Graph::print() const {
    for (int i = 0; i < n; i++) {
        std::cout << i << ": ";
        for (int v : adj[i]) std::cout << v <<", ";
        std::cout << "\n";
    }
}


// 1. Load Target Read K-mers (Sk) as an unordered_set (ignoring frequencies)
std::unordered_set<std::string> loadSk(const std::string& sk_file) {
    std::unordered_set<std::string> Sk;
    std::ifstream file(sk_file);
    if (!file.is_open()) throw std::runtime_error("Cannot open Sk file: " + sk_file);

    std::string line;
    while(std::getline(file, line)){
        if (line.empty()) continue;
        if (line.back() == '\r') line.pop_back();

        // Target kmers are stored as <kmer>\t<number_of_occurrences>
        // We only extract the kmer and ignore the occurrences
        size_t tab = line.find('\t');
        if (tab != std::string::npos) {
            line = line.substr(0, tab);
        }
        Sk.insert(line);
    }
    return Sk;
}

// 2. Process Graph K-mers in format: <node_id>\t<multiset of kmers>
void processGraphKmers(const std::string& graph_file,
                       const std::unordered_set<std::string>& Sk,
                       std::map<std::string, std::vector<int>>& targetKmerToNodes,
                       std::vector<uint32_t>& novel_kmer_count,
                       int num_nodes)
{
    std::ifstream file(graph_file);
    if (!file.is_open()) throw std::runtime_error("Cannot open Graph Kmers: " + graph_file);

    // Temporary map to collect unique node IDs for each target k-mer across different lines
    std::map<std::string, std::set<int>> temp_target_nodes;

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        if (line.back() == '\r') line.pop_back();

        size_t tab_pos = line.find('\t');
        if (tab_pos == std::string::npos) {
            tab_pos = line.find(' ');
            if (tab_pos == std::string::npos) continue;
        }

        int node_id = -1;
        try {
            node_id = std::stoi(line.substr(0, tab_pos));
        } catch (...) {
            continue;
        }

        if (node_id < 0 || node_id >= num_nodes) continue;

        std::string kmers_str = line.substr(tab_pos + 1);

        // Optional: remove surrounding brackets if present in file format
        if (!kmers_str.empty() && kmers_str.front() == '[') kmers_str.erase(0, 1);
        if (!kmers_str.empty() && kmers_str.back() == ']') kmers_str.pop_back();

        size_t start = 0;
        while (start < kmers_str.size()) {
            size_t comma_pos = kmers_str.find(',', start);
            if (comma_pos == std::string::npos) comma_pos = kmers_str.size();

            if (comma_pos > start) {
                std::string kmer = kmers_str.substr(start, comma_pos - start);
                if (!kmer.empty()) {
                    if (Sk.find(kmer) != Sk.end()) {
                        temp_target_nodes[kmer].insert(node_id);
                    } else {
                        novel_kmer_count[node_id]++; // =1 if no frqeuency needed
                    }
                }
            }
            start = comma_pos + 1;
        }
    }

    // Convert sets to vectors for the final targetKmerToNodes map
    for (const auto& pair : temp_target_nodes) {
        targetKmerToNodes[pair.first].assign(pair.second.begin(), pair.second.end());
    }
}

// ===================== HAPLOTYPE NODE-COORDINATE LOADER =====================
// Parses a debug-log-style file where each "row" of the lattice (one candidate
// reference haplotype) is introduced by a "Processing ... (FASTA header: NAME)"
// line, followed by exactly r+1 "Node k: start=S, end=E, length=L" lines giving
// that haplotype's own coordinate range for lattice column k-1 (0-indexed).
// Haplotype blocks are assumed to appear in the same order as lattice rows
// 0..H-1. This lets us translate a solved path (a sequence of lattice node ids)
// back into "which named haplotype was this segment copied from, and what
// coordinate range on that haplotype does it span".
struct HaplotypeCoords {
    std::string name;                                  // e.g. "HG00096#1#MHC"
    std::vector<std::pair<long long,long long>> seg;    // seg[col] = (start,end) for lattice column `col`
};

std::vector<HaplotypeCoords> loadNodeCoords(const std::string& filename) {
    std::vector<HaplotypeCoords> haps;
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cerr << "Warning: could not open node-coords file: " << filename
                   << " (haplotype-region output will be skipped)\n";
        return haps;
    }

    std::string line;
    while (std::getline(fin, line)) {
        // New haplotype block, e.g.: "... Processing HG00096_1 (FASTA header: HG00096#1#MHC)..."
        size_t proc_pos = line.find("Processing ");
        if (proc_pos != std::string::npos) {
            size_t hdr_pos = line.find("(FASTA header: ", proc_pos);
            if (hdr_pos != std::string::npos) {
                size_t name_start = hdr_pos + std::string("(FASTA header: ").size();
                size_t name_end = line.find(')', name_start);
                if (name_end != std::string::npos) {
                    haps.push_back(HaplotypeCoords{line.substr(name_start, name_end - name_start), {}});
                }
            }
            continue;
        }

        // Node line, e.g.: "... Node 1: start=0, end=20459, length=20459"
        size_t node_pos = line.find("Node ");
        if (node_pos != std::string::npos && !haps.empty()) {
            size_t start_pos = line.find("start=", node_pos);
            size_t end_pos   = line.find("end=", node_pos);
            if (start_pos != std::string::npos && end_pos != std::string::npos) {
                try {
                    long long s = std::stoll(line.substr(start_pos + 6));
                    long long e = std::stoll(line.substr(end_pos + 4));
                    haps.back().seg.emplace_back(s, e);
                } catch (...) {
                    // malformed line; skip it rather than aborting the whole load
                }
            }
        }
    }
    return haps;
}

// Maps a lattice node id back to (row, col). Returns (-1,-1) for the source,
// sink, or a recombination hub node, none of which belong to a single haplotype.
std::pair<int,int> latticeRowCol(int node_id, int H, int r, int sink) {
    if (node_id == 0 || node_id == sink) return {-1, -1};
    int hub_start = H * (r + 1) + 1;
    int hub_end   = H * (r + 1) + r; // inclusive
    if (node_id >= hub_start && node_id <= hub_end) return {-1, -1};
    int idx = node_id - 1;
    int row = idx / (r + 1);
    int col = idx % (r + 1);
    if (row >= H) return {-1, -1};
    return {row, col};
}

// ===================== MAIN =====================
int main(int argc, char *argv[]){
    auto start = std::chrono::high_resolution_clock::now();
    std::string GRAPH_KMERS_FILE = "";
    std::string SK_FILE          = "";
    // c1 = missing target k-mer penalty, c2 = novel k-mer penalty.
    // Recombinations carry no objective cost; they are only limited by the --R budget.
    double c1 = 1.0, c2 = 1.0;
    int H = -1, r = -1;
    int R = -1;
    std::string LABEL_INPUT_FILE = "ref_graph.txt";
    std::string sample_name      = "Sample";
    std::string NODE_COORDS_FILE = "";

    std::string out_dir      = ".";
    std::string cov_label    = "";
    std::string rel_weights  = "";
    std::string lockfile_path = ""; // optional; serializes WLS checkout across concurrent jobs

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--graph_kmers" && i + 1 < argc) GRAPH_KMERS_FILE = argv[++i];
        else if (arg == "--sk" && i + 1 < argc)     SK_FILE = argv[++i];
        else if (arg == "--c1" && i + 1 < argc)     c1 = std::stod(argv[++i]);
        else if (arg == "--c2" && i + 1 < argc)     c2 = std::stod(argv[++i]);
        else if (arg == "--c3") {
            std::cerr << "Error: --c3 was removed; use --c1 (missing k-mer), --c2 (novel k-mer)\n";
            return 1;
        }
        else if (arg == "--cr") {
            std::cerr << "Error: --cr was removed; recombinations are not penalised, only capped by --R\n";
            return 1;
        }
        else if (arg == "--H" && i + 1 < argc)      H = std::stoi(argv[++i]);
        else if (arg == "--r" && i + 1 < argc)      r = std::stoi(argv[++i]);
        else if (arg == "--R" && i + 1 < argc)      R = std::stoi(argv[++i]);           // === NEW ===
        else if (arg == "--labels" && i + 1 < argc) LABEL_INPUT_FILE = argv[++i];
        else if (arg == "--sample" && i + 1 < argc) sample_name = argv[++i];
        else if (arg == "--node_coords" && i + 1 < argc) NODE_COORDS_FILE = argv[++i]; // === NEW ===
        else if (arg == "--out_dir" && i + 1 < argc) out_dir = argv[++i];         // === NEW ===
        else if (arg == "--cov" && i + 1 < argc)     cov_label = argv[++i];        // === NEW ===
        else if (arg == "--rel_weights" && i + 1 < argc) rel_weights = argv[++i];  // === NEW ===
        else if (arg == "--lockfile" && i + 1 < argc) lockfile_path = argv[++i];   // === NEW ===
    }

    if (GRAPH_KMERS_FILE.empty() || SK_FILE.empty() || H < 0 || r < 0) {
        std::cerr << "Usage: " << argv[0] << " --graph_kmers <file> --sk <file> --H <int> --r <int> [--R <int>] [options]\n";
        return 1;
    }

    const int n = H * (r + 1) + r + 2;
    int sink = H * (r + 1) + r + 1;
    Graph g(n);

    // Build Graph
    for(int row = 0; row < H; row++){
        g.addEdge(0, row*(r+1) + 1);
        g.addEdge(row*(r+1) + r + 1, sink);
    }
    for(int row = 0; row < H; row++){
        for(int col = 0; col < r; col++) g.addEdge(row*(r+1) + col + 1, row*(r+1) + col + 2);
    }
    for(int recomb = 0; recomb < r; ++recomb){
        int w_node = H*(r+1) + recomb + 1;
        for(int row = 0; row < H; row++){
            g.addEdge(row*(r+1) + recomb + 1, w_node);
            g.addEdge(w_node, row*(r+1) + recomb + 2);
        }
    }

    g.loadLabels(LABEL_INPUT_FILE);
    std::vector<HaplotypeCoords> node_coords;
    if (!NODE_COORDS_FILE.empty()) {
        node_coords = loadNodeCoords(NODE_COORDS_FILE);
        std::cout << "Loaded node coordinates for " << node_coords.size() << " haplotypes";
        if ((int)node_coords.size() != H) {
            std::cout << " (warning: this does not match H=" << H
                       << " -- rows beyond the loaded haplotypes will be labeled UNKNOWN)";
        }
        std::cout << "\n";
    }

    // --- PIPELINE EXECUTION ---
    std::cout << "Loading Target Kmers...\n";
    auto Sk = loadSk(SK_FILE);

    std::cout << "Processing Graph Kmers in memory...\n";
    std::map<std::string, std::vector<int>> kmerToNodes;
    std::vector<uint32_t> novel_kmer_count(g.n, 0);

    processGraphKmers(GRAPH_KMERS_FILE, Sk, kmerToNodes, novel_kmer_count, g.n);

    // --- GUROBI MODEL ---
    try {
        int lockfd = -1;
        if (!lockfile_path.empty()) {
            std::filesystem::path lock_dir = std::filesystem::path(lockfile_path).parent_path();
            if (!lock_dir.empty()) std::filesystem::create_directories(lock_dir);
            lockfd = open(lockfile_path.c_str(), O_CREAT | O_RDWR, 0666);
            if (lockfd < 0) std::cerr << "Warning: cannot open lockfile " << lockfile_path
                                      << ", continuing without lock\n";
        }

        GRBEnv env(true);
        env.set("LogFile", "dip_ilp.log");

        const int max_attempts = 3;
        for (int attempt = 1; attempt <= max_attempts; ++attempt) {
            if (lockfd >= 0) flock(lockfd, LOCK_EX);
            try {
                env.start();
                if (lockfd >= 0) flock(lockfd, LOCK_UN);
                break;
            } catch (GRBException& e) {
                if (lockfd >= 0) flock(lockfd, LOCK_UN);
                if (attempt == max_attempts) { if (lockfd >= 0) close(lockfd); throw; }
                std::cerr << "WLS checkout failed (attempt " << attempt << "/"
                          << max_attempts << "): " << e.getMessage()
                          << " -- retrying...\n";
                std::this_thread::sleep_for(std::chrono::seconds(5 * attempt * attempt));
            }
        }
        if (lockfd >= 0) close(lockfd);

        GRBModel model(env);
        model.set(GRB_StringAttr_ModelName, "Diploid Path Inference");
        model.set(GRB_IntParam_Threads, 1); // serial solve per run

        std::vector<GRBVar> x_1_vars, x_2_vars;
        std::vector<std::pair<int,int>> edge_list;
        std::vector<std::vector<int>> out_edge_idx1(g.n), in_edge_idx1(g.n);
        std::vector<std::vector<int>> out_edge_idx2(g.n), in_edge_idx2(g.n);
        auto is_hub = [&](int node) { return node > H*(r+1) && node != sink; };
        auto recomb_weight = [&](int u, int v) { return (is_hub(u) || is_hub(v)) ? 0.5 : 0.0; };
        std::vector<double> edge_w;

        GRBQuadExpr obj = 0;

        for(int u = 0; u < g.n; ++u){
            for(int v : g.adj[u]){
                GRBVar var1 = model.addVar(0, 1, 0, GRB_BINARY, "x1_" + std::to_string(u) + "_" + std::to_string(v));
                GRBVar var2 = model.addVar(0, 1, 0, GRB_BINARY, "x2_" + std::to_string(u) + "_" + std::to_string(v));
                x_1_vars.push_back(var1); out_edge_idx1[u].push_back(x_1_vars.size()-1); in_edge_idx1[v].push_back(x_1_vars.size()-1);
                x_2_vars.push_back(var2); out_edge_idx2[u].push_back(x_2_vars.size()-1); in_edge_idx2[v].push_back(x_2_vars.size()-1);
                edge_list.emplace_back(u,v);
                edge_w.push_back(recomb_weight(u, v));

                // Penalty naturally scales with the multiset frequency stored here
                double cost = c2 * novel_kmer_count[v];

                obj += cost * (var1 + var2);
            }
        }

        std::vector<std::string> kmer_list;
        for (const auto& pair : kmerToNodes) kmer_list.push_back(pair.first);

        std::vector<GRBVar> y_1_vars(kmer_list.size()), y_2_vars(kmer_list.size());
        for(size_t i = 0; i < kmer_list.size(); ++i){
            y_1_vars[i] = model.addVar(0, 1, 0, GRB_BINARY, "y1_" + kmer_list[i]);
            y_2_vars[i] = model.addVar(0, 1, 0, GRB_BINARY, "y2_" + kmer_list[i]);

            // Flat penalty c1 for any missing target k-mer
            obj += c1 * (1.0 - y_1_vars[i]) * (1.0 - y_2_vars[i]);
            // obj += c1 * (2.0 - y_1_vars[i] - y_2_vars[i]); // simpler equivalent form
            GRBLinExpr cover1 = 0, cover2 = 0;
            for(int v : kmerToNodes[kmer_list[i]]){
                for(int idx1 : in_edge_idx1[v]) cover1 += x_1_vars[idx1];
                for(int idx2 : in_edge_idx2[v]) cover2 += x_2_vars[idx2];
            }
            model.addConstr(cover1 >= y_1_vars[i]);
            model.addConstr(cover2 >= y_2_vars[i]);
        }

        model.setObjective(obj, GRB_MINIMIZE);

        for(int i = 0; i < g.n; i++){
            GRBLinExpr in1 = 0, out1 = 0, in2 = 0, out2 = 0;
            for(int idx : in_edge_idx1[i]) in1 += x_1_vars[idx];
            for(int idx : out_edge_idx1[i]) out1 += x_1_vars[idx];
            for(int idx : in_edge_idx2[i]) in2 += x_2_vars[idx];
            for(int idx : out_edge_idx2[i]) out2 += x_2_vars[idx];

            if(i == 0) { model.addConstr(out1 - in1 == 1); model.addConstr(out2 - in2 == 1); }
            else if(i == sink) { model.addConstr(out1 - in1 == -1); model.addConstr(out2 - in2 == -1); }
            else { model.addConstr(out1 - in1 == 0); model.addConstr(out2 - in2 == 0); }
        }

        if (R >= 0) {
            GRBLinExpr nrec1 = 0, nrec2 = 0;
            for (size_t idx = 0; idx < edge_list.size(); ++idx) {
                if (edge_w[idx] == 0.0) continue;
                nrec1 += edge_w[idx] * x_1_vars[idx];
                nrec2 += edge_w[idx] * x_2_vars[idx];
            }
            model.addConstr(nrec1 <= R, "max_recomb_1");
            model.addConstr(nrec2 <= R, "max_recomb_2");
            std::cout << "Recombination cap: at most " << R << " per haplotype (W = 1/2 per hub edge)\n";
        }

        std::cout << "Optimizing ILP Model...\n";
        model.set(GRB_IntParam_PoolSearchMode, 0);
        model.set(GRB_IntParam_PoolSolutions, 1);
        model.optimize();
        int status = model.get(GRB_IntAttr_Status);

        if (status != GRB_OPTIMAL && status != GRB_SUBOPTIMAL) {
            std::cout << "No feasible solution found." << std::endl;
        } else {
            int nsol = model.get(GRB_IntAttr_SolCount);
            std::cout << "Found " << nsol << " optimal solutions\n";

            for (int s = 0; s < nsol; ++s) {
                model.set(GRB_IntParam_SolutionNumber, s);
                double obj_val = model.get(GRB_DoubleAttr_PoolObjVal);

                std::cout << "\n========================================\n";
                std::cout << "Solution " << s << "   Objective = " << obj_val << "\n";
                std::cout << "========================================\n";

                std::map<int,int> next_node_1;
                std::map<int,int> next_node_2;

                for (size_t e = 0; e < x_1_vars.size(); ++e) {
                    if (x_1_vars[e].get(GRB_DoubleAttr_Xn) > 0.5) {
                        auto [u, v] = edge_list[e];
                        next_node_1[u] = v;
                    }
                }

                for (size_t e = 0; e < x_2_vars.size(); ++e) {
                    if (x_2_vars[e].get(GRB_DoubleAttr_Xn) > 0.5) {
                        auto [u, v] = edge_list[e];
                        next_node_2[u] = v;
                    }
                }

                auto process_and_save_path = [&](const std::map<int,int>& next_node,
                                                 const std::string& tag,
                                                 const std::string& out_filename) -> int {
                    std::cout << "\n--- Final Path Traversal (" << tag << ") ---\n";

                    int src = 0;
                    int t = sink;
                    int curr = src;
                    int n_recomb = 0;

                    std::vector<int> path_nodes;
                    std::set<int> visited;

                    visited.insert(curr);
                    path_nodes.push_back(curr);

                    if (!next_node.count(src)) {
                        std::cout << "Error: No path from source\n";
                        return -1;
                    }

                    std::cout << curr;

                    while (curr != t) {
                        auto it = next_node.find(curr);
                        if (it == next_node.end()) {
                            std::cout << " -> [BROKEN PATH]";
                            break;
                        }

                        int next = it->second;

                        if (visited.count(next)) {
                            std::cout << " -> " << next << " [CYCLE DETECTED]";
                            break;
                        }

                        std::cout << " -> " << next;

                        if (next > H*(r+1) && next != t) {
                            n_recomb++;
                        }

                        curr = next;
                        visited.insert(curr);
                        path_nodes.push_back(curr);
                    }

                    if (curr == t) {
                        std::cout << "\nPath verification: Reached Sink (" << t << ")";
                        std::cout << "\nNumber of recombinations: " << n_recomb << "\n";

                        std::string sequence;
                        for (int node_id : path_nodes) {
                            sequence += g.label[node_id];
                        }

                        std::ofstream fout(out_filename);
                        if (fout.is_open()) {
                            fout << ">" << sample_name;
                            if (nsol > 1) fout << "_sol" << s;
                            fout << "_" << tag << "_recomb" << n_recomb << "\n";

                            const size_t line_width = 60; // standard FASTA line wrap
                            for (size_t i = 0; i < sequence.size(); i += line_width) {
                                fout << sequence.substr(i, line_width) << "\n";
                            }
                            fout.close();
                            std::cout << "Saved path FASTA to: " << out_filename << "\n";
                        } else {
                            std::cerr << "Error writing to file: " << out_filename << "\n";
                        }
                        if (!node_coords.empty()) {
                            std::string region_filename = out_filename;
                            size_t dot = region_filename.find_last_of('.');
                            if (dot != std::string::npos) region_filename = region_filename.substr(0, dot);
                            region_filename += "_haplotype_regions.tsv";

                            std::ofstream rout(region_filename);
                            if (rout.is_open()) {
                                rout << "start\tend\thaplotype\tgraph_row\tn_columns\n";

                                int seg_row = -2;      // -2 = no open segment yet
                                long long seg_start = -1, seg_end = -1;
                                int seg_cols = 0;

                                auto flush_segment = [&]() {
                                    if (seg_row >= 0) {
                                        std::string hap_name = (seg_row < (int)node_coords.size())
                                            ? node_coords[seg_row].name
                                            : ("row_" + std::to_string(seg_row) + "_UNKNOWN");
                                        rout << seg_start << "\t" << seg_end << "\t" << hap_name
                                             << "\t" << seg_row << "\t" << seg_cols << "\n";
                                    }
                                };

                                for (int node_id : path_nodes) {
                                    auto [row, col] = latticeRowCol(node_id, H, r, t);
                                    if (row < 0) continue; // source / sink / recombination hub

                                    long long node_start = -1, node_end = -1;
                                    if (row < (int)node_coords.size() &&
                                        col < (int)node_coords[row].seg.size()) {
                                        node_start = node_coords[row].seg[col].first;
                                        node_end   = node_coords[row].seg[col].second;
                                    }

                                    if (row != seg_row) {
                                        flush_segment();
                                        seg_row = row;
                                        seg_start = node_start;
                                        seg_cols = 0;
                                    }
                                    seg_end = node_end;
                                    seg_cols++;
                                }
                                flush_segment();

                                rout.close();
                                std::cout << "Saved haplotype-region breakdown to: " << region_filename << "\n";
                            } else {
                                std::cerr << "Error writing to file: " << region_filename << "\n";
                            }
                        }
                        // === END NEW ===

                        return n_recomb;
                    } else {
                        std::cout << "\nWarning: Path did not reach Sink. File not written.\n";
                        return -1;
                    }
                };

                std::filesystem::create_directories(out_dir);
                std::string sol_tag = (nsol > 1 && s > 0) ? "_sol" + std::to_string(s) : "";
                std::string file_path_1 = out_dir + "/" + sample_name + sol_tag + "_1_" + rel_weights + "_" + cov_label + ".fa";
                std::string file_path_2 = out_dir + "/" + sample_name + sol_tag + "_2_" + rel_weights + "_" + cov_label + ".fa";

                int n_recomb_1 = process_and_save_path(next_node_1, "x1", file_path_1);
                int n_recomb_2 = process_and_save_path(next_node_2, "x2", file_path_2);

                double novel_cost = 0.0;
                for (size_t e = 0; e < x_1_vars.size(); ++e) {
                    auto [u,v] = edge_list[e];
                    if (x_1_vars[e].get(GRB_DoubleAttr_Xn) > 0.5) novel_cost += c2 * novel_kmer_count[v];
                    if (x_2_vars[e].get(GRB_DoubleAttr_Xn) > 0.5) novel_cost += c2 * novel_kmer_count[v];
                }

                double missing_cost = 0.0;
                for (size_t i = 0; i < y_1_vars.size(); ++i) {
                    bool y1 = y_1_vars[i].get(GRB_DoubleAttr_Xn) > 0.5;
                    bool y2 = y_2_vars[i].get(GRB_DoubleAttr_Xn) > 0.5;
                    if (!y1 && !y2) {
                        missing_cost += c1;
                    }
                }

                std::cout << "\n-----------------------------------------\n";
                std::cout << "Recombinations (not in objective): " << n_recomb_1 << " + " << n_recomb_2 << "\n";
                std::cout << "Objective breakdown:\n";
                std::cout << "Novel kmer penalty: "      << novel_cost  << "\n";
                std::cout << "Missing kmer penalty: "    << missing_cost<< "\n";
                std::cout << "Total objective: "        << obj_val     << "\n";
            }
        }

    } catch(GRBException& e){
        std::cerr << "Gurobi Error: " << e.getMessage() << std::endl;
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "\nElapsed time: " << elapsed.count() << " seconds\n";

    return 0;
}