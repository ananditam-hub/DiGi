import os
import yaml
import argparse

# Usage: python graph_builder.py --config graph.yml [OPTIONAL] --sample <sample_name>
#        [OPTIONAL] --samples_included <path/to/samples_included.txt>
# sample_name is especially crucial if the sample is a part of the reference graph
# samples_included.txt lists (one per line) the ONLY samples that should be
# included when building the graph

def load_config(config_path):
    with open(config_path) as f:
        config = yaml.safe_load(f)

    fmt = {
        "base_dir": config["paths"]["base_dir"]
    }

    for key, value in config["graph"].items():
        if isinstance(value, str):
            config["graph"][key] = value.format(**fmt)

    return config

def load_single_fasta(path):
    """Load the first (or only) sequence from a single FASTA file."""
    name = None
    buf = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name:
                    break  # only read the first sequence
                name = line[1:].split()[0]
            elif line:
                buf.append(line)
    return name, "".join(buf)

def write_node(out, node_start_id, written, seq, s, e):
    if s is None or e is None:
        return 0
    
    seq_len = len(seq)
    s = max(s, 0)
    e = min(e, seq_len) # Standard exclusive end
    
    # Skip if boundaries overlap or sequence is empty
    if e <= s:
        return 0
        
    node_id = node_start_id + written + 1
    out.write(f"{node_id}\t{seq[s:e]}\n")
    print(f">> Node {node_id}: start={s}, end={e}, length={e-s}")
    return 1

def extract_nodes_generic(sequence, bed_path, out, node_start):
    seq_len = len(sequence)
    written = 0
    start = 0
    
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip():
                continue
                
            parts = line.strip().split("\t")
            if len(parts) > 6:
                boundary = int(parts[6])  # col7 (6th index) acts as the continuous breakpoint/recombination hotspot
                
                # Create node from current position up to the boundary (seq[s:e])
                written += write_node(out, node_start, written, sequence, start, boundary)
                
                # The next node will automatically start exactly at this coordinate (seq[e:e'])
                start = boundary
                
    # Flush the final remainder of the sequence into the last node
    if start < seq_len:
        written += write_node(out, node_start, written, sequence, start, seq_len)

    return written

def filter_headers(headers, exclude_sample=None):
    if not exclude_sample:
        return headers

    prefix = f"{exclude_sample}"
    return [h for h in headers if not h.startswith(prefix)]

def load_included_samples(path):
    """Read a text file of sample names (one per line) to include.
    Blank lines and lines starting with '#' are ignored."""
    included = set()
    with open(path) as f:
        for line in f:
            name = line.strip()
            if name and not name.startswith("#"):
                included.add(name)
    return included

if __name__ == "__main__":

    # load YAML config
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, help="Path to YAML config file")
    parser.add_argument("--sample", default=None, help="Sample to be inferred")

    parser.add_argument("--samples_included", default=None,
                         help="Optional path to a text file listing sample names "
                              "(one per line) to include; graph is built only for these")
    args = parser.parse_args()

    config = load_config(args.config)
    exclude_sample = args.sample

    # paths from config
    fasta_folder = config["graph"]["haplotype_fasta_dir"]
    out_path = config["graph"]["graph_file"]
    bed_folder = config["graph"]["rh_bed_file"]

# discover haplotypes from FASTA files
    headers = sorted(
        os.path.splitext(f)[0]
        for f in os.listdir(fasta_folder)
        if f.endswith((".fa", ".fasta"))
    )
    # apply exclusions
    headers = filter_headers(headers, exclude_sample=exclude_sample)

    if args.samples_included:
        included = load_included_samples(args.samples_included)
        before = len(headers)
        headers = [h for h in headers if h in included]
        print(f"📋 samples_included filter: kept {len(headers)}/{before} haplotypes "
              f"based on {args.samples_included}")

        missing = included - set(headers)
        if missing:
            print(f"!! {len(missing)} sample(s) from samples_included not found "
                  f"in {fasta_folder}: {sorted(missing)}")

    print(f"Using {len(headers)} haplotypes after filtering")
    total = 0

    with open(out_path, "w") as out:
        for h in headers:
            fasta_file = os.path.join(fasta_folder, f"{h}.fa")
            bed_file = os.path.join(bed_folder, f"{h}/{h}_hotspots_lifted_final.bed") # *** extract bed files in one folder

            if not os.path.exists(fasta_file):
                print(f"!! Missing FASTA file: {fasta_file}")
                continue

            print(f"Loading {fasta_file}…")
            seq_name, seq = load_single_fasta(fasta_file)

            if not seq:
                print(f"!! Empty or unreadable FASTA: {fasta_file}")
                continue

            print(f"Processing {h} (FASTA header: {seq_name})…")

            try:
                # Note: removed the unused `is_grch38` argument from this call
                w = extract_nodes_generic(seq, bed_file, out, total)
                total += w
                print(f"   → wrote {w} nodes (running total={total})")
            except FileNotFoundError:
                print(f"Missing BED: {bed_file}")
            except Exception as e:
                print(f"Error in {h}: {e}")

    print(f"\n🏁 Done. Total nodes written: {total}")