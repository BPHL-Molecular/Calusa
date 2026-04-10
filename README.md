<img width="2156" height="1365" alt="screenshot" src="https://github.com/user-attachments/assets/0021caff-62ef-4ad5-adf8-178f37b5bf72" />
# Calusa

**Calusa** is an interactive tool for visualizing pathogen transmission networks from genomic sequence data. It computes pairwise Hamming distances from FASTA-formatted haplotype sequences, identifies transmission clusters, and renders a force-directed network graph in the browser using D3.js.

Named after [Calusa Beach](https://www.floridastateparks.org/BahiaHonda) on Bahia Honda Key, Florida — one of the most iconic beaches in the Florida Keys.

Developed for genomic surveillance of bloodborne and foodborne pathogens (HCV, HAV, HBV) in collaboration with CDC and the Florida Department of Health Bureau of Public Health Laboratories ([BPHL-Molecular](https://github.com/BPHL-Molecular)).

---

<img width="2156" height="1365" alt="screenshot" src="https://github.com/user-attachments/assets/1b578bed-92f2-4915-a042-20eb75e1aca4" />


*Interactive D3.js force-directed network showing 3 transmission clusters (10 samples, 12 links, threshold = 0.0) with cluster-based coloring, sample labels, and summary statistics.*

---

## Overview

Calusa consists of two components:

**`calusa.py`** — A Python pipeline that reads multi-sequence FASTA files, calculates minimum pairwise Hamming distances between samples, identifies transmission clusters via depth-first search, and exports results as CSV and JSON.

**`calusa.html`** — A standalone, single-file HTML/JavaScript application that loads the exported JSON and renders an interactive D3.js force-directed network graph with zoom, pan, cluster coloring, tooltips, label toggling, and SVG export.

## Features

- Pairwise Hamming distance calculation with ambiguous-base handling (only A/T/C/G positions scored)
- Configurable distance threshold (default 0.037 for HCV; set 0 for HAV-VPB)
- Transmission cluster identification using connected-component DFS
- Interactive D3.js force-directed graph with cluster-based coloring
- Node sizing scaled by haplotype sequence count
- Distance-weighted link styling (thicker = closer)
- Hover tooltips showing sample ID, cluster, sequence count, and link distances
- Zoom, pan, label toggle, and SVG export controls
- Fully client-side visualization — no server required

## Requirements

**Python pipeline:**

- Python 3.8+
- pandas

Install dependencies:

```bash
pip install pandas
```

**Visualization:**

- Any modern web browser (Chrome, Firefox, Edge, Safari)
- D3.js v7 is loaded from CDN automatically

## Quick Start

### 1. Generate Network Data

Process a FASTA file with the default HCV threshold (0.037):

```bash
python calusa.py --input sequences.fasta
```

Specify a custom threshold and output directory:

```bash
python calusa.py --input sequences.fasta --threshold 0.05 --output ./results
```

For HAV with a zero-distance threshold:

```bash
python calusa.py --input sequences.fasta --threshold 0 --output ./results
```

Create a sample FASTA file for testing:

```bash
python calusa.py --create-sample
```

### 2. Visualize the Network

1. Open `calusa.html` in a web browser.
2. Click **Upload JSON** and select the generated `network.json`.
3. Interact with the network: zoom, pan, hover for details, toggle labels, and export as SVG.

## Input Format

The pipeline expects a standard multi-sequence FASTA file. Each header should follow the pattern `>SampleID_seqN` where the trailing `_seqN` (or `_hapN`) segment identifies individual haplotype sequences within a sample. All sequences sharing the same sample ID prefix are grouped together.

```
>Sample_001_seq_1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>Sample_001_seq_2
ATCGATCGATCGAACGATCGATCGATCGATCGATCGATCG
>Sample_002_seq_1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>Sample_003_seq_1
GGCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
```

## Output Files

The pipeline produces three output files:

| File | Description |
|------|-------------|
| `threshold_links.csv` | All sample pairs with minimum Hamming distance ≤ threshold, including cluster assignment |
| `cluster_summary.csv` | Per-sample summary with cluster ID and haplotype sequence count |
| `network.json` | JSON file for the D3.js visualization containing nodes, links, and metadata |

### JSON Schema

```json
{
  "metadata": {
    "threshold": 0.037,
    "total_samples": 100,
    "total_links": 250,
    "total_clusters": 8
  },
  "nodes": [
    { "id": "Sample_001", "num_sequences": 3, "cluster": 0 }
  ],
  "links": [
    { "source": "Sample_001", "target": "Sample_002", "distance": 0.012,
      "num_seqs_source": 3, "num_seqs_target": 2 }
  ]
}
```

## Distance Calculation

Hamming distance is computed as the proportion of differing nucleotides at valid (A/T/C/G) positions between two sequences. Gaps, N's, and other ambiguous bases are excluded from both the numerator and denominator. When comparing two samples with multiple haplotype sequences each, the **minimum** distance across all sequence pairs is used as the inter-sample distance.

## Command-Line Reference

```
usage: calusa.py [-h] [--input INPUT] [--threshold THRESHOLD]
                 [--output OUTPUT] [--create-sample]

optional arguments:
  -i, --input         Input FASTA file
  -t, --threshold     Transmission distance threshold (default: 0.037)
  -o, --output        Output directory (default: current directory)
  --create-sample     Create a sample FASTA file and exit
```

## Recommended Thresholds

| Pathogen | Threshold | Reference |
|----------|-----------|-----------|
| HCV      | 0.037     | Default   |
| HAV-VPB  | 0.000     | Exact match |

## License

This project is part of the BPHL-Molecular bioinformatics pipeline collection for public health genomic surveillance.
