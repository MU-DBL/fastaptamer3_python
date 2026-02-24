# Fastaptamer3 Backend API Documentation

**Base URL:** `/api/v1`
**Framework:** FastAPI
**Version:** 1.0.0

All endpoints accept and return JSON unless otherwise noted.

---

## Table of Contents

1. [File Handler](#1-file-handler)
2. [Count](#2-count)
3. [Recount](#3-recount)
4. [Preprocess](#4-preprocess)
5. [Progress (SSE)](#5-progress-sse)
6. [Cluster (LED)](#6-cluster-led)
7. [Cluster Diversity](#7-cluster-diversity)
8. [Cluster MSA](#8-cluster-msa)
9. [Cluster PHMM](#9-cluster-phmm)
10. [Cluster List](#10-cluster-list)
11. [Recluster](#11-recluster)
12. [Position Enrichment](#12-position-enrichment)
13. [Motif Search](#13-motif-search)
14. [Motif Discovery](#14-motif-discovery)
15. [Motif Tracker](#15-motif-tracker)
16. [Sequence Enrichment](#16-sequence-enrichment)
17. [Translate](#17-translate)
18. [Sequence Distance](#18-sequence-distance)
19. [Mutation Network](#19-mutation-network)
20. [Data Merge](#20-data-merge)
21. [Differential Expression](#21-differential-expression)

---

## 1. File Handler

Manages file upload, download, listing, and deletion in the server's upload directory.

### `POST /upload`

Upload a file to the server.

**Request:** `multipart/form-data`

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `file` | `UploadFile` | Yes | The file to upload |

**Response:**
```json
{
  "original_filename": "sample.fastq",
  "saved_filename": "sample.fastq",
  "message": "File uploaded successfully"
}
```

---

### `GET /download/{filename}`

Download a file from the server by filename.

**Path Parameter:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `filename` | `string` | Name of the file to download |

**Response:** The file as a binary stream with appropriate `Content-Disposition` and `Content-Type` headers.

**Errors:** `404` if file not found.

---

### `GET /files/`

List all files currently stored on the server.

**Response:**
```json
{
  "files": ["sample.fastq", "sample_count.csv"],
  "count": 2
}
```

---

### `DELETE /delete/{filename}`

Delete a file from the server.

**Path Parameter:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `filename` | `string` | Name of the file to delete |

**Response:**
```json
{
  "message": "File 'sample.fastq' deleted successfully"
}
```

**Errors:** `404` if file not found, `500` on deletion failure.

---

## 2. Count

### `POST /count`

Count sequences from a FASTQ or FASTA input file. Produces a ranked sequence table with read counts and RPU (Reads Per Unit) values.

**Request Body:**

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `input_path` | `string` | `""` | Filename of the input file (relative to upload directory) |
| `reverseComplement` | `boolean` | `false` | If `true`, also consider reverse complement of sequences |
| `scaling_factor` | `float` | `1e6` | Scaling factor used to compute RPU |
| `output_format` | `string` | `""` | Output format (`csv`, `tsv`, `fasta`) |

**Response:**
```json
{
  "status": "ok",
  "result": "sample_count.csv"
}
```

The output file contains ranked unique sequences with their read counts and RPU values.

---

## 3. Recount

### `POST /recount`

Combine two existing count files and produce a recombined count file. Used to merge sequence data from two experimental runs/replicates.

**Request Body:**

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `input_path_1` | `string` | `""` | First count file (relative to upload directory) |
| `input_path_2` | `string` | `""` | Second count file (relative to upload directory) |
| `scaling_factor` | `float` | `1e6` | Scaling factor for RPU calculation |
| `output_format` | `string` | `""` | Output format (`csv`, `tsv`, `fasta`) |

**Response:**
```json
{
  "status": "ok",
  "result": "combine_file1_file2.csv"
}
```

---

## 4. Preprocess

### `POST /preprocess`

Preprocess a raw FASTQ sequencing file. This is an **async** operation that runs in the background. The endpoint immediately returns a `job_id` that can be polled via the [Progress](#5-progress-sse) endpoint.

Operations performed:
- Trim 5' and 3' constant regions (primers)
- Filter sequences by length range
- Filter sequences by maximum error rate

**Request Body:**

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `input_path` | `string` | `null` | Input FASTQ filename |
| `const5p` | `string` | `""` | 5' constant/primer sequence to trim |
| `const3p` | `string` | `""` | 3' constant/primer sequence to trim |
| `min_length` | `integer` | `0` | Minimum sequence length after trimming |
| `max_length` | `integer` | `100` | Maximum sequence length after trimming |
| `max_error` | `float` | `0.005` | Maximum allowed error rate for primer matching |
| `output_format` | `string` | `"fasta"` | Output format (`fasta`, `fastq`) |

**Response:**
```json
{
  "status": "ok",
  "result": "abc123-job-id",
  "output_path": "sample_preprocess.fasta"
}
```

Use `result` (the `job_id`) with the `/progress/{job_id}` endpoint to track progress.

---

## 5. Progress (SSE)

### `GET /progress/{job_id}`

Stream real-time progress updates for a background job using **Server-Sent Events (SSE)**.

**Path Parameter:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `job_id` | `string` | The job ID returned by a background operation (e.g., `/preprocess`) |

**Response:** `text/event-stream` — streams JSON events:

```
data: {"stage": "connected", "message": "Processing large files may take some time..."}

data: {"stage": "processing", "message": "Reading file...", "progress": 20, "timestamp": 1700000000.0, "data": null}

data: {"stage": "complete", "message": "Done", "progress": 100, "timestamp": 1700000001.0, "data": {...}}
```

The stream closes after receiving a `stage` of `"complete"` or `"error"`. A `": keepalive"` line is sent every second if no data is available.

---

## 6. Cluster (LED)

### `POST /clusterled`

Cluster sequences using **Levenshtein Edit Distance (LED)**. Sequences within a defined edit distance threshold of a cluster seed are grouped together.

**Request Body:**

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `input_path` | `string` | `""` | Input file (count output, relative to upload dir) |
| `output_format` | `string` | `"fasta"` | Output format (`fasta`, `csv`) |
| `min_reads` | `integer` | `10` | Minimum read count for a sequence to become a cluster seed |
| `max_led` | `integer` | `7` | Maximum Levenshtein distance to assign a sequence to a cluster |
| `total_clusters` | `integer` | `30` | Maximum number of clusters to generate |
| `keep_nc` | `boolean` | `true` | Whether to keep non-clustered sequences in output |

**Response:**
```json
{
  "status": "ok",
  "result": "sample_cluster.fasta"
}
```

---

## 7. Cluster Diversity

### `POST /clusterdiversity`

Compute diversity statistics for clustered sequences.

**Request Body:**

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `input_path` | `string` | `""` | Clustered file from `/clusterled` |
| `output_format` | `string` | `"csv"` | Output format |

**Response:**
```json
{
  "status": "ok",
  "result": "sample_cluster_diversity.csv"
}
```

---

### `POST /cluster-kmer-analysis`

Perform **k-mer frequency analysis** on selected clusters with PCA or UMAP dimensionality reduction for visualization.

**Request Body:**

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `input_path` | `string` | _(required)_ | Clustered file |
| `selected_clusters` | `List[int]` | _(required)_ | List of cluster IDs to analyze |
| `k_size` | `integer` | `4` | K-mer size (2–8) |
| `method` | `string` | `"pca"` | Dimensionality reduction method: `"pca"` or `"umap"` |

**Response:**
```json
{
  "status": "ok",
  "data": { "coordinates": [...], "cluster_labels": [...], "colors": [...] }
}
```

---

## 8. Cluster MSA

### `POST /clustermsa`

Perform **Multiple Sequence Alignment (MSA)** on sequences from a single selected cluster using MUSCLE or Clustal.

**Request Body:**

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `input_path` | `string` | `""` | Clustered file |
| `output_format` | `string` | `"fasta"` | Output format |
| `seq_type` | `string` | `"dna"` | Sequence type: `"dna"` or `"protein"` |
| `cluster_selected` | `integer` | `1` | Cluster ID to align |

**Response:**
```json
{
  "status": "ok",
  "result": "sample_cluster1_msa.fasta",
  "num_sequences": 42
}
```

---

### `POST /cluster-msa-entropy`

Calculate **Shannon entropy** at each position of an MSA file. Higher entropy indicates greater variability at that position.

**Request Body:**

| Field | Type | Description |
|-------|------|-------------|
| `input_path` | `string` | MSA output file from `/clustermsa` |

**Response:**
```json
{
  "status": "ok",
  "positions": [1, 2, 3, ...],
  "entropy_values": [0.12, 0.85, 0.03, ...],
  "total_positions": 40
}
```

---

### `POST /cluster-msa-mutinfo`

Calculate **pairwise mutual information** between all positions of an MSA. Reveals co-varying positions that may indicate structural or functional constraints.

**Request Body:**

| Field | Type | Description |
|-------|------|-------------|
| `input_path` | `string` | MSA output file from `/clustermsa` |

**Response:**
```json
{
  "status": "ok",
  "mutual_info_matrix": [[...], [...]],
  "positions": [1, 2, 3, ...],
  "total_positions": 40
}
```

---

## 9. Cluster PHMM

### `POST /cluster-phmm-simulate`

Build a **Profile Hidden Markov Model (PHMM)** from an MSA and simulate new sequences. Useful for generating synthetic aptamer sequences that follow the same statistical profile as a cluster.

**Request Body:**

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `input_path` | `string` | `""` | MSA file from `/clustermsa` |
| `num_sequences` | `integer` | `10` | Number of sequences to simulate |
| `sequence_length` | `integer` | `70` | Target length of simulated sequences |
| `output_format_phmm` | `string` | `"txt"` | Format for saving the PHMM model |
| `output_format_simulation` | `string` | `"fasta"` | Format for saving simulated sequences |
| `pseudocount_method` | `string` | `"laplace"` | Pseudocount strategy: `"laplace"` or `"background"` |

**Response:**
```json
{
  "status": "ok",
  "simulated_sequences": ["ACGT...", "TGCA..."],
  "num_unique_sequences": 10,
  "output_file_phmm": "sample_phmm_model.txt",
  "output_file_simulation": "sample_phmm_simulated.fasta",
  "phmm_stats": {
    "num_sequences": 42,
    "alignment_length": 40,
    "num_match_states": 38,
    "pseudocount_method": "laplace"
  }
}
```

---

## 10. Cluster List

### `POST /cluster-list`

Retrieve the list of unique cluster IDs present in a clustered file.

**Request Body:**

| Field | Type | Description |
|-------|------|-------------|
| `input_path` | `string` | Clustered file from `/clusterled` |

**Response:**
```json
{
  "clusters": [1, 2, 3, 4, 5]
}
```

---

## 11. Recluster

### `POST /recluster`

Merge clusters from two experimental populations based on **Levenshtein Edit Distance** between cluster seeds. Sequences from both populations are combined into a unified cluster set and enrichment scores are calculated.

**Request Body:**

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `fadf1_cluster_path` | `string` | `""` | Clustered file for Population 1 |
| `fadf2_cluster_path` | `string` | `""` | Clustered file for Population 2 |
| `led_threshold` | `integer` | `7` | Maximum LED to consider two cluster seeds as the same cluster |
| `output_format` | `string` | `"csv"` | Output format |

**Response:**
```json
{
  "status": "ok",
  "result": "reclustered_led_7.csv",
  "num_sequences": 500,
  "num_clusters": 30
}
```

---

### `POST /recluster-led-matrix`

Compute the full **LED matrix** between cluster seeds of two populations. Returns data formatted for heatmap visualization.

**Request Body:**

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `fadf1_cluster_path` | `string` | `""` | Clustered file for Population 1 |
| `fadf2_cluster_path` | `string` | `""` | Clustered file for Population 2 |
| `led_threshold` | `integer` | `null` | Optional threshold for highlighting in the frontend |
| `use_parallel` | `boolean` | `true` | Whether to use parallel computation |
| `n_jobs` | `integer` | `null` | Number of parallel jobs (defaults to CPU count) |

**Response:**
```json
{
  "status": "ok",
  "led_matrix": [[0, 3, 7], [2, 0, 5]],
  "p1_cluster_ids": [1, 2],
  "p2_cluster_ids": [1, 2, 3],
  "p1_seeds": ["ACGT...", "TGCA..."],
  "p2_seeds": ["ACGT...", "TGCA...", "GGCC..."],
  "threshold": 7,
  "matrix_shape": {"rows": 2, "cols": 3},
  "statistics": {"min": 0, "max": 7, "mean": 3.1, "median": 3.0, "std": 2.1}
}
```

---

## 12. Position Enrichment

### `POST /position-enrichment`

Calculate **positional enrichment** for sequences in a specific cluster from a reclustered file. Each residue (nucleotide or amino acid) at each alignment position is scored by its average enrichment value.

**Request Body:**

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `fadf_recluster_path` | `string` | `""` | Reclustered file from `/recluster` |
| `cluster_selection` | `integer` | `1` | Cluster ID to analyze |
| `seq_type` | `string` | `"dna"` | Sequence type: `"dna"` or `"protein"` |
| `output_format` | `string` | `"csv"` | Output format |
| `max_sequences` | `integer` | `500` | Max sequences for MSA (prevents timeout) |

**Response:**
```json
{
  "status": "ok",
  "result": "pos_enrich_cluster_1.csv",
  "enrichment_matrix": [[1.2, 0.8, ...], [0.5, 2.1, ...]],
  "residues": ["A", "C", "G", "T", "-"],
  "positions": [1, 2, 3, ...],
  "avg_enrichment": [1.05, 1.31, ...]
}
```

---

## 13. Motif Search

### `POST /motif-search`

Search sequences for one or more user-defined motifs. Returns sequences that **contain** the motifs.

**Request Body:**

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `input_path` | `string` | `""` | Input FASTA/CSV file |
| `motif` | `string` | `""` | Comma-separated list of motif patterns |
| `highlight` | `boolean` | `false` | Whether to highlight matched motifs in output |
| `partial` | `boolean` | `false` | If `true`, return sequences matching ANY motif (OR); if `false`, ALL motifs (AND) |
| `motif_type` | `string` | `"Nucleotide"` | Motif type: `"Nucleotide"`, `"AminoAcid"`, or `"String"` |
| `output_format` | `string` | `"fasta"` | Output format |

**Response:**
```json
{
  "status": "ok",
  "result": "sample_motif_search.fasta"
}
```

---

### `POST /motif-omit`

Omit sequences that contain one or more user-defined motifs. The inverse of `/motif-search`.

**Request Body:**

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `input_path` | `string` | `""` | Input FASTA/CSV file |
| `motif` | `string` | `""` | Comma-separated list of motif patterns |
| `partial` | `boolean` | `false` | If `true`, omit sequences matching ANY motif (OR); if `false`, ALL motifs (AND) |
| `motif_type` | `string` | `"Nucleotide"` | Motif type |
| `output_format` | `string` | `"fasta"` | Output format |

**Response:**
```json
{
  "status": "ok",
  "result": "sample_motif_omit.fasta"
}
```

---

## 14. Motif Discovery

### `POST /motif-discovery`

Discover statistically **enriched sub-sequence motifs** from a set of sequences. For each candidate motif of a given length range, computes observed vs. expected frequency, Z-score, and p-value using a background nucleotide model.

**Request Body:**

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `input_path` | `string` | _(required)_ | Input file with sequences |
| `min_reads` | `integer` | `10` | Minimum read count for a sequence to be included |
| `length_range` | `[int, int]` | `[5, 10]` | `[min_motif_length, max_motif_length]` |
| `output_format` | `string` | `"csv"` | Output format |

**Response:**
```json
{
  "status": "ok",
  "result": "sample_motif_discovery.csv"
}
```

The output file contains columns: `Motif`, `P_Value`, `ZZ_Score`, `MotifLength`, `Rank`, `SeqCount`.

---

## 15. Motif Tracker

### `POST /motif-tracker`

Track the occurrence of one or more **motif patterns** across multiple populations/rounds. Reports total read counts and percentage for each motif in each population, plus enrichment between consecutive populations.

**Request Body:**

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `input_paths` | `List[string]` | _(required)_ | Ordered list of FASTA/CSV files (one per population) |
| `population_names` | `List[string]` | _(required)_ | Display name for each population (must match length of `input_paths`) |
| `query_list` | `List[string]` | _(required)_ | Motif patterns to search for (IUPAC codes supported for Nucleotide type) |
| `query_aliases` | `List[string]` | `null` | Optional display aliases for each motif |
| `motif_type` | `string` | `"Nucleotide"` | Motif type: `"Nucleotide"`, `"AminoAcid"`, or `"String"` |
| `output_format` | `string` | `"csv"` | Output format |

**Response:**
```json
{
  "status": "ok",
  "result": "sample_motif_tracker.csv",
  "enrichment": "sample_motif_tracker_enrichment.csv"
}
```

---

### `POST /sequence-tracker`

Track exact **full-sequence matches** across multiple populations. Unlike `/motif-tracker`, this performs exact string matching rather than regex/motif matching.

**Request Body:** Same structure as `/motif-tracker`. The `query_list` field contains exact sequences to track instead of motif patterns.

**Response:**
```json
{
  "status": "ok",
  "result": "sample_sequence_tracker.csv",
  "enrichment": "sample_sequence_tracker_enrichment.csv"
}
```

---

## 16. Sequence Enrichment

### `POST /sequence-enrich`

Calculate **sequence-level enrichment** between two populations by comparing RPU values. Produces a merged table with Enrichment, log2(Enrichment), and MA-plot values (R and A).

**Request Body:**

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `fadf1_cluster_path` | `string` | `""` | Population 1 file (baseline) |
| `fadf2_cluster_path` | `string` | `""` | Population 2 file (comparison) |
| `keep_na` | `boolean` | `false` | If `true`, include sequences present in only one population (outer join); if `false`, only shared sequences (inner join) |
| `output_format` | `string` | `"csv"` | Output format |

**Response:**
```json
{
  "status": "ok",
  "result": "enriched_inner.csv",
  "num_sequences": 312,
  "enrichment_stats": {
    "mean_enrichment": 1.45,
    "median_enrichment": 1.12,
    "mean_log2E": 0.53,
    "median_log2E": 0.16
  }
}
```

---

## 17. Translate

### `POST /translate`

Translate DNA sequences to **amino acid sequences** using a selected genetic code. Optionally merges non-unique amino acid sequences (convergence).

**Request Body:**

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `input_path` | `string` | _(required)_ | Input DNA sequence file |
| `orf` | `integer` | `1` | Open reading frame offset (1 = start from nucleotide 1) |
| `converge` | `boolean` | `true` | If `true`, merge sequences that translate to the same amino acid sequence and sum their reads |
| `input_changes` | `List[CodonChange]` | `null` | Custom codon overrides: `[{"Codon": "TGA", "Translation": "W"}]` |
| `translate_selection` | `string` | `"Standard"` | Genetic code to use (e.g., `"Standard"`, `"Vertebrate Mitochondrial"`) |
| `output_format` | `string` | `"csv"` | Output format |

**Response:**
```json
{
  "status": "ok",
  "result": "sample_translated.csv"
}
```

The output file adds a `Length` column and optionally a `Unique_Nts` column (when `converge=true`) indicating how many distinct DNA sequences encode each amino acid sequence.

---

## 18. Sequence Distance

### `POST /sequence-distance`

Calculate the **Levenshtein edit distance** between a query sequence and all sequences in a file. Results are sorted by distance (ascending).

**Request Body:**

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `input_path` | `string` | _(required)_ | Input FASTA/CSV file |
| `query_sequence` | `string` | _(required)_ | The reference sequence to compare against |
| `output_format` | `string` | `"csv"` | Output format |

**Response:**
```json
{
  "status": "ok",
  "result": "sample_distance.csv"
}
```

The output file includes all original columns plus a `Distance` column.

---

## 19. Mutation Network

### `POST /mutation-network`

Find the **shortest mutation path** between two sequences through a network of sequences. Builds a graph where sequences within `max_cost` edit distance are connected, then computes the shortest path (Dijkstra's algorithm) from `start_node` to `end_node`.

**Request Body:**

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `input_path` | `string` | _(required)_ | Input FASTA/CSV file containing the sequence network |
| `start_node` | `string` | _(required)_ | Starting sequence |
| `end_node` | `string` | _(required)_ | Target sequence |
| `max_cost` | `integer` | `1` | Maximum edit distance to draw an edge between two sequences |
| `output_format` | `string` | `"csv"` | Output format (`csv` or `tsv`) |

**Response:**
```json
{
  "status": "ok",
  "result": "sample_mutation_network.csv",
  "path_length": 4
}
```

The output file contains columns: `From_Sequence`, `To_Sequence`, `Transition_Cost`.

---

## 20. Data Merge

### `POST /data-merge`

Merge multiple sequence files (2–26 files) on the `Sequences` column. Each population's columns are renamed with a suffix (`.a`, `.b`, `.c`, ...).

**Request Body:**

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `input_paths` | `List[string]` | _(required)_ | List of 2–26 files to merge |
| `merge_type` | `string` | `"outer"` | SQL join type: `"outer"`, `"inner"`, `"left"`, `"right"` |
| `output_format` | `string` | `"csv"` | Output format |

**Response:**
```json
{
  "status": "ok",
  "result": "sample_merged.csv"
}
```

---

### `POST /sequence-persistence`

Analyze **sequence persistence** from a merged file. Counts how many populations each sequence appears in and returns a frequency distribution.

**Request Body:**

| Field | Type | Description |
|-------|------|-------------|
| `merged_file_path` | `string` | Merged file from `/data-merge` |

**Response:**
```json
{
  "status": "ok",
  "data": [
    {"freq": 1, "seqCount": 150},
    {"freq": 2, "seqCount": 80},
    {"freq": 3, "seqCount": 20}
  ]
}
```

`freq` = number of populations a sequence appears in; `seqCount` = number of sequences with that persistence.

---

### `POST /upset-data`

Compute **UpSet intersection data** from a merged file. Calculates the size of every exclusive set intersection for visualization in an UpSet plot.

**Request Body:**

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `merged_file_path` | `string` | _(required)_ | Merged file from `/data-merge` |
| `fasta_names` | `List[string]` | `null` | Optional custom display names for each population |

**Response:**
```json
{
  "status": "ok",
  "sets": ["Population_a", "Population_b", "Population_c"],
  "set_sizes": {"Population_a": 200, "Population_b": 180, "Population_c": 160},
  "intersections": [
    {"sets": ["Population_a"], "size": 50, "sequences": ["ACGT...", ...]},
    {"sets": ["Population_a", "Population_b"], "size": 30, "sequences": [...]}
  ],
  "total_unique_sequences": 350
}
```

Intersections are exclusive (sequences in the listed sets only) and sorted by size descending. Up to 100 example sequences are included per intersection.

---

## 21. Differential Expression

### `POST /differential-expression`

Perform **differential expression analysis** between two experimental conditions (each with one or more replicate files). Uses CPM normalization and either Welch's t-test (≥3 replicates) or Mann-Whitney U test (fewer replicates), with Benjamini-Hochberg FDR correction.

**Request Body:**

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `cond1_paths` | `List[string]` | _(required)_ | List of FASTA files for Condition 1 (baseline) |
| `cond2_paths` | `List[string]` | _(required)_ | List of FASTA files for Condition 2 (treatment) |
| `p_cutoff` | `float` | `0.1` | FDR significance threshold (between 0 and 1) |
| `output_format` | `string` | `"csv"` | Output format |

Analysis pipeline:
1. Merge replicates within each condition (inner join on sequences present in all replicates)
2. Normalize by library size to CPM (Counts Per Million)
3. Filter low-count sequences (min 5 reads in at least 2 samples)
4. Calculate log fold change (logFC) and average logCPM
5. Apply statistical test and FDR correction

**Response:**
```json
{
  "status": "ok",
  "result": "sample_diff_expr.csv"
}
```

The output file contains: `Sequence`, `logFC`, `logCPM`, `PValue` (adjusted), `PClass` (`"Sig."` or `"Not Sig."`).

---

## Common Error Responses

| Status Code | Meaning |
|-------------|---------|
| `400` | Bad request — missing or invalid input parameters |
| `404` | File not found |
| `500` | Internal server error — processing failed |

All error responses follow the format:
```json
{
  "detail": "Descriptive error message"
}
```

---

## Root Endpoints

### `GET /`
Returns a welcome message.
```json
{"message": "Welcome to Fastaptamer3 Project"}
```

### `GET /health`
Health check endpoint.
```json
{"status": "healthy"}
```
