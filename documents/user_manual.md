# FASTAptameR-3 User Manual

---

## Table of Contents

- [3.1 FASTAptameR-Start](#31-fastapt ameR-start)
  - [3.1.1 Preprocess](#311-preprocess)
  - [3.1.2 Count](#312-count)
  - [3.1.3 Recount](#313-recount)
- [3.2 FASTAptameR-Cluster](#32-fastaptamer-cluster)
  - [3.2.1 Cluster](#321-cluster)
  - [3.2.2 Diversity](#322-diversity)
  - [3.2.3 MSA](#323-msa-multiple-sequence-alignment)
  - [3.2.4 PHMM](#324-phmm-profile-hidden-markov-model)
  - [3.2.5 Recluster](#325-recluster)
  - [3.2.6 Position Enrichment](#326-position-enrichment)
- [3.3 FASTAptameR-Motif](#33-fastaptamer-motif)
  - [3.3.1 Motif Discovery](#331-motif-discovery)
  - [3.3.2 Motif Search](#332-motif-search)
  - [3.3.3 Motif Omit](#333-motif-omit)
  - [3.3.4 Motif Tracker](#334-motif-tracker)
- [3.4 FASTAptameR-Translate](#34-fastaptamer-translate)
- [3.5 FASTAptameR-Sequence Enrichment](#35-fastaptamer-sequence-enrichment)
- [3.6 FASTAptameR-Data Merge](#36-fastaptamer-data-merge)
- [3.7 FASTAptameR-Differential Analysis](#37-fastaptamer-differential-analysis)
- [3.8 FASTAptameR-Distance](#38-fastaptamer-distance)
- [3.9 FASTAptameR-Mutation Network](#39-fastaptamer-mutation-network)
- [3.10 FASTAptameR-Pipeline](#310-fastaptamer-pipeline)

---

## 3.1 FASTAptameR-Start

The Start page is the entry point for processing raw sequencing data. It contains three tabs: **Preprocess**, **Count**, and **Recount**.

---

### 3.1.1 Preprocess

#### Description

Trims constant 5′ and 3′ flanking regions from sequences in a FASTQ or FASTA file, then filters sequences by length range and maximum average sequencing error probability. Outputs a cleaned FASTA file ready for counting. Processing progress is shown in real time via a log panel and progress bar on the right side of the screen.

#### Usage

1. Click **Browse...** to upload a FASTQ or FASTA file (`.fasta`, `.fastq`, `.fa`, `.fq`).
2. Enter the constant 5′ and 3′ primer/flanking sequences to be trimmed.
3. Set the sequence length range and maximum allowed error using the sliders.
4. Click **Start**. Progress logs appear on the right panel in real time.
5. Click **Download** to save the cleaned FASTA output.

#### Parameters

| Parameter | Description |
|---|---|
| Constant 5′ region | Sequence string to trim from the 5′ end of each read. Leave blank to skip 5′ trimming. |
| Constant 3′ region | Sequence string to trim from the 3′ end of each read. Leave blank to skip 3′ trimming. |
| Sequence length range (Min / Max) | After trimming, sequences outside this length window are discarded. Range: 1 to 300 nt. |
| Maximum allowed average error probability | Maximum mean error probability per sequence, derived from Phred quality scores. Sequences exceeding this threshold are discarded. Range: 0.001–1.0. Lower values retain only higher-quality reads. |

---

### 3.1.2 Count

#### Description

Counts and ranks unique sequences from a FASTA or FASTQ file. Each unique sequence is assigned a read count, a rank by abundance, and a Reads Per Unit (RPU) value normalized to a user-selected scale. Outputs a ranked sequence table that is the standard input for most downstream modules.

#### Usage

1. Click **Browse...** to upload a FASTA or FASTQ file.
2. If uploading a previously generated count result, click **Visualize Counting Results** to load it directly into the table and plots without re-running.
3. Set normalization scale, reverse complement option, and download format.
4. Click **Start**. Total and unique sequence counts are displayed after completion.
5. Click **Download** to save the output.

#### Parameters

| Parameter | Default | Description |
|---|---|---|
| Normalization value (RPU scale) | — | Denominator for RPU calculation. Options: 1, 10, 100, 1,000, 10,000, 1×10⁵, 1×10⁶. Select to match your downstream analysis convention. |
| Return reverse complement? | No | If Yes, the reverse complement of each sequence is reported instead of the original strand. |
| FASTA or CSV download? | FASTA | Output format. FASTA is required for most downstream modules. |

#### Output Table Columns

| Column | Description |
|---|---|
| Rank | Global rank of the sequence by read count (1 = most abundant) |
| Reads | Raw read count |
| RPU | Reads per unit, normalized to the selected scale |
| Sequence | Nucleotide sequence |

#### Plotting

Three plots are available after a successful run. All plots share two pre-filter sliders:

| Slider | Range | Description |
|---|---|---|
| Min. number of reads to plot | 0–1,000 (step 10) | Exclude sequences with fewer reads from the plot |
| Max. rank to plot | 10–1,000 (step 10) | Limit the plot to sequences up to this rank |

**Reads per Rank Plot**

Line chart of sequence abundance (y-axis) vs. rank (x-axis). Shows the abundance distribution of the population.

- Click **Reads per rank** to generate.
- Set **Adjust default reads-per-rank plot?** to **Yes** for customization.

Customization options: X-axis label, Y-axis label, plot title, line color.

**Sequence-Length Histogram**

Stacked bar chart showing the number of unique sequences (top) and total reads (bottom) per sequence length bin.

- Click **Sequence-length histogram** to generate.
- Set **Adjust default sequence-length histogram?** to **Yes** for customization.

Customization options: X-axis label, Y-axis 1 label (unique), Y-axis 2 label (total), plot title, bar fill colors for unique and total.

**Abundance Plot**

Bar chart grouping sequences by read-count bins, showing what fraction of the population falls in each abundance tier. A gradient color encodes the number of unique sequences per bin.

- Click **Abundance plot** to generate.
- Set **Adjust default abundance plot?** to **Yes** for customization.

Customization options: singleton bin toggle, comma-separated breakpoints, X-axis label, Y-axis label, plot title, gradient start/end colors.

---

### 3.1.3 Recount

#### Description

Integrates read counts from 2 to 5 counted FASTA files, consolidates unique sequences across all files, and recomputes RPU values normalized to a selected scale. Useful for pooling technical replicates or combining data from multiple sequencing runs into a single ranked table.

#### Usage

1. Click **Browse...** and select 2–5 FASTA files at once (or upload them one at a time).
2. After all files appear in the upload list, set normalization scale and download format.
3. Click **Start**. Total and unique sequence counts are shown after completion.
4. Click **Download** to save the recounted output.

> **Note:** Files must all be counted FASTA files (output from Count or Preprocess → Count).

#### Parameters

| Parameter | Default | Description |
|---|---|---|
| Normalization value (RPU scale) | — | Same scaling options as Count: 1, 10, 100, 1,000, 10,000, 1×10⁵, 1×10⁶. |
| FASTA or CSV download? | FASTA | Output format. |

#### Plotting

Same plots as Count (Reads per Rank and Sequence-Length Histogram) with identical customization options, available after processing completes.

---

## 3.2 FASTAptameR-Cluster

The Cluster page groups sequences into clusters and provides tools for characterizing cluster structure. It contains six tabs.

---

### 3.2.1 Cluster

#### Description

Groups sequences from a counted FASTA file into clusters using a greedy Levenshtein Edit Distance (LED) algorithm. The most abundant unclustered sequence seeds each new cluster; all remaining unclustered sequences within a defined edit distance are assigned to that cluster. Sequences not assigned are labeled NC (Not Clustered).

#### Usage

1. Click **Browse...** to upload a counted FASTA file.
2. Set parameters using sliders or number inputs.
3. Click **Start**. Results appear in the table on the right.
4. Click **Download** to save the output as FASTA or CSV.

> FASTA format is required for downstream Cluster modules.

#### Output Table Columns

| Column | Description |
|---|---|
| ID | Sequence identifier |
| Cluster | Assigned cluster number |
| Rank In Cluster | Rank within the cluster by read count |
| LED | Levenshtein edit distance from the cluster seed |
| Reads | Raw read count |
| Rank | Global rank across the population |
| RPU | Reads per unit |
| Sequence | Nucleotide or amino acid sequence |

#### Parameters

| Parameter | Default | Range | Description |
|---|---|---|---|
| Min. number of reads to cluster | 10 | 0–1,000 (step 5) | Sequences with fewer reads are excluded before clustering. Higher values shorten runtime. |
| Max. LED | 7 | 1–20 (step 1) | Maximum edit distance between a sequence and a cluster seed for cluster membership. Use FASTAptameR-Distance to guide this threshold. |
| Max. number of clusters to generate | 20 | 5–1,000 (step 5) | Caps the number of clusters produced. Remaining sequences are labeled NC. |
| Keep non-clustered sequences? | No | Yes / No | If Yes, NC-labeled sequences are retained in the output. |
| FASTA or CSV download? | FASTA | FASTA / CSV | Output format. |

---

### 3.2.2 Diversity

#### Description

Computes per-cluster diversity statistics from a clustered FASTA file: unique sequence count, total reads, total RPU, and average LED from the seed. Results are shown in a summary table and visualized through cluster metaplots and k-mer composition plots.

#### Usage

1. Upload a clustered FASTA or previously generated diversity CSV file.
2. To load an existing diversity result directly, click **Visualize Results**.
3. Click **Start** to compute statistics. Results appear in the table.
4. Click **Download** to save the diversity CSV.

#### Output Table Columns

| Column | Description |
|---|---|
| Cluster | Cluster number |
| Total Sequences | Number of unique sequences |
| Total Reads | Sum of all read counts |
| Total RPU | Sum of RPU values |
| Average LED | Mean edit distance from the seed |
| SeedID | ID of the cluster's seed sequence |

#### Plotting

**Cluster Metaplots**

Three stacked line charts per cluster: unique sequence count (top), total reads (middle), average LED (bottom).

- Click **Cluster Metaplots** after running analysis.
- Set **Adjust default plots of cluster metadata?** to **Yes** for customization.

| Customization Option | Description |
|---|---|
| X-axis label | Label for the cluster axis |
| Y-axis 1 label | Label for unique sequence count trace |
| Color for Y-axis 1 | Line color for unique sequence count |
| Y-axis 2 label | Label for total read count trace |
| Color for Y-axis 2 | Line color for total read count |
| Y-axis 3 label | Label for average LED trace |
| Color for Y-axis 3 | Line color for average LED |
| Plot title | Title above the chart |

**k-mer Plot (PCA / UMAP)**

Projects sequences from selected clusters into 2D space using k-mer composition features.

1. Select clusters from the **Available** chip list.
2. Choose **k-mer size**: 3, 4, or 5.
3. Choose **Plot type**: PCA or UMAP.
4. Click **k-mer plot**.

> Characters outside `[A, C, G, T, U]` are converted to `X` before k-mer computation.

Set **Adjust default k-mer plot?** to **Yes** for customization.

| Customization Option | Description |
|---|---|
| X-axis label | Label for the first reduced dimension |
| Y-axis label | Label for the second reduced dimension |
| Legend title | Label for the cluster color legend |
| Plot title | Title above the scatter plot |
| Color palette | Qualitative palette: Accent, Dark2, Paired, Pastel1, Pastel2, Set1, Set2 |

---

### 3.2.3 MSA (Multiple Sequence Alignment)

#### Description

Performs Multiple Sequence Alignment on all sequences within a single selected cluster from a clustered FASTA file. Outputs an aligned FASTA. Provides entropy and mutual information plots to characterize positional variability.

> MSA is computationally intensive. Runtime scales with cluster size.

#### Usage

1. Upload a clustered FASTA file. To load an existing MSA result directly, click **Visualize Results**.
2. Select the **cluster number** from the dropdown.
3. Choose **Type of sequences**: Nucleotide or AminoAcid.
4. Choose **download format**: FASTA or CSV.
5. Click **Start**. Aligned sequences appear in the table.
6. Click **Download** to save the output.

#### Parameters

| Parameter | Description |
|---|---|
| Cluster | Cluster number to align (populated from uploaded file) |
| Type of sequences | Nucleotide or AminoAcid — selects the appropriate substitution model |
| FASTA or CSV download? | Output format for the aligned sequences |

#### Plotting

**MSA Alignment Chart**

Stacked bar chart showing nucleotide/amino acid composition at each alignment position.

- Click **MSA alignment chart** (requires completed run with data in table).

**MSA Entropy Plot**

Bar chart of Shannon entropy (in nats) per alignment position.

- Click **MSA entropy** after a successful run.
- Set **Adjust default entropy plot?** to **Yes** for customization.

| Customization Option | Description |
|---|---|
| X-axis label | Label for alignment position axis |
| Y-axis label | Label for entropy axis (nats) |
| Legend title | Label for the legend |
| Plot title | Title above the chart |
| Bar outline color | Color of bar borders |
| Bar fill color | Color of bar interiors |

**MSA Mutual Information Plot**

Heatmap of pairwise mutual information between all alignment positions, revealing co-varying positions.

- Click **MSA mutual information** after a successful run.
- Set **Adjust default mutual information plot?** to **Yes** for customization.

| Customization Option | Description |
|---|---|
| X-axis label | Label for positions on x-axis |
| Y-axis label | Label for positions on y-axis |
| Legend title | Label for the MI color scale |
| Plot title | Title above the heatmap |
| Fill palette | Color palette for the heatmap |

---

### 3.2.4 PHMM (Profile Hidden Markov Model)

#### Description

Builds a Profile Hidden Markov Model (PHMM) from a cluster's MSA output, then simulates new sequences that statistically match the profile. Enables in-silico exploration of the cluster's sequence space. Outputs the PHMM model file and a FASTA of simulated sequences.

#### Usage

1. Upload a **cluster MSA FASTA file** (output from the MSA tab).
2. Set the number of sequences to simulate and the target sequence length.
3. Click **Start**.
4. Click **Download PHMM** to save the model file, or **Download sequences** to save the simulated FASTA.

#### Parameters

| Parameter | Range | Description |
|---|---|---|
| Simulate how many sequences? | 10–1,000 (step 1) | Number of new sequences to generate from the PHMM |
| Length of simulated sequence (including gaps) | 10–200 (step 1) | Target length of each simulated sequence, counting gap characters |

---

### 3.2.5 Recluster

#### Description

Merges cluster families from separately clustered FASTA populations by comparing seed sequences using pairwise Levenshtein Edit Distance. Cluster families whose seeds fall within the LED threshold are merged into a single super-cluster. Two modes are available:

- **2 Populations** — cross-compare two populations (e.g., pre- and post-selection). Outputs a cluster-level summary table and a sequence-level table with per-sequence enrichment.
- **Multi-round (3 files)** — track cluster family convergence across three consecutive selection rounds. Outputs a single merged table with abundance and enrichment columns for all three rounds.

> **Avg RPU** — mean RPU across all sequences assigned to the super-cluster in that population. Reflects overall cluster family abundance but is pulled down by low-count peripheral variants.
>
> **Seed RPU** — RPU of the cluster seed (the most abundant sequence in the original cluster). More robust than Avg RPU because it represents the dominant variant only, unaffected by rare sequences at the cluster periphery.
>
> Use Seed RPU to compare dominant aptamers between rounds. Use Avg RPU to capture the enrichment of the entire cluster family.

#### Usage — 2 Populations Mode

1. Select **2 Populations** mode.
2. Upload **Cluster data 1** (first population clustered FASTA or CSV).
3. Upload **Cluster data 2** (second population clustered FASTA or CSV).
4. Set the **Max. LED** threshold.
5. Click **Recluster**. The cluster-level summary table appears on the right; if sequence-level data is available, it appears below.
6. Click **Download Summary** to save the cluster-level CSV.
7. Click **Download Reclustered Sequences** to save the sequence-level CSV.

#### Usage — Multi-round Mode

1. Select **Multi-round (3 files)** mode.
2. For each round (Round 1, Round 2, Round 3), enter a **Round label** (e.g., LR8, LR10, LR15) and upload the corresponding clustered FASTA or CSV.
3. Set the **Max. LED** threshold.
4. Click **Recluster**. Results appear in the table.
5. Click **Download Summary** to save the output.

#### Parameters

| Parameter | Range | Description |
|---|---|---|
| Max. LED | 1–20 (step 1) | Maximum edit distance between seeds from different populations for them to be merged into the same super-cluster. |
| Round label (Multi-round only) | — | Display name for each round, used as column suffixes in the output (e.g., "LR8" → `AvgRPU.LR8`). |

#### Output Table Columns — 2 Populations

**Cluster-level summary:**

| Column | Description |
|---|---|
| Super-cluster | Merged cluster identifier |
| Seed | Seed sequence of the super-cluster |
| Size (Pop 1 / Pop 2) | Number of unique sequences in each population |
| Avg RPU (Pop 1 / Pop 2) | Mean RPU across all sequences in the super-cluster per population |
| Enrichment (Avg) | Avg RPU (Pop 2) / Avg RPU (Pop 1) |
| log2E (Avg) | log₂ of the average RPU enrichment |
| Seed RPU (Pop 1 / Pop 2) | RPU of the cluster seed in each population |
| Enrichment (Seed) | Seed RPU (Pop 2) / Seed RPU (Pop 1) |
| log2E (Seed) | log₂ of the seed RPU enrichment |
| Status | Whether the super-cluster was found in both populations or only one |

**Sequence-level data:**

| Column | Description |
|---|---|
| Sequence | Nucleotide or amino acid sequence |
| Super-cluster | Merged cluster assignment |
| Orig. Cluster (Pop 1 / Pop 2) | Original cluster number from each input file |
| Rank In Cluster | Rank within the super-cluster by read count |
| LED | Edit distance from the super-cluster seed |
| RPU (Pop 1 / Pop 2) | RPU of the sequence in each population |
| Enrichment | RPU (Pop 2) / RPU (Pop 1) per sequence |
| log2E | log₂ per-sequence enrichment |

#### Output Table Columns — Multi-round

Columns are named using the round labels you supplied (e.g., R1, R2, R3):

| Column | Description |
|---|---|
| Super-cluster | Merged cluster identifier |
| Seed | Seed sequence of the super-cluster |
| AvgRPU.R1 / R2 / R3 | Mean RPU across all super-cluster sequences per round |
| SeedRPU.R1 / R2 / R3 | RPU of the cluster seed per round |
| E.R1.R2 | AvgRPU.R2 / AvgRPU.R1 — average RPU enrichment from round 1 to round 2 |
| E.R2.R3 | AvgRPU.R3 / AvgRPU.R2 — average RPU enrichment from round 2 to round 3 |
| SeedE.R1.R2 | SeedRPU.R2 / SeedRPU.R1 — seed RPU enrichment from round 1 to round 2 |
| SeedE.R2.R3 | SeedRPU.R3 / SeedRPU.R2 — seed RPU enrichment from round 2 to round 3 |

#### Plotting

**2 Populations mode** — four plots, each with an optional customization panel.

**LED Heatmap**

Pairwise LED between all cluster seeds of population 1 (x-axis) vs. population 2 (y-axis). Click **LED Heatmap** (both files must be uploaded before running).

| Customization Option | Description |
|---|---|
| X-axis label | Population 1 cluster axis label |
| Y-axis label | Population 2 cluster axis label |
| Legend title | Label for the LED color scale |
| Plot title | Title above the heatmap |
| Color palette | Sequential palette: Magma, Inferno, Plasma, Viridis, Cividis, Hot, Jet, Turbo |

**Population Size Plot**

Side-by-side bar chart comparing unique sequence counts per super-cluster between the two populations. Click **Population Plot** (requires recluster results in the table).

| Customization Option | Description |
|---|---|
| X-axis label | Cluster axis label |
| Y-axis label | Unique sequence count label |
| Legend title | Population legend label |
| Plot title | Title above the chart |
| Population 1 / 2 bar color | Bar color per population |

**RPU Plot**

Compares mean Avg RPU per super-cluster between the two populations. Click **RPU Plot**.

Customization options: same as Population Size Plot.

**Enrichment Plot**

Bar chart of enrichment (Avg RPU Pop 2 / Pop 1) per super-cluster. Click **Enrichment Plot**.

| Customization Option | Description |
|---|---|
| X-axis label | Cluster axis label |
| Y-axis label | Enrichment label |
| Plot title | Title above the chart |
| Bar fill / outline color | Bar colors |

**Multi-round mode** — one plot.

**Trajectory Plot**

Line chart showing how Avg RPU evolves across the three rounds for each super-cluster. Click **Trajectory Plot** (requires recluster results in the table).

| Customization Option | Description |
|---|---|
| Y-axis label | Label for the RPU axis |
| Plot title | Title above the chart |

---

### 3.2.6 Position Enrichment

#### Description

Calculates per-position nucleotide or amino acid enrichment by aligning sequences within a cluster and computing the mean enrichment value at each alignment position. Reveals which positions are most conserved (high enrichment) or most variable. Accepts a Recluster CSV, a Sequence Enrichment CSV, or any CSV that contains an `Enrichment` column alongside a `Sequence` column.

> Position enrichment performs MSA internally and can be slow for large clusters. Wait on the page until processing completes.

#### Usage

1. Upload an **Input CSV file** — a Recluster CSV (from the Recluster tab), a Sequence Enrichment CSV, or any compatible CSV with `Sequence` and `Enrichment` columns.
2. If the file contains a cluster column, a **cluster selector** appears — choose which cluster to analyze. If no cluster column is detected, all sequences in the file are analyzed together.
3. Choose **Type of sequences**: Nucleotide or AminoAcid.
4. Click **Start**. Results appear in the download-ready output.
5. Click **Download** to save the positional enrichment CSV.

#### Plotting

**Position Enrichment Bar Plot**

Bar chart with alignment positions on the x-axis and mean enrichment at each position on the y-axis. Click **Position Enrichment** after a successful run.

Set **Adjust default bar plot?** to **Yes** for customization.

| Customization Option | Description |
|---|---|
| X-axis label | Label for the alignment position axis |
| Y-axis label | Label for the mean enrichment axis |
| Plot title | Title above the chart |
| Bar outline color | Color of bar borders |
| Bar fill color | Color of bar interiors |

**Position Enrichment Heatmap**

Heatmap with alignment positions on the x-axis, residue identities on the y-axis, and color intensity encoding mean enrichment at each position-residue combination. Click **Heatmap** after a successful run.

Set **Adjust default heat plot?** to **Yes** for customization.

| Customization Option | Description |
|---|---|
| X-axis label | Label for alignment positions |
| Y-axis label | Label for residue identities |
| Legend title | Label for the enrichment color scale |
| Plot title | Title above the heatmap |
| Fill palette | Color palette for the heatmap |

---

## 3.3 FASTAptameR-Motif

The Motif page contains four tabs: **Motif Discovery**, **Motif Search**, **Motif Omit**, and **Motif Tracker**.

---

### 3.3.1 Motif Discovery

#### Description

Scans a counted FASTA file to find over-represented sub-sequences using Z-score and p-value statistics. Enumerates all possible sub-sequences within a user-defined length range and tests each for enrichment relative to a background model. Outputs a ranked motif table and a significance scatter plot.

#### Usage

1. Upload a counted FASTA or previously generated motif discovery CSV.
2. To load an existing result directly, click **Visualize Results**.
3. Set sequence alphabet, minimum reads, and motif length range.
4. Click **Start**. Results appear in the table.
5. Click **Download** to save the output.

#### Parameters

| Parameter | Default | Range | Description |
|---|---|---|---|
| Sequence alphabet | — | Nucleotide / AminoAcid | Selects the character set for sub-sequence enumeration. |
| Min. number of reads to consider | — | 0–1,000 (step 5) | Sequences with fewer reads are excluded from the analysis. |
| Motif length range (Min / Max) | — | 3–15 (step 1) | Range of sub-sequence lengths to enumerate and test. |

#### Plotting

**Motif Discovery Plot**

Scatter plot with motif rank on the x-axis and −log₁₀(p-value) on the y-axis, colored by motif length. Click **Motif Discovery Plot** (requires results in the table).

Set **Adjust default plots?** to **Yes** for customization.

| Customization Option | Description |
|---|---|
| X-axis label | Label for motif rank axis |
| Y-axis label | Label for −log₁₀(p-value) axis |
| Legend title | Label for motif-length color legend |
| Plot title | Title above the plot |
| Color palette | Qualitative palette: Accent, Dark2, Paired, Pastel1, Pastel2, Set1, Set2 |

---

### 3.3.2 Motif Search

#### Description

Searches a counted FASTA file for sequences containing one or more user-specified motifs (nucleotide, amino acid, or arbitrary string), with AND/OR logic for multiple patterns. Outputs a filtered FASTA or CSV containing only the matching sequences.

#### Usage

1. Upload a counted FASTA file.
2. Enter one or more comma-separated patterns in the **Comma-separated patterns** field (e.g., `AAA,GTG,UUCG`).
3. Set pattern type, partial match behavior, and optional parenthesization.
4. Choose download format (FASTA or CSV).
5. Click **Start**. Matching sequences appear in the table.
6. Click **Download** to save the output.

#### Parameters

| Parameter | Default | Description |
|---|---|---|
| Comma-separated patterns | — | One or more motif strings to search for, separated by commas. Supports nucleotide (IUPAC), amino acid, or plain string patterns. |
| Place patterns in parentheses in output? | — | If Yes, matched motif occurrences in the output sequence string are wrapped in parentheses for visual identification. |
| If multiple patterns, return partial matches? | — | If Yes (OR logic): return sequences matching at least one pattern. If No (AND logic): return only sequences matching all patterns. |
| Type of pattern? | — | Nucleotide / AminoAcid / String — selects the character matching rules applied. |
| FASTA or CSV download? | FASTA | Output format. |

---

### 3.3.3 Motif Omit

#### Description

Removes sequences containing specified motifs from a counted FASTA file. The inverse of Motif Search — sequences that match any of the provided patterns are excluded from the output. Accepts comma-separated motif patterns.

#### Usage

1. Upload a counted FASTA file.
2. Enter comma-separated patterns to exclude.
3. Set partial match behavior and pattern type.
4. Choose download format.
5. Click **Start**. Filtered (non-matching) sequences appear in the table.
6. Click **Download** to save the output.

#### Parameters

| Parameter | Default | Description |
|---|---|---|
| Comma-separated patterns | — | Motif strings to use as exclusion filters, separated by commas. |
| If multiple patterns, return partial matches? | — | If Yes (OR logic): remove sequences matching at least one pattern. If No (AND logic): remove sequences matching all patterns. |
| Type of pattern? | — | Nucleotide / AminoAcid / String. |
| FASTA or CSV download? | FASTA | Output format. |

---

### 3.3.4 Motif Tracker

#### Description

Monitors how motifs or exact sequences rise or fall in frequency across multiple counted FASTA populations (representing selection rounds or time points). Upload files in order, specify your motifs or sequences, and get a line plot of frequency per population plus a summary table with enrichment values.

#### Usage

1. Click **Browse...** and select multiple counted FASTA files.
2. After upload, reorder files using drag-and-drop or arrow buttons. Use checkboxes to exclude specific files.
3. Enter motifs or sequences (one per line) in the **Motif or Sequence List** text area.
4. Optionally provide display aliases (one per line) in the **Alias List**.
5. Choose **Motif** or **Sequence** search mode.
6. If Motif mode, choose pattern type (Nucleotide / AminoAcid / String).
7. Click **Start Tracking**. Results appear in both the Tracker Results and Enrichment Results tables.
8. Click **Tracker Plot** to visualize frequency trajectories.
9. Click **Download Summary** or **Download Enrichment** to save each table.

#### Parameters

| Parameter | Description |
|---|---|
| File order | Drag and reorder files to define the population sequence (e.g., round 1 → round N). Unchecked files are excluded. |
| Motif or Sequence List | One query per line. Motifs are partial-match sub-sequences; sequences are full exact matches. |
| Alias List | Optional friendly names for each query, corresponding to the same line order. Used as labels in the plot and table. |
| Search for motifs or whole sequences? | Motif: reports frequency of each sub-sequence. Sequence: reports total RPU for each exact match. |
| Type of pattern? (Motif mode only) | Nucleotide / AminoAcid / String. |

#### Plotting

**Tracker Plot**

Line chart with populations on the x-axis and frequency (% for motifs, total RPU for sequences) on the y-axis. One line per query.

- Click **Tracker Plot** (requires results in the table).
- Set **Adjust default tracker plot?** to **Yes** for customization.

| Customization Option | Description |
|---|---|
| X-axis label | Label for the population (round) axis |
| Y-axis label | Label for frequency or RPU axis |
| Plot title | Title above the chart |
| Color palette | Qualitative palette: Accent, Dark2, Paired, Pastel1, Pastel2, Set1, Set2 |

---

## 3.4 FASTAptameR-Translate

### Description

Translates DNA sequences from a counted FASTA file into amino acid sequences using a selected standard genetic code or a customized codon table. Optionally merges identical translation products and re-ranks by total read count. Outputs a ranked table of translated amino acid sequences.

### Usage

1. Upload a counted FASTA file. To load an existing result, click **Visualize Results**.
2. Select open reading frame, merge option, and genetic code.
3. If using a non-standard codon table, enable customization and enter modified codons.
4. Choose download format.
5. Click **Start**. Translated sequences appear in the table.
6. Click **Download** to save the output.

### Parameters

| Parameter | Default | Description |
|---|---|---|
| Open reading frame | 1 | Reading frame offset: 1, 2, or 3. Determines which base to start translation from. |
| Should non-unique sequences be merged? | — | If Yes, identical translation products from different DNA sequences are combined; their read counts are summed and the sequence is re-ranked. |
| Genetic code selection | Standard | Dropdown of standard NCBI genetic codes (e.g., Standard, Vertebrate Mitochondrial, Bacterial). |
| Customize translation table? | No | If Yes, enter comma-separated codons and corresponding amino acid outputs to override the selected code. |
| Codons to modify | — | Comma-separated codon strings (e.g., `AAA,TTT`). Active only when customization is enabled. |
| Custom translations | — | Comma-separated single-letter amino acid outputs, in the same order as the codons above (e.g., `Z,X`). |
| FASTA or CSV download? | FASTA | Output format. |

### Plotting

Two plots are available after a successful run. Both share two pre-filter sliders:

| Slider | Range | Description |
|---|---|---|
| Min. number of reads to plot | 0–1,000 (step 10) | Exclude low-abundance sequences from plot |
| Max. rank to plot | 10–1,000 (step 10) | Limit plot to top N ranked sequences |

**Reads per Rank Plot**

Line chart of read abundance vs. rank for translated sequences. Click **Reads per rank**. Customization: X-axis label, Y-axis label, plot title, line color.

**Sequence-Length Histogram**

Histogram of unique sequence count (top) and total read count (bottom) by translated sequence length. Click **Sequence-length histogram**. Customization: axis labels, plot title, bar outline color, bar fill colors (unique and total).

---

## 3.5 FASTAptameR-Sequence Enrichment

### Description

Compares RPU values for sequences shared between two counted FASTA populations and computes log₂ enrichment (fold change in normalized abundance from population 1 to population 2). Provides a quantitative measure of which sequences are selected for or against. Outputs an enrichment table with log₂ enrichment scores.

### Usage

1. Upload **Input data 1** (earlier population, e.g., pre-selection or earlier round).
2. Upload **Input data 2** (later population, e.g., post-selection or later round).
3. Choose whether to keep sequences not found in both files.
4. Click **Start**. Results appear in the table.
5. Click **Download** to save the output CSV.

### Parameters

| Parameter | Default | Description |
|---|---|---|
| Keep missing sequences? | — | If Yes (union): include sequences found in only one population, with NA for missing values. If No (intersection): include only sequences present in both files. |

### Output Table Columns

| Column | Description |
|---|---|
| Sequence | Amino acid or nucleotide sequence |
| Reads 1 / RPU 1 | Read count and RPU from population 1 |
| Reads 2 / RPU 2 | Read count and RPU from population 2 |
| log2 Enrichment | log₂(RPU₂ / RPU₁) — positive values indicate selection, negative indicate depletion |

### Plotting

Three plots are available after a successful run.

**log₂(Enrichment) Histogram**

Histogram showing the distribution of log₂ enrichment values across sequences. Click **log2(enrichment) histogram**.

Set **Adjust default enrichment histogram?** to **Yes** for customization.

| Customization Option | Description |
|---|---|
| X-axis label | Label for log₂ enrichment bins |
| Y-axis label | Label for unique sequence count |
| Plot title | Title above the histogram |
| Bar outline color | Color of bar borders |
| Bar fill color | Color of bar interiors |

**RPU Scatter Plot**

Scatter plot with RPU from population 1 on the x-axis and RPU from population 2 on the y-axis. Each point is a sequence. Click **RPU scatter plot**.

| Customization Option | Description |
|---|---|
| X-axis label | Label for RPU (population 1) axis |
| Y-axis label | Label for RPU (population 2) axis |
| Plot title | Title above the scatter plot |
| Point color | Color of scatter points |

**RA Plot (MA Plot)**

Scatter plot with mean log₂(RPU) on the x-axis (A) and log₂ fold change on the y-axis (R). A classical representation for identifying differentially enriched sequences. Click **RA plot**.

| Customization Option | Description |
|---|---|
| X-axis label | Label for mean log₂(RPU) axis |
| Y-axis label | Label for fold change axis |
| Plot title | Title above the plot |
| Point color | Color of scatter points |

---

## 3.6 FASTAptameR-Data Merge

### Description

Combines 2 to 26 counted FASTA files into a single merged table, tracking sequence presence and abundance across all populations. Three join modes control which sequences are included. Outputs a merged CSV with one row per sequence and one column per population. Useful as a prerequisite for cross-round analysis and enrichment calculations.

### Usage

1. Click **Browse...** and select 2–26 counted FASTA files.
2. After upload, reorder files by dragging cards or using arrow buttons. Uncheck files to exclude them.
3. Select the **join type**: Union, Intersection, or Left.
4. Click **Start**. The merged table appears on the right.
5. Click **Download** to save the merged CSV.

### Parameters

| Parameter | Description |
|---|---|
| File order | Order of files determines column order in the output. Drag or use arrows to reorder; uncheck to exclude from the merge. |
| How should the data be joined? | **Union**: include all sequences from all files (missing values filled with NA). **Intersection**: include only sequences present in every file. **Left**: include only sequences from the first file (left join). |

### Plotting

**Persistence Plot**

Bar chart showing the number of unique sequences found in exactly 1, 2, 3, … N populations (where N = number of files). Shows sequence persistence across rounds. Click **Persistence plot** after running.

Set **Adjust default persistence plot?** to **Yes** for customization.

| Customization Option | Description |
|---|---|
| X-axis label | Label for "number of populations" axis |
| Y-axis label | Label for unique sequence count axis |
| Plot title | Title above the chart |
| Bar outline color | Color of bar borders |
| Bar fill color | Color of bar interiors |

**UpSet Plot**

Visualizes set intersections between all uploaded populations in an UpSet diagram format. Shows the number of unique sequences shared by specific combinations of populations. Click **UpSet plot** after running.

Set **Adjust default UpSet plot?** to **Yes** for customization.

| Customization Option | Description |
|---|---|
| X-axis label | Label for unique sequence count axis |
| Y-axis label | Label for the population sets axis |
| Bar fill color | Color of intersection bars |

---

## 3.7 FASTAptameR-Differential Analysis

### Description

Performs differential expression analysis between two experimental conditions, each represented by 2 or more counted FASTA files (replicates). Uses the edgeR statistical framework to compute log₂ fold change (logFC), log₂ counts per million (logCPM), and a Benjamini–Hochberg corrected p-value for each sequence. Identifies sequences that are significantly more or less abundant in one condition vs. another.

### Usage

1. Upload 2 or more FASTA files for **Condition 1** (e.g., treatment replicates).
2. Upload 2 or more FASTA files for **Condition 2** (e.g., control replicates).
3. Set the p-value cutoff threshold.
4. Click **Start**. Results appear in the table.
5. Click **Download** to save the output CSV.

> The Start button is disabled until at least 2 files per condition are uploaded.

### Parameters

| Parameter | Default | Range | Description |
|---|---|---|---|
| Cutoff for p-value | — | 0.01–0.50 (step 0.01) | Maximum BH-corrected p-value to flag as significant in the edgeR plot. |

### Output Table Columns

| Column | Description |
|---|---|
| Sequence | Nucleotide or amino acid sequence |
| logFC | Log₂ fold change between condition 2 and condition 1 |
| logCPM | Log₂ counts per million — overall abundance measure |
| p-value | BH-corrected p-value for differential abundance |

### Plotting

**edgeR Plot (MA-style)**

Scatter plot with logCPM on the x-axis and logFC on the y-axis. Significant sequences (below p-value cutoff) are colored differently from non-significant ones. Click **edgeR plot** (requires results in table).

Set **Adjust default edgeR plot?** to **Yes** for customization.

| Customization Option | Description |
|---|---|
| X-axis label | Label for log(counts per million) axis |
| Y-axis label | Label for log(fold change) axis |
| Plot title | Title above the scatter plot |
| Color for significant points | Color for sequences meeting the p-value cutoff |
| Color for insignificant points | Color for sequences above the p-value cutoff |

---

## 3.8 FASTAptameR-Distance

### Description

Computes the Levenshtein Edit Distance between a user-supplied query sequence and every sequence in a counted FASTA file. Outputs a complete table of all sequences ranked by their distance to the query. This module is useful for identifying sequences similar to a known candidate and for guiding LED threshold selection in FASTAptameR-Cluster.

### Usage

1. Upload a counted FASTA or previously generated distance CSV.
2. To load an existing result directly, click **Visualize Results**.
3. Type the **query sequence** into the text field.
4. Choose download format (FASTA or CSV).
5. Click **Start**. Results appear ranked by distance in the table.
6. Click **Download** to save the output.

> FASTA output is required for subsequent modules; CSV retains all columns from the data table.

### Parameters

| Parameter | Description |
|---|---|
| Query sequence | The reference sequence against which all population sequences are compared. Enter as a plain nucleotide or amino acid string. |
| Download format | FASTA or CSV. |

### Output Table Columns

| Column | Description |
|---|---|
| ID | Sequence identifier |
| LED | Levenshtein edit distance from the query sequence |
| Reads | Raw read count |
| Rank | Global rank in the original counted file |
| RPU | Reads per unit |
| Sequence | Nucleotide or amino acid sequence |

### Plotting

**Distance Histogram**

Dual bar chart showing unique sequence count (top) and total sequence count (bottom) binned by LED from the query. Click **Distance histogram** (requires results in table).

Set **Adjust default distance histograms?** to **Yes** for customization.

| Customization Option | Description |
|---|---|
| X-axis label | Label for LED from query axis |
| Y-axis 1 label | Label for unique sequence count |
| Y-axis 2 label | Label for total sequence count |
| Plot title | Title above the histogram |
| Bar outline color | Color of bar borders |
| Bar fill color | Color of bar interiors |

---

## 3.9 FASTAptameR-Mutation Network

### Description

Finds the shortest mutational path between a start sequence and an end sequence through the sequences present in a FASTA file. At each step of the path, the algorithm selects a sequence in the file that is within a maximum allowed edit distance from the current node and reduces the distance to the end sequence. Outputs the ordered sequence path as a table.

> This analysis may take a few moments. Remain on the page while it is running.

### Usage

1. Upload a counted FASTA file.
2. Enter the **Start sequence** (the sequence to begin the path from).
3. Enter the **End sequence** (the target sequence to reach).
4. Set the **Maximum allowable distance** per step using the slider.
5. Click **Start**. The path is shown in the table when complete.
6. Click **Download** to save the output.

### Parameters

| Parameter | Default | Range | Description |
|---|---|---|---|
| Start sequence | — | — | The source sequence for the mutational path. Must be a valid nucleotide or amino acid sequence string. |
| End sequence | — | — | The destination sequence for the mutational path. |
| Maximum allowable distance | — | 1–5 (step 1) | Maximum Levenshtein edit distance allowed per step along the path. Higher values allow larger jumps between steps. |

### Output Table Columns

| Column | Description |
|---|---|
| Step | Step number along the mutational path (1 = start) |
| Sequence | Sequence at this step |
| LED to Next | Edit distance to the next sequence in the path |
| Reads | Read count of the sequence in the population |

---

## 3.10 FASTAptameR-Pipeline

### Description

Chains multiple FASTAptameR operations — including Preprocess, Count, Translate, Motif Search, Motif Omit, Motif Discovery, Mutation Network, Sequence Distance, Cluster, Cluster MSA, Cluster PHMM, and Cluster Diversity — into a single automated workflow. Upload an initial file, define and configure steps, then execute the entire sequence with one click. Each step's output is automatically passed as input to the next. Step status, input/output filenames, and error messages are shown in real time in the right panel.

> Stay on this page while the pipeline is running. All temporary files are deleted when the session is closed or refreshed.

### Usage

1. Click **Browse...** to upload an initial FASTA or FASTQ file.
2. Click **Add Step** to add operations to the pipeline. Each step appears as a card.
3. For each step, select an operation from the dropdown and configure parameters.
4. Reorder or remove steps as needed.
5. Click **Run Pipeline**. Steps execute sequentially; the right panel shows real-time status.
6. Click **Download** on any completed step to save its output.
7. Click **Cancel** to stop the pipeline mid-run.

### Available Operations and Parameters

| Operation | Key Parameters | Valid Inputs |
|---|---|---|
| Preprocess | 5′ region, 3′ region, min/max length, max error | FASTA / FASTQ |
| Count | Reverse complement, scaling factor | Preprocessed or raw FASTA/FASTQ |
| Translate | Open reading frame (1/2/3), merge identical (yes/no), genetic code | Counted FASTA |
| Motif Search | Motifs (comma-separated), highlight, partial match, motif type | Counted FASTA |
| Motif Omit | Motifs (comma-separated), partial match, motif type | Counted FASTA |
| Motif Discovery | Min reads, min/max motif length, alphabet (nucleotide/protein) | Counted FASTA |
| Mutation Network | Start sequence, end sequence, max edit distance | Counted FASTA |
| Sequence Distance | Query sequence | Counted FASTA |
| Cluster | Min reads, max LED, max clusters, keep NC (yes/no) | Counted FASTA |
| Cluster MSA | Cluster ID, sequence type (nucleotide/protein) | Clustered FASTA |
| Cluster PHMM | Number of sequences to simulate, target sequence length | Cluster MSA FASTA |
| Cluster Diversity | (no parameters) | Clustered FASTA |

### Step Ordering Rules

The pipeline enforces valid step sequences. Some operations are only valid following specific prior steps:

- **Count** must follow **Preprocess** (or be the first step if starting from a counted file is not intended).
- **Cluster** accepts output from **Count**, **Translate**, **Motif Search**, or **Motif Omit**.
- **Cluster MSA** and **Cluster Diversity** must follow **Cluster**.
- **Cluster PHMM** must follow **Cluster MSA**.

Invalid step combinations are flagged with a warning icon and an error message on the step card.

### Pipeline Parameters (by operation)

**Preprocess**

| Parameter | Default | Description |
|---|---|---|
| Constant 5′ Region | (blank) | Sequence to trim from 5′ end |
| Constant 3′ Region | (blank) | Sequence to trim from 3′ end |
| Min Length | 10 | Minimum retained sequence length after trimming |
| Max Length | 100 | Maximum retained sequence length after trimming |
| Max Allowed Error | 0.005 | Maximum mean error probability per sequence (0–1) |

**Count**

| Parameter | Default | Options | Description |
|---|---|---|---|
| Reverse Complement | false | true / false | Report reverse complement of sequences |
| Scaling Factor | 1,000,000 | 1, 10, 100, 1,000, 10,000, 100,000, 1,000,000 | RPU normalization denominator |

**Translate**

| Parameter | Default | Options | Description |
|---|---|---|---|
| Open Reading Frame | 1 | 1 / 2 / 3 | Translation start offset |
| Merge Identical Translations | true | true / false | Combine identical protein products |
| Genetic Code | Standard | NCBI genetic code list | Codon table for translation |

**Cluster**

| Parameter | Default | Description |
|---|---|---|
| Min Reads | 10 | Minimum reads for a sequence to enter clustering |
| Max LED | 7 | Maximum edit distance for cluster membership |
| Max Clusters to Generate | 20 | Maximum number of clusters to produce |
| Keep Non-Clustered | false | Whether to retain NC-labeled sequences |

**Cluster MSA**

| Parameter | Default | Options | Description |
|---|---|---|---|
| Cluster ID | 1 | — | The cluster number to align |
| Sequence Type | dna | dna / protein | Alignment substitution model |

**Cluster PHMM**

| Parameter | Default | Description |
|---|---|---|
| Number of Sequences to Simulate | 100 | How many sequences to generate from the PHMM |
| Target Sequence Length | 50 | Target length including gaps |
