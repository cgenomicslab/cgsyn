# CGSyn
***A Comprehensive, Robust Tool for a Step-By-Step Comparative Macrosynteny Analysis***

Comparative macrosynteny analysis is the study of conserved gene order and chromosomal organization across species. While individual genes mutate and diverge rapidly, the large-scale arrangement of genes along chromosomes can remain stable over hundreds of millions of years of evolution. By comparing which genes co-localize on the same chromosomes across distantly related species, we can reconstruct the chromosomal architecture of common ancestors and identify large-scale genomic rearrangements.

The approach has gained significant traction in phylogenomics following two landmark studies. [Simakov et al. (2022)](https://doi.org/10.1126/sciadv.abi5884) identified Ancestral Linkage Groups (ALGs) conserved across bilaterians, cnidarians and sponges, demonstrating that chromosomal organization can be traced across over 600 million years of animal evolution. [Schultz et al. (2023)](https://doi.org/10.1038/s41586-023-05936-6) extended this framework to early-branching animal lineages, using macrosynteny patterns to address the contested phylogenetic position of ctenophores - providing genomic evidence that sequence-based phylogenetics alone had struggled to resolve.

These studies established macrosynteny as a genuinely independent line of phylogenomic evidence. CGSyn was built to make this type of analysis accessible to any research group working with chromosome-level genome assemblies.

## Getting Started
CGSyn currently only works on Unix-based systems (Linux/macOS).

#### Install via command line
```bash
git clone https://github.com/cgenomicslab/CGSyn.git
cd CGSyn
```

#### Install via browser
1. Click on the blue "<> Code" button at the top of the repository page
2. Select "Download ZIP"
3. Extract the ZIP file to your desired location


### Set Up Environment
1. Open a terminal and navigate to the project directory:
```bash
cd CGSyn
```
2. Create a conda environement with all the required software for the tool and activate it
```bash
conda env create -f workflow/envs/synteny.yaml -n <name>
conda activate <name>
```
3. To see available options and workflows, run
```bash
./synteny.sh --help
```

## Repo Tree Organization

```md
.
├── config
│   └── config.yaml
├── intermediates
│   ├── filtered_proteomes
│   └── tsv
├── logs
├── resources
│   ├── gff
│   └── proteomes
├── results
│   ├── alg_orthofinder
│   ├── alg_rbh
│   ├── dotplots_orthofinder
│   ├── dotplots_rbh
│   ├── gene_analysis
│   ├── orthofinder
│   ├── rbh
│   ├── ribbons_orthofinder
│   ├── ribbons_rbh
│   ├── ribbons_multi_orthofinder
│   └── ribbons_multi_rbh
├── synteny.sh
└── workflow
    ├── envs
    │   └── synteny.yaml
    ├── scripts
    │   ├── ncbi_download.py
    │   ├── run_alg_discovery_orthofinder.py
    │   ├── run_alg_discovery_rbh.py
    │   ├── run_dotplot_orthofinder.py
    │   ├── run_dotplot_rbh.py
    │   ├── run_gene_analysis.py
    │   ├── run_parsing.py
    │   ├── run_ribbon_orthofinder.py
    │   ├── run_ribbon_rbh.py
    │   ├── run_ribbons_multi_orthofinder.py
    │   ├── run_ribbons_multi_rbh.py
    │   └── synteny.py
    └── Snakefile
```

## Snakemake Workflow

<img width="572" height="900" alt="SyntenyAnalysisWorkflow" src="https://github.com/user-attachments/assets/f12e27c3-d2f2-4fe5-835c-12aaaf91d1a4" />

## All Usage Functions

### Step 1: Downloading Resources

This tool uses proteome (.fasta) and functional annotation (.gff) files as primary resources for all downstream analysis.

You can either download these files manually from [NCBI Genome](https://www.ncbi.nlm.nih.gov/datasets/genome/) or let the tool do it for you.

- Example:
```bash
./synteny.sh --download --species-queries "9606,Mus musculus" --species-labels "Hsap,Mmus"
```
The "species-queries" flag can take either the species' Tax ID or its scientific name. The "species-labels" flag renames the proteome and annotation 
files with your preferred labels (e.g. Hsap.faa.gz, Hsap.gff.gz). While optional, it is highly recommended you utilize this flag to set easily distinguishable, 
as well as publishable species names, since those labels will be used in the tables and figures created by all downstream analyses.

Note❗: The tool will NOT redownload already existing files, unless you download the same assembly with a different species label.

The tool follows a hierarchical search order (see below), prioritizing Complete and Chromosome-Level Assemblies in RefSeq (preferred) or Genbank, WITH existing annotation files.
If none exist, it will move to Scaffold- and Contig-Level ones, but will ask for user permission after informing of the amount of scaffolds/contigs in the assembly.

WARNING ⚠️: You CANNOT perform a successful comparative synteny analysis without an annotated genome with well-assembled scaffolds. If the scaffolds intuitively
seem too many, it is very likely that the assembly has not successfully recreated the organism's chromosomes, creating a false image of the syntenic conservation.

**Hierarchical Search Order:**
1. Complete genome in RefSeq
2. Chromosome-level in RefSeq
3. Complete genome in GenBank (with user confirmation)
4. Chromosome-level in GenBank (with user confirmation)
5. Scaffold-level in RefSeq (with user confirmation + warning)
6. Scaffold-level in GenBank (with user confirmation + warning)
7. Contig-level in RefSeq (with user confirmation + warning)
8. Contig-level in GenBank (with user confirmation + warning)

All downloaded files can be found in ```./resources/proteomes``` and ```./resources/gff```. If you download your own assemblies manually, make sure to move them in the
correct folders and rename them with a readable, publication-ready species label.

### Step 2: File Parsing and Data Cleaning

By running:
```bash
./synteny.sh --species <sp1_label,sp2_label,sp3_label,...> --parse
```

CGSyn will parse the gff and proteome files and
1. extract the genome coordinates (chromosome/scaffold, start coord, end coord, strand) of every gene annotated in the assembly
2. filter through the isoforms produced by it, in order to keep the ProteinID of its longest isoform
3. clean the proteome files to remove all isoforms apart from the longest for every gene

The new filtered proteomes will be put in ```./intermediates/filtered_proteomes```, while .tsv files containing the following columns: ```GeneID, ProteinID, chr, start, end, strand```
will be created and put in ```./intermediates/tsv```. The former will be used in Orthology Inference, while the latter will be instrumental in the post-inference analyses.

### Step 3: Orthology Inference

There are 2 alternatives for inferring gene orthology:
1. Using the [Orthofinder](https://github.com/OrthoFinder/OrthoFinder) software by David Emms
```bash
./synteny.sh --species <sp1,sp2,sp3,...> --orthofinder --aligner <diamond or blastp, default: diamond>
```
2. Using a Reciprocal Best Hits (RBH) algorithm
```bash
./synteny.sh --species <sp1,sp2,sp3,...> --rbh --aligner <diamond or blastp, default: diamond>
```

**Orthofinder** runs an all-vs-all similarity search between the proteins of the query organisms, produces a sequence similarity graph and runs Markov Clustering to 
infer Orthogroups (Gene Families across all species which seem to be descended from the same common ancestor), as well as orthologous pairs of genes between pairs of Organisms.

**RBH** runs a similarity search between all proteins of species A against all proteins of species B and vice versa and keeps the reciprocal best hits as orthologous pairs.
It then does this for all possible pairs of species (Multiple Reciprocal Best Hits - MBH) to create a cluster of proteins that supposedly belong to the same gene family.

RBH is much faster, even with many species, but Orthofinder finds more orthologous pairs and produces a lot more outputs for further use, including gene trees and a species tree.
Therefore, Orthofinder is suggested.

Similarly, using Diamond instead of BLASTP as an aligner for both inference methods is suggested due to its significantly quicker running time.

Orthology Inference results will be saved in ```./results/orthofinder``` and ```./results/rbh``` respectively.

### Step 4: Dot Plots and Ribbon Diagrams Creation

**Oxford Dot Plots** are a standard tool in comparative genomics for visualizing synteny between two species. Each dot represents a pair of orthologous genes, placed according to their genomic position in each species - species 1 on the X axis and species 2 on the Y axis. Chromosomes are displayed sequentially along each axis, separated by gridlines.

Only genes from chromosome pairs that pass a Fisher's exact test for significant ortholog enrichment (Bonferroni-corrected, α = 0.01 by default - can be changed in config file) are colored; 
all other dots are shown in gray. We assume these chromosome pairs have a conserved synteny, since they share more orthologous pairs than they would by random chance.
The colors themselves represent the chromosomes of species 1 and have no further significance. 

In Dot Plots, chromosomes are ordered nominally (e.g. 1-10, I-VII, A-F etc).

To create dot plots, you can run either of these, depending on if you want to use the inference results from Orthofinder or RBH.

```bash
./synteny.sh --species <sp1,sp2,sp3,...> --dotplots-orthofinder
```
or
```bash
./synteny.sh --species <sp1,sp2,sp3,...> --dotplots-rbh
``` 
The plots are saved in .png format in the ```./results/dotplots_orthofinder``` and ```./results/dotplots_rbh``` directories.

**Synteny Ribbon Diagrams** are visual representations used in comparative genomics to illustrate the conservation of gene order and large-scale evolutionary relationships across multiple genomes. They highlight structural rearrangements—such as inversions and translocations—by connecting homologous chromosomal regions with colored, curved "ribbons"

Similarly to Oxford Dot Plots, only ribbons connecting genes from chromosome pairs that pass a Fisher's exact test for significant ortholog enrichment
(Bonferroni-corrected, α = 0.01 by default - can be changed in config file) are colored, with each color simply representing one of the chromosomes of species 1.

In Ribbon Plots, chromosomes of species 1 are ordered nominally, while chromosomes of species 2 are ordered in a way that will create the least amount of curved (Bézier) ribbons.

To create ribbon diagrams, you can run either of these, depending on if you want to use the inference results from Orthofinder or RBH:

```bash
./synteny.sh --species <sp1,sp2,sp3,...> --ribbons-orthofinder
```
or
```bash
./synteny.sh --species <sp1,sp2,sp3,...> --ribbons-rbh
```
The plots are saved in .png format in the ```./results/ribbons_orthofinder``` and ```./results/ribbons_rbh``` directories.

### Step 5: Ancestral Linkage Group Discovery and Multi-Species Ribbon Diagram

**Ancestral Linkage Groups (ALGs)** are sets of genes that were physically linked on the same chromosome in the last common ancestor of the species being compared, 
and have remained co-localized across evolutionary time. Identifying ALGs allows us to reconstruct the ancestral chromosomal architecture of a lineage 
and understand how chromosomes have been broken, fused or rearranged since that ancestor.

CGSyn's ALG discovery algorithm takes a multi-species approach to identifying these conserved chromosomal units. It works as follows:

1. Multi-species Fisher's test: Fisher's exact test is run for every possible pair of species simultaneously, identifying which chromosome-to-chromosome relationships share significantly more orthologs than expected by chance. This produces a filtered ortholog map tracking which species pairs each gene is significant in.
2. Synteny similarity matrix: For each pair of species, the fraction of their shared orthologs that fall in statistically significant chromosome pairs is computed, producing an N×N similarity matrix.
3. Species clustering: Species are grouped into synteny clusters via hierarchical clustering on the similarity matrix. Species with strong conserved synteny (e.g. two species from the same phylum) will cluster together, while distantly related species will form separate clusters (default similarity threshold = 0.3, can be changed in config).
4. Chain building: For each cluster, all possible chromosome chains are constructed by following significant adjacent pairwise matches across species (e.g. Hsap1 → Mmus1 → Ggal5 → Bflor17).
5. Chain verification: Each candidate chain is verified by checking that every pairwise combination of species in the chain has a statistically significant chromosome match - not just adjacent ones. Chains that fail any pairwise check are rejected.

You can run the ALG discovery algorithm with:

```bash
./synteny.sh --species <sp1,sp2,sp3,...> --alg-discovery-orthofinder
```
or
```bash
./synteny.sh --species <sp1,sp2,sp3,...> --alg-discovery-rbh
```
The output is a set of ALG assignments per species per chromosome, prefixed by cluster (e.g. C2_ALG1), saved as both a machine-readable .pkl file and a human-readable summary .txt file, saved in the ```./results/alg_orthofinder``` and ```./results/alg_rbh``` directories.

**Multi-Species Ribbon Diagrams** extend the pairwise ribbon plot concept to N species simultaneously. Species are arranged vertically, with chromosomes displayed as horizontal lines for each species. Individual Bezier curve ribbons connect each orthologous gene across adjacent species pairs.
If ALG discovery has been run, ribbons are colored by ALG identity, making it visually straightforward to trace ancestral chromosomal units across all species in the analysis. Ribbons between species in different synteny clusters are suppressed, and a red dashed line marks cluster boundaries. Colored rectangles below the bottom species of each cluster serve as an ALG legend.
If ALG discovery has not been run (not suggested), ribbons are colored by the chromosomes of the first species in the --species list.

You can create multi-species ribbons diagrams with:

```bash
./synteny.sh --species <sp1,sp2,sp3,...> --ribbons-multi-orthofinder
```
or
```bash
./synteny.sh --species <sp1,sp2,sp3,...> --ribbons-multi-rbh
```
The plots are saved in .png format in the ```./results/ribbons_multi_orthofinder``` and ```./results/ribbons_multi_rbh``` directories.

### Extra: Anchor Gene Family Analysis

**Anchor Gene Families** are sets of genes that are present as single-copy orthologs across all species being compared and show conserved synteny. The concept was introduced by [Simakov et al. (2022)](https://doi.org/10.1126/sciadv.abi5884) as a metric to support the validity of the ALGs he discovered and anchor synteny comparisons across deeply diverged species - hence the name.

To qualify as an anchor gene family, a gene must satisfy two strict criteria:
* Single-copy: the gene must appear exactly once in each species - no duplications allowed
* Mutual best hit: each gene in the family must be the reciprocal best BLAST/DIAMOND hit of its ortholog in every other species.

These two filters together ensure that the genes are truly orthologous (not paralogs) and have not undergone lineage-specific expansions, making them reliable landmarks for synteny comparison.

CGSyn implements two modes:
* **Triangle mode (default)**: A set of 3 core species (the "triangle") is defined first. Only genes that form a 3-way mutual best hit across all three are kept as the core set. Additional species are then added one by one - a gene joins a family only if it is a mutual best hit with all 3 core species. This mode is more lenient and recommended when the species are distantly related.

```bash
./synteny.sh --triangle <sp1,sp2,sp3> --rest <sp4,sp5> --gene-analysis
```

* **Full MBH mode**: All species are treated equally from the start, requiring mutual best hits across every possible pair simultaneously. More stringent than triangle mode as the number of species increases.


```bash
./synteny.sh --rest <sp1,sp2,sp3,sp4,sp5> --mbh --gene-analysis
```

Results are saved as a .tsv table with one row per family and one column per species, and a human-readable .txt file listing the protein IDs and chromosomal locations of each family member in the ```results/gene_analysis``` directory.

## Usage Examples
A user can run multiple flags at the same time, as long as they all belong to the same Orthology Inference Pathway (Orthofinder vs RBH)

Orthofinder Full Analysis Example:

```bash
./synteny.sh --download --species-queries <"species1,species2,species3,..."> --species-labels <"sp1,sp2,sp3,..."> --species <sp1,sp2,sp3,...> --parse --orthofinder --dotplots-orthofinder --ribbons-orthofinder --ribbons-multi-orthofinder --alg-discovery-orthofinder --triangle <sp1,sp2,sp3> --rest <sp4,sp5> --gene-analysis --cores 32
```

RBH Full Analysis Example:

```bash
./synteny.sh --download --species-queries <"species1,species2,species3,..."> --species-labels <"sp1,sp2,sp3,..."> --species <sp1,sp2,sp3,...> --parse --rbh --dotplots-rbh --ribbons-rbh --ribbons-multi-rbh --alg-discovery-rbh --triangle <sp1,sp2,sp3> --rest <sp4,sp5> --gene-analysis --cores 32
```

## Citations:
1. Simakov, O. et al. Deeply conserved synteny and the evolution of metazoan chromosomes. Sci. Adv. 8, eabi5884 (2022). https://doi.org/10.1126/sciadv.abi5884
2. Schultz, D.T., Haddock, S.H.D., Bredeson, J.V. et al. Ancient gene linkages support ctenophores as sister to other animals. Nature 618, 110–117 (2023). https://doi.org/10.1038/s41586-023-05936-6
3. Emms, D.M. & Kelly, S. OrthoFinder: phylogenetic orthology inference for comparative genomics. Genome Biology 20, 238 (2019). https://doi.org/10.1186/s13059-019-1832-y 
4. Emms, D.M. & Kelly, S. OrthoFinder: solving fundamental biases in whole genome comparisons dramatically improves orthogroup inference accuracy. Genome Biology 16, 157 (2015). https://doi.org/10.1186/s13059-015-0721-2
