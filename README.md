# CGSyn
A Comprehensive, Robust Tool for a Step-By-Step Comparative Macrosynteny Analysis

[CGSyn was built with the purpose of giving all biologists the tools to run a comparative macrosynteny analysis with one line of code directly from their terminal]

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

<img width="1270" height="2000" alt="SyntenyAnalysisWorkflow" src="https://github.com/user-attachments/assets/f12e27c3-d2f2-4fe5-835c-12aaaf91d1a4" />

## All Usage Functions

### Step 1: Downloading Resources

This tool uses proteome (.fasta) and functional annotation (.gff) files as primary resources for all downstream analysis.

You can either download these files manually from [NCBI Genome](https://www.ncbi.nlm.nih.gov/datasets/genome/) or let the tool do it for you.

- Example:
```bash
./synteny.sh --download --species-queries "9606,Mus musculus” --species-labels "Hsap,Mmus"
```
The "species-queries" flag can take either the species' Tax ID or its scientific name. The "species-labels" flag renames the proteome and annotation 
files with your preferred labels (e.g. Hsap.faa.gz, Hsap.gff.gz). While optional, it is highly recommended you utilize this flag to set easily distinguishable, 
as well as publishable species names, since those labels will be used in the tables and figures created by all downstream analyses.

Note: The tool will NOT redownload already existing files, unless you download the same assembly with a different species label.

The tool follows a hierarchical search order (see below), prioritizing Complete and Chromosome-Level Assemblies in RefSeq (preferred) or Genbank, WITH existing annotation files.
If none exist, it will move to Scaffold- and Contig-Level ones, but will ask for user permission after informing of the amount of scaffolds/contigs in the assembly.

WARNING: You CANNOT perform a successful comparative synteny analysis without an annotated genome with well-assembled scaffolds. If the scaffolds intuitively
seem too many, it is very likely that the assembly has not successfully recreated the organism's chromosomes, creating a false image of in the syntenic conservation.

Hiearchical Search Order:
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

Orthofinder runs an all-vs-all similarity search between the proteins of the query organisms, produces a sequence similarity graph and runs Markov Clustering to 
infer Orthogroups (Gene Families across all species which seem to be descended from the same common ancestor), as well as orthologous pairs of genes between pairs of Organisms.

RBH runs a similarity search between all proteins of species A against all proteins of species B and vice versa and keeps the reciprocal best hits as orthologous pairs.
It then does this for all possible pairs of species (Multiple Reciprocal Best Hits - MBH) to create a cluster of proteins that supposedly belong to the same gene family.

RBH is much faster, even with many species, but Orthofinder finds more orthologous pairs and produces a lot more outputs for further use, including gene trees and a species tree.
Therefore, Orthofinder is suggested.

Similarly, using Diamond instead of BLASTP as an aligner for both inference methods is suggested due to its significantly quicker running time.

Orthology Inference results will be saved in ```./results/orthofinder``` and ```./results/rbh``` respectively.

### Step 4. 
