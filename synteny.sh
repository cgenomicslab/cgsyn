#!/usr/bin/env bash

# synteny - Command-line wrapper for Synteny Pipeline
# Usage: ./synteny [OPTIONS]

PYTHON=$(which python3)
set -e

# Default config file
CONFIG_FILE="config/config.yaml"
TEMP_CONFIG="/tmp/synteny_config_$$.yaml"

# Copy original config
cp "$CONFIG_FILE" "$TEMP_CONFIG"

# Parse command-line arguments
SNAKEMAKE_ARGS=()
CORES=1

show_help() {
    cat << EOF
Synteny Analysis Pipeline

Usage: ./synteny [OPTIONS]

Genome Download (NCBI):
  --download                      Enable download mode
  --species-queries QUERIES       Comma-separated species names or taxids
                                  Example: "9606,Mus musculus,Drosophila melanogaster"
  --species-labels LABELS         Comma-separated output labels (must match queries)
                                  Example: "Hsapiens,Mmusculus,Dmelanogaster"
  --non-interactive               Skip user prompts (auto-decisions)

Download Priority:
  1. Complete genome (RefSeq)
  2. Chromosome-level (RefSeq)
  3. Complete genome (GenBank) - requires confirmation
  4. Chromosome-level (GenBank) - requires confirmation
  5. Scaffold-level (RefSeq) - requires confirmation + warning
  6. Scaffold-level (GenBank) - requires confirmation + warning
  7. Contig-level (RefSeq) - requires confirmation + STRONG warning
  8. Contig-level (GenBank) - requires confirmation + STRONG warning
  
  Note: Contig-level assemblies are NOT recommended for synteny analysis
      due to extreme fragmentation. Use only as a last resort. Scaffold-level
      assemblies are not ideal either. Always check the number of scaffolds!

Examples:
  # Download single genome
  ./synteny --download --species-queries "Homo sapiens" --species-labels Hsapiens

  # Download multiple genomes (mix names and taxids)
  ./synteny --download \\
    --species-queries "9606,Mus musculus,7227" \\
    --species-labels human,mouse,fly

  # Non-interactive batch download
  ./synteny --download --non-interactive \\
    --species-queries "Homo sapiens,Mus musculus" \\
    --species-labels Hsapiens,Mmusculus

  # Then run pipeline
  ./synteny --species Hsapiens,Mmusculus,fly --orthofinder --ribbons-multi --cores 16

Species Selection:
  --species SP1,SP2,SP3   Comma-separated species list

Workflow Control:
  --parse                         Run GFF/proteome parsing
  --no-parse                      Skip parsing (use existing intermediates)
  --project NAME                  Use a separate results_NAME/ directory for
                                  this analysis, keeping it isolated from other
                                  projects. Resources and intermediates are
                                  shared across projects. If not set, defaults
                                  to "results/".
  --orthofinder                   Run OrthoFinder
  --rbh                           Run RBH/MBH orthology inference
  --compare-methods       	  Compare OrthoFinder vs RBH orthology inference results
  --dotplots-orthofinder          Generate Oxford dot plots from OrthoFinder
  --dotplots-rbh                  Generate Oxford dot plots from RBH
  --ribbons-orthofinder           Generate pairwise ribbon plots from OrthoFinder
  --ribbons-rbh                   Generate pairwise ribbon plots from RBH
  --ribbons-multi-orthofinder     Generate multi-species ribbon plot from OrthoFinder
                                  (ribbons colored by ALG if --alg-discovery-orthofinder is set,
                                  otherwise colored by first species chromosomes)
  --ribbons-multi-rbh             Generate multi-species ribbon plot from RBH
                                  (ribbons colored by ALG if --alg-discovery-rbh is set,
                                  otherwise colored by first species chromosomes)
  --alg-discovery-orthofinder     Run ALG discovery using OrthoFinder results
  --alg-discovery-rbh             Run ALG discovery using RBH results


Parameters:
  --threads N          Number of threads (default: from config)
  --aligner TOOL       diamond or blast (default: diamond)
  --inflation VALUE               MCL inflation parameter for OrthoFinder (default: 1.2)
                                  Higher = more, smaller orthogroups
                                  Lower  = fewer, larger orthogroups
  --tree-method METHOD            Gene tree inference method (default: msa)
                                  Options: msa, dendroblast
  --msa-program PROGRAM           MSA program, only used with --tree-method msa (default: famsa)
                                  Options: famsa, muscle, mafft
  --tree-inference METHOD         Tree inference method, only used with --tree-method msa (default: fasttree)
                                  Options: fasttree, fasttree_fastest, raxml, iqtree3
  --alpha VALUE        		  Fisher's test alpha (default: 0.01)
  --color-nonsignificant          Color non-significant dots by chromosome instead of grey in dot plots
  --no-cluster                    Treat all species as one group, skip synteny clustering
  --similarity-threshold VALUE    Similarity threshold for species clustering (default: 0.3)
                                  Species pairs below this threshold go into different clusters
  --shared-ogs                    Use shared orthogroup counts instead of strict 1-to-1
                                  ortholog pairs for Fisher's exact test. Increases
                                  sensitivity for distantly related species. OrthoFinder
                                  only - incompatible with --rbh.
  --cb-colors                     Use colorblind-safe color palette for all plots
                                  Based on Wong (2011) and Paul Tol's color schemes

Snakemake Options:
  --cores N            Number of cores for Snakemake (default: 1)
  --dry-run            Show what would be executed
  --forcerun           Force rerun all rules
  --unlock             Unlock working directory

Examples:
  # Run full pipeline with OrthoFinder
  ./synteny --orthofinder --dotplots --threads 16 --cores 16

  # Run RBH with Simakov analysis
  ./synteny --rbh --simakov --threads 12 --cores 12

  # Dry run to see what will execute
  ./synteny --orthofinder --dotplots --dry-run

  # Use specific species
  ./synteny --species Bfloridae,Pmaximus,Emuelleri --rbh --cores 8

EOF
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        --download)
            DOWNLOAD_MODE=true
            shift
            ;;
        --species-queries)
            SPECIES_QUERIES="$2"
            shift 2
            ;;
        --species-labels)
            SPECIES_LABELS="$2"
            shift 2
            ;;
        --non-interactive)
            NON_INTERACTIVE="--non-interactive"
            shift
            ;;
        --parse)
            sed -i 's/run_parsing: .*/run_parsing: true/' "$TEMP_CONFIG"
            shift
            ;;
        --no-parse)
            sed -i 's/run_parsing: .*/run_parsing: false/' "$TEMP_CONFIG"
            shift
            ;;
        --project)
            PROJECT="$2"
            sed -i "s|results_dir: .*|results_dir: \"results_${PROJECT}\"|" "$TEMP_CONFIG"
            shift 2
            ;;
        --orthofinder)
            sed -i 's/run_orthofinder: .*/run_orthofinder: true/' "$TEMP_CONFIG"
            sed -i 's/run_rbh: .*/run_rbh: false/' "$TEMP_CONFIG"
            ORTHOFINDER_REQUESTED=true
            shift
            ;;
        --inflation)
            INFLATION="$2"
            sed -i "s/inflation: .*/inflation: $INFLATION/" "$TEMP_CONFIG"
            shift 2
            ;;
        --tree-method)
            TREE_METHOD="$2"
            sed -i "s/tree_method: .*/tree_method: \"$TREE_METHOD\"/" "$TEMP_CONFIG"
            shift 2
            ;;
        --msa-program)
            MSA_PROGRAM="$2"
            sed -i "s/msa_program: .*/msa_program: \"$MSA_PROGRAM\"/" "$TEMP_CONFIG"
            shift 2
            ;;
        --tree-inference)
            TREE_INFERENCE="$2"
            sed -i "s/tree_inference: .*/tree_inference: \"$TREE_INFERENCE\"/" "$TEMP_CONFIG"
            shift 2
            ;;
        --rbh)
            sed -i 's/run_rbh: .*/run_rbh: true/' "$TEMP_CONFIG"
            sed -i 's/run_orthofinder: .*/run_orthofinder: false/' "$TEMP_CONFIG"
            RBH_REQUESTED=true
            shift
            ;;
        --compare-methods)
            sed -i 's/run_compare_methods: .*/run_compare_methods: true/' "$TEMP_CONFIG"
            shift
            ;;
        --dotplots-orthofinder)
            sed -i 's/run_dotplots_orthofinder: .*/run_dotplots_orthofinder: true/' "$TEMP_CONFIG"
            shift
            ;;
        --dotplots-rbh)
            sed -i 's/run_dotplots_rbh: .*/run_dotplots_rbh: true/' "$TEMP_CONFIG"
            shift
            ;;
        --color-nonsignificant)
            sed -i 's/color_nonsignificant: .*/color_nonsignificant: true/' "$TEMP_CONFIG"
            shift
            ;;
        --ribbons-orthofinder)
            sed -i 's/run_ribbons_orthofinder: .*/run_ribbons_orthofinder: true/' "$TEMP_CONFIG"
            shift
            ;;
        --ribbons-rbh)
            sed -i 's/run_ribbons_rbh: .*/run_ribbons_rbh: true/' "$TEMP_CONFIG"
            shift
            ;;
        --alg-discovery-orthofinder)
            sed -i 's/run_alg_discovery_orthofinder: .*/run_alg_discovery_orthofinder: true/' "$TEMP_CONFIG"
            shift
            ;;
        --ribbons-multi-orthofinder)
            sed -i 's/run_ribbons_multi_orthofinder: .*/run_ribbons_multi_orthofinder: true/' "$TEMP_CONFIG"
            shift
            ;;
        --alg-discovery-rbh)
            sed -i 's/run_alg_discovery_rbh: .*/run_alg_discovery_rbh: true/' "$TEMP_CONFIG"
            shift
            ;;
        --ribbons-multi-rbh)
            sed -i 's/run_ribbons_multi_rbh: .*/run_ribbons_multi_rbh: true/' "$TEMP_CONFIG"
            shift
            ;;
        --no-cluster)
            sed -i 's/cluster_species: .*/cluster_species: false/' "$TEMP_CONFIG"
            shift
            ;;
        --similarity-threshold)
            SIM_THRESHOLD="$2"
            sed -i "s/similarity_threshold: .*/similarity_threshold: $SIM_THRESHOLD/" "$TEMP_CONFIG"
            shift 2
            ;;
        --shared-ogs)
            sed -i 's/shared_ogs: .*/shared_ogs: true/' "$TEMP_CONFIG"
            shift
            ;;
        --cb-colors)
            sed -i 's/cb_colors: .*/cb_colors: true/' "$TEMP_CONFIG"
            shift
            ;;
        --gene-analysis)
            sed -i 's/run_gene_analysis: .*/run_gene_analysis: true/' "$TEMP_CONFIG"
            shift
            ;;
        --threads)
            THREADS="$2"
            sed -i "s/threads: .*/threads: $THREADS/" "$TEMP_CONFIG"
            shift 2
            ;;
        --alpha)
            ALPHA="$2"
            sed -i "s/alpha: .*/alpha: $ALPHA/" "$TEMP_CONFIG"
            shift 2
            ;;
        --aligner)
            ALIGNER="$2"
            sed -i "s/aligner: .*/aligner: \"$ALIGNER\"/" "$TEMP_CONFIG"
            shift 2
            ;;
        --min-matches)
            sed -i 's/min_matches: .*/min_matches: true/' "$TEMP_CONFIG"
            shift
            ;;
        --mbh)
            sed -i 's/mbh: .*/mbh: true/' "$TEMP_CONFIG"
            shift
            ;;
        --species)
            SPECIES="$2"
            # Convert comma-separated to YAML list
            YAML_SPECIES=$(echo "$SPECIES" | sed 's/,/\n  - /g')
            # This is complex for sed, better to use Python
            $PYTHON << EOF
import yaml
with open('$TEMP_CONFIG', 'r') as f:
    config = yaml.safe_load(f)
config['species'] = '$SPECIES'.split(',')
with open('$TEMP_CONFIG', 'w') as f:
    yaml.dump(config, f)
EOF
            shift 2
            ;;
        --triangle)
            TRIANGLE="$2"
            $PYTHON << EOF

with open('$TEMP_CONFIG', 'r') as f:
    config = yaml.safe_load(f)
config['gene_analysis']['triangle_species'] = '$TRIANGLE'.split(',')
with open('$TEMP_CONFIG', 'w') as f:
    yaml.dump(config, f)
EOF
            shift 2
            ;;
        --rest)
            REST="$2"
            $PYTHON << EOF

with open('$TEMP_CONFIG', 'r') as f:
    config = yaml.safe_load(f)
config['gene_analysis']['rest_species'] = '$REST'.split(',')
with open('$TEMP_CONFIG', 'w') as f:
    yaml.dump(config, f)
EOF
            shift 2
            ;;
        --cores)
            CORES="$2"
            shift 2
            ;;
        --dry-run)
            SNAKEMAKE_ARGS+=("--dry-run")
            shift
            ;;
        --forcerun)
            SNAKEMAKE_ARGS+=("--forcerun")
            shift
            ;;
        --unlock)
            SNAKEMAKE_ARGS+=("--unlock")
            shift
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

if [ "$DOWNLOAD_MODE" = true ]; then
    echo "Downloading genomes from NCBI..."
    
    if [ -z "$SPECIES_QUERIES" ]; then
        echo "Error: --species-queries required with --download"
        echo "Example: --download --species-queries 'Homo sapiens,9606' --species-labels Hsapiens,mouse"
        exit 1
    fi
    
    if [ -z "$SPECIES_LABELS" ]; then
        echo "Error: --species-labels required with --download"
        exit 1
    fi
    
    python3 workflow/scripts/ncbi_download.py \
        --species "$SPECIES_QUERIES" \
        --labels "$SPECIES_LABELS" \
        --output-dir resources \
        ${NON_INTERACTIVE}
    
    if [ $? -ne 0 ]; then
        echo "Download failed. Exiting."
        exit 1
    fi
    
    echo ""
    echo "Download complete! You can now run the pipeline with:"
    echo "  ./synteny --species $SPECIES_LABELS --parse --orthofinder --dotplots-orthofinder --ribbons-orthofinder --alg-discovery-orthofinder --ribbons-multi-orthofinder --threads 16"
    echo "or"
    echo "  ./synteny --species $SPECIES_LABELS --parse --rbh --dotplots-rbh --ribbons-rbh --alg-discovery-rbh --ribbons-multi-rbh --threads 16"
    exit 0
fi

# Run Snakemake with modified config
echo "Running Synteny Pipeline..."
echo "Using config: $TEMP_CONFIG"
echo ""

RESULTS_DIR=$($PYTHON -c "import yaml; print(yaml.safe_load(open('$TEMP_CONFIG'))['results_dir'])")

# Phase 1: run OrthoFinder/RBH alone first if requested, to avoid race
# conditions with downstream rules that read their output directly
# without a Snakemake dependency edge.
if [ "$ORTHOFINDER_REQUESTED" = true ]; then
    echo "Phase 1: Running OrthoFinder first..."
    snakemake \
        --configfile "$TEMP_CONFIG" \
        --cores "$CORES" \
        "${RESULTS_DIR}/orthofinder/done.txt" \
        "${SNAKEMAKE_ARGS[@]}"
    echo ""
    echo "Phase 1 complete. Proceeding to remaining analyses..."
    echo ""
fi

if [ "$RBH_REQUESTED" = true ]; then
    echo "Phase 1: Running RBH first..."
    snakemake \
        --configfile "$TEMP_CONFIG" \
        --cores "$CORES" \
        "${RESULTS_DIR}/rbh/MBH_Orthogroups.tsv" \
        "${SNAKEMAKE_ARGS[@]}"
    echo ""
    echo "Phase 1 complete. Proceeding to remaining analyses..."
    echo ""
fi

# Phase 2: run everything else (rule all picks up all other requested
# outputs from the config; OrthoFinder/RBH outputs already exist so
# Snakemake won't rerun them)
snakemake \
    --configfile "$TEMP_CONFIG" \
    --cores "$CORES" \
    "${SNAKEMAKE_ARGS[@]}"

# Cleanup
rm "$TEMP_CONFIG"
