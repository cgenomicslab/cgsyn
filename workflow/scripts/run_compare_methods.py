import sys
import os
sys.path.insert(0, "workflow/scripts")
from synteny import compare_orthology_methods

ortho_dir = snakemake.params.ortho_dir
rbh_dir = snakemake.params.rbh_dir

# Check OrthoFinder results exist
if not os.path.exists(ortho_dir):
    raise FileNotFoundError(
        f"OrthoFinder results not found at {ortho_dir}. "
        f"Please run --orthofinder first."
    )

# Check RBH results exist
mbh_path = os.path.join(rbh_dir, "MBH_Orthogroups.tsv")
if not os.path.exists(mbh_path):
    raise FileNotFoundError(
        f"RBH/MBH results not found at {mbh_path}. "
        f"Please run --rbh first."
    )

ortho_dir = open(ortho_dir).read().strip()

# Check species consistency between the two methods
import pandas as pd

mbh_df = pd.read_csv(mbh_path, sep='\t', index_col=0)
mbh_species = set(mbh_df.columns.tolist())

# Read OrthoFinder species from Orthogroups.tsv
of_orthogroups_path = os.path.join(
    os.path.dirname(ortho_dir),
    "Orthogroups", "Orthogroups.tsv"
)
if not os.path.exists(of_orthogroups_path):
    raise FileNotFoundError(
        f"OrthoFinder Orthogroups.tsv not found at {of_orthogroups_path}."
    )

of_df = pd.read_csv(of_orthogroups_path, sep='\t', index_col=0)
of_species = set(of_df.columns.tolist())

if of_species != mbh_species:
    raise ValueError(
        f"Species mismatch between OrthoFinder and RBH results!\n"
        f"OrthoFinder species: {sorted(of_species)}\n"
        f"RBH species:         {sorted(mbh_species)}\n"
        f"Please ensure both methods were run on the same set of species."
    )

# Infer species list from results instead of config
species_list = sorted(mbh_species)

compare_orthology_methods(
    orthofinder_dir=ortho_dir,
    rbh_dir=rbh_dir,
    output_path=snakemake.output.report
)
