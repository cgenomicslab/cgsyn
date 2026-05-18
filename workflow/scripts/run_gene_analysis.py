import sys
sys.path.insert(0, "workflow/scripts")
from synteny import simakov_gene_families
import pandas as pd

families, details = simakov_gene_families(
    triangle_species=snakemake.params.triangle,
    rest_species=snakemake.params.rest,
    rbh_dir=snakemake.params.rbh_dir,
    tsv_dir=snakemake.params.tsv_dir,
    mbh=snakemake.params.mbh
)

# Save families as TSV
all_species = snakemake.params.triangle + snakemake.params.rest
df = pd.DataFrame(families, columns=all_species)
df.to_csv(snakemake.output.families, sep='\t', index=False)

# Save details
with open(snakemake.output.details, 'w') as f:
    for i, detail in enumerate(details, 1):
        f.write(f"Family {i}:\n")
        f.write(f"  Proteins: {detail['proteins']}\n")
        f.write(f"  Chromosomes: {detail['chromosomes']}\n\n")
