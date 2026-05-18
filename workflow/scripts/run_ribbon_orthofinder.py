import sys
sys.path.insert(0, "workflow/scripts")
from synteny import *
import matplotlib.pyplot as plt

# Parse Orthofinder data
df = df_parsing(
    f"{snakemake.params.ortho_dir}/Orthologues_{snakemake.params.sp1}/{snakemake.params.sp1}__v__{snakemake.params.sp2}.tsv",
    snakemake.params.sp1,
    snakemake.params.sp2
)

# Create Maps
sp1_map = synteny_map_creator(df, snakemake.params.sp1, snakemake.params.tsv_dir)
sp2_map = synteny_map_creator(df, snakemake.params.sp2, snakemake.params.tsv_dir)
comparison_map = create_comparison_map(sp1_map, sp2_map, snakemake.params.sp1, snakemake.params.sp2)

# Filter
filtered_map = fishers(
    comparison_map, 
    alpha=snakemake.params.alpha,
    min_matches=snakemake.params.min_matches,
    gene_filtering=True
)

# Plot
fig, ax = plot_synteny_ribbons(
    filtered_map, sp1_map, sp2_map,
    species1=snakemake.params.sp1,
    species2=snakemake.params.sp2,
    ribbon_alpha=snakemake.params.ribbon_alpha,
    curve_style=snakemake.params.curve_style
)

plt.savefig(snakemake.output.plot, dpi=300, bbox_inches='tight')
plt.close()
