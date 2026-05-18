import sys
sys.path.insert(0, "workflow/scripts")
from synteny import *
import matplotlib.pyplot as plt

# Parse OrthoFinder data
df = df_parsing(
    f"{snakemake.params.ortho_dir}/Orthologues_{snakemake.params.sp1}/{snakemake.params.sp1}__v__{snakemake.params.sp2}.tsv",
    snakemake.params.sp1,
    snakemake.params.sp2
)

sp1_map = synteny_map_creator(df, snakemake.params.sp1, snakemake.params.tsv_dir)
sp2_map = synteny_map_creator(df, snakemake.params.sp2, snakemake.params.tsv_dir)
comparison_map = create_comparison_map(sp1_map, sp2_map, snakemake.params.sp1, snakemake.params.sp2)

# Run Fisher's test
counts_df, results_df, significant_pairs = fishers(
    comparison_map,
    alpha=snakemake.params.alpha,
    min_matches=snakemake.params.min_matches
)

# Create plot
fig, ax = plot_synteny_dotplot(
    comparison_map, sp1_map, sp2_map, significant_pairs,
    species1=snakemake.params.sp1,
    species2=snakemake.params.sp2,
    dot_size=snakemake.params.dot_size,
    dot_alpha=snakemake.params.dot_alpha
)

plt.savefig(snakemake.output.plot, dpi=300, bbox_inches='tight')
plt.close()
