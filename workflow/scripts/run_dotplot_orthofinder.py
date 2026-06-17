import sys
sys.path.insert(0, "workflow/scripts")
from synteny import *
import matplotlib.pyplot as plt
import os

results_dir = snakemake.params.results_dir
orthologues_path_file = f"{results_dir}/orthofinder/orthologues_path.txt"
if not os.path.exists(orthologues_path_file):
    raise FileNotFoundError(
        f"OrthoFinder results not found at {orthologues_path_file}. Please run --orthofinder first."
    )
ortho_dir = open(orthologues_path_file).read().strip()

color_palette = custom_colors_cb if snakemake.params.cb_colors else custom_colors
sp1 = snakemake.params.sp1
sp2 = snakemake.params.sp2

# Always build sp1_map and sp2_map from 1-to-1 pairs for chromosome
# length scaling in the dotplot
df = df_parsing(
    f"{ortho_dir}/Orthologues_{sp1}/{sp1}__v__{sp2}.tsv",
    sp1, sp2
)
sp1_map = synteny_map_creator(df, sp1, snakemake.params.tsv_dir)
sp2_map = synteny_map_creator(df, sp2, snakemake.params.tsv_dir)

if snakemake.params.shared_ogs:
    comparison_map = create_comparison_map_shared_ogs(
        species1_name=sp1,
        species2_name=sp2,
        orthogroups_tsv_path=snakemake.params.orthogroups_tsv,
        tsv_dir=snakemake.params.tsv_dir
    )
    counts_df, results_df, significant_pairs = fishers_shared_ogs(
        comparison_map_shared=comparison_map,
        species1_name=sp1,
        species2_name=sp2,
        orthogroups_tsv_path=snakemake.params.orthogroups_tsv,
        alpha=snakemake.params.alpha,
        min_matches=snakemake.params.min_matches,
        gene_filtering=False
    )
    comparison_map_all_pairs = create_comparison_map_shared_ogs(
        species1_name=sp1,
        species2_name=sp2,
        orthogroups_tsv_path=snakemake.params.orthogroups_tsv,
        tsv_dir=snakemake.params.tsv_dir,
        all_pairs=True
    )
    fig, ax = plot_synteny_dotplot(
        comparison_map_all_pairs, sp1_map, sp2_map, significant_pairs,
        species1=sp1,
        species2=sp2,
        dot_size=snakemake.params.dot_size,
        dot_alpha=snakemake.params.dot_alpha,
        color_nonsignificant=snakemake.params.color_nonsignificant,
        color_palette=color_palette
    )
else:
    comparison_map = create_comparison_map(sp1_map, sp2_map, sp1, sp2)
    counts_df, results_df, significant_pairs = fishers(
        comparison_map,
        alpha=snakemake.params.alpha,
        min_matches=snakemake.params.min_matches
    )

    fig, ax = plot_synteny_dotplot(
        comparison_map, sp1_map, sp2_map, significant_pairs,
        species1=sp1,
        species2=sp2,
        dot_size=snakemake.params.dot_size,
        dot_alpha=snakemake.params.dot_alpha,
        color_nonsignificant=snakemake.params.color_nonsignificant,
        color_palette=color_palette
    )

plt.savefig(snakemake.output.plot, dpi=300, bbox_inches='tight')
plt.close()
