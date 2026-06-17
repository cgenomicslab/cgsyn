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
    filtered_map = fishers_shared_ogs(
        comparison_map_shared=comparison_map,
        species1_name=sp1,
        species2_name=sp2,
        orthogroups_tsv_path=snakemake.params.orthogroups_tsv,
        alpha=snakemake.params.alpha,
        min_matches=snakemake.params.min_matches,
        gene_filtering=True
    )
    # Build all-pairs version of filtered_map for visualization
    comparison_map_all_pairs = create_comparison_map_shared_ogs(
        species1_name=sp1,
        species2_name=sp2,
        orthogroups_tsv_path=snakemake.params.orthogroups_tsv,
        tsv_dir=snakemake.params.tsv_dir,
        all_pairs=True
    )
    # Keep only pairs whose chromosomes are in filtered_map
    significant_chroms = {
        (v[0][1], v[1][1]) for v in filtered_map.values()
    }
    filtered_map_all_pairs = {
        k: v for k, v in comparison_map_all_pairs.items()
        if (v[0][1], v[1][1]) in significant_chroms
    }
    fig, ax = plot_synteny_ribbons(
        filtered_map_all_pairs, sp1_map, sp2_map,
        species1=sp1,
        species2=sp2,
        ribbon_alpha=snakemake.params.ribbon_alpha,
        curve_style=snakemake.params.curve_style,
        color_palette=color_palette
    )
else:
    comparison_map = create_comparison_map(sp1_map, sp2_map, sp1, sp2)
    filtered_map = fishers(
        comparison_map,
        alpha=snakemake.params.alpha,
        min_matches=snakemake.params.min_matches,
        gene_filtering=True
    )
    print(f"Normal mode - filtered_map size: {len(filtered_map)}")
    print(f"Normal mode SP1 chroms: {sorted(set(v[0][1] for v in filtered_map.values()))}")
    print(f"Normal mode SP2 chroms: {sorted(set(v[1][1] for v in filtered_map.values()))}")

    fig, ax = plot_synteny_ribbons(
        filtered_map, sp1_map, sp2_map,
        species1=sp1,
        species2=sp2,
        ribbon_alpha=snakemake.params.ribbon_alpha,
        curve_style=snakemake.params.curve_style,
        color_palette=color_palette
    )

plt.savefig(snakemake.output.plot, dpi=300, bbox_inches='tight')
plt.close()
