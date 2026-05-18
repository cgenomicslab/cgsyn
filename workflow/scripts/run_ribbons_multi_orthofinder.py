import sys
import os
sys.path.insert(0, "workflow/scripts")
from synteny import *
import matplotlib.pyplot as plt
import pickle

species_list = snakemake.params.species_list

print(f"\n{'='*60}")
print(f"MULTI-SPECIES RIBBON PLOT")
print(f"{'='*60}")
print(f"Species: {', '.join(species_list)}")
print(f"Total species: {len(species_list)}")

# Load species maps
species_maps = []

for idx, species in enumerate(species_list):
    print(f"\nLoading data for {species}...")
    
    if idx == 0:
        df = df_parsing(
            f"{snakemake.params.ortho_dir}/Orthologues_{species_list[0]}/{species_list[0]}__v__{species_list[1]}.tsv",
            species_list[0],
            species_list[1]
        )
        sp_map = synteny_map_creator(df, species, snakemake.params.tsv_dir)
    else:
        df = df_parsing(
            f"{snakemake.params.ortho_dir}/Orthologues_{species_list[0]}/{species_list[0]}__v__{species}.tsv",
            species_list[0],
            species
        )
        sp_map = synteny_map_creator(df, species, snakemake.params.tsv_dir)
    
    species_maps.append(sp_map)
    print(f"  Loaded {len(sp_map)} chromosomes")

# Run Fisher's multi to get filtered map
print(f"\nRunning Fisher's exact tests...")
filtered_map = fishers_multi(
    species_maps=species_maps,
    species_names=species_list,
    alpha=snakemake.params.alpha,
    min_matches=snakemake.params.min_matches
)

# Run ALG discovery if flag is set
alg_results = None
if snakemake.params.alg_discovery:
    print(f"\nRunning ALG discovery...")
    alg_results = compute_algs_full_pipeline(
        species_maps=species_maps,
        species_names=species_list,
        fishers_multi_results=filtered_map,
        similarity_threshold=snakemake.params.similarity_threshold,
        min_orthologs=snakemake.params.min_orthologs
    )
    
    # Save ALG results as pickle for future use
    os.makedirs(os.path.dirname(snakemake.params.alg_results_path), exist_ok=True)
    with open(snakemake.params.alg_results_path, 'wb') as f:
        pickle.dump({
            'alg_results': alg_results,
            'filtered_map': filtered_map,
            'species_maps': species_maps,
            'species_list': species_list
        }, f)
    print(f"✅ ALG results saved to: {snakemake.params.alg_results_path}")
    
    # Save human-readable summary
    with open(snakemake.params.alg_summary_path, 'w') as f:
        f.write(f"ALG DISCOVERY SUMMARY\n")
        f.write(f"{'='*60}\n")
        f.write(f"Species: {', '.join(species_list)}\n\n")
        
        clusters_info = alg_results['clusters_info']
        alg_assignments = alg_results['alg_assignments']
        
        f.write(f"CLUSTERS:\n")
        for cluster_id, cluster_species in sorted(clusters_info.items()):
            f.write(f"  Cluster {cluster_id}: {', '.join(cluster_species)}\n")
        
        f.write(f"\nALG ASSIGNMENTS:\n")
        for cluster_id, cluster_species in sorted(clusters_info.items()):
            f.write(f"\n  Cluster {cluster_id}:\n")
            
            cluster_algs = {}
            for sp in cluster_species:
                if sp in alg_assignments:
                    for chrom, alg_list in alg_assignments[sp].items():
                        if not isinstance(alg_list, list):
                            alg_list = [alg_list]
                        for alg_id in alg_list:
                            if alg_id not in cluster_algs:
                                cluster_algs[alg_id] = {}
                            cluster_algs[alg_id][sp] = chrom
            
            for alg_id in sorted(cluster_algs.keys()):
                f.write(f"    {alg_id}:\n")
                for sp in cluster_species:
                    if sp in cluster_algs[alg_id]:
                        f.write(f"      {sp}: {cluster_algs[alg_id][sp]}\n")
    
    print(f"✅ ALG summary saved to: {snakemake.params.alg_summary_path}")
else:
    print(f"\nSkipping ALG discovery (use --alg-discovery flag to enable)")

# Plot
print(f"\nGenerating multi-species ribbon plot...")
fig, ax = plot_synteny_ribbons_multi(
    filtered_map=filtered_map,
    species_maps=species_maps,
    species_names=species_list,
    alg_results=alg_results,  # None if no ALG discovery = original coloring
    ribbon_alpha=snakemake.params.ribbon_alpha,
    figsize=tuple(snakemake.params.figsize),
    curve_style=snakemake.params.curve_style
)

plt.savefig(snakemake.output.plot, dpi=300, bbox_inches='tight')
plt.close()
print(f"\n✅ Multi-species ribbon plot saved to: {snakemake.output.plot}")
