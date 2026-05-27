import sys
sys.path.insert(0, "workflow/scripts")
from synteny import *
import pickle

if snakemake.params.get('shared_ogs', False):
    raise ValueError(
        "--shared-ogs is incompatible with RBH orthology inference. "
        "RBH already produces strict 1-to-1 pairs by definition. "
        "Please use --orthofinder with --shared-ogs."
    )

species_list = snakemake.params.species_list

print(f"\n{'='*60}")
print(f"ALG DISCOVERY (RBH)")
print(f"{'='*60}")
print(f"Species: {', '.join(species_list)}")

# Load species maps from MBH pairwise TSV files
species_maps = []

for idx, species in enumerate(species_list):
    print(f"\nLoading data for {species}...")
    
    if idx == 0:
        df = df_parsing(
            f"{snakemake.params.rbh_dir}/Orthologues_{species_list[0]}/{species_list[0]}__v__{species_list[1]}.tsv",
            species_list[0],
            species_list[1]
        )
    else:
        df = df_parsing(
            f"{snakemake.params.rbh_dir}/Orthologues_{species_list[0]}/{species_list[0]}__v__{species}.tsv",
            species_list[0],
            species
        )
    
    sp_map = synteny_map_creator(df, species, snakemake.params.tsv_dir)
    species_maps.append(sp_map)
    print(f"  Loaded {len(sp_map)} chromosomes")

# Run Fisher's multi to get filtered map
print(f"\nRunning Fisher's exact tests (all pairs)...")
filtered_map = fishers_multi(
    species_maps=species_maps,
    species_names=species_list,
    alpha=snakemake.params.alpha,
    min_matches=snakemake.params.min_matches
)

# Run ALG discovery
print(f"\nRunning ALG discovery...")
alg_results = compute_algs_full_pipeline(
    species_maps=species_maps,
    species_names=species_list,
    fishers_multi_results=filtered_map,
    similarity_threshold=snakemake.params.similarity_threshold,
    min_orthologs=snakemake.params.min_orthologs,
    cluster_species=snakemake.params.cluster_species
)

# Save similarity heatmap
heatmap_path = os.path.join(
    os.path.dirname(snakemake.output.alg_results),
    "synteny_similarity_heatmap.png"
)
save_similarity_heatmap(
    similarity_matrix=alg_results['similarity_matrix'],
    species_names=species_list,
    output_path=heatmap_path
)

# Save results using pickle
print(f"\nSaving results...")
with open(snakemake.output.alg_results, 'wb') as f:
    pickle.dump({
        'alg_results': alg_results,
        'filtered_map': filtered_map,
        'species_maps': species_maps,
        'species_list': species_list
    }, f)

print(f"\n✅ ALG results saved to: {snakemake.output.alg_results}")

# Save human-readable summary
with open(snakemake.output.alg_summary, 'w') as f:
    f.write(f"ALG DISCOVERY SUMMARY (RBH)\n")
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

print(f"✅ Summary saved to: {snakemake.output.alg_summary}")
