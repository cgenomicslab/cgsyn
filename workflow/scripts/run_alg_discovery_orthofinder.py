import sys
import os
sys.path.insert(0, "workflow/scripts")
from synteny import *
import pickle
import os

results_dir = snakemake.params.results_dir
orthologues_path_file = f"{results_dir}/orthofinder/orthologues_path.txt"
if not os.path.exists(orthologues_path_file):
    raise FileNotFoundError(
        f"OrthoFinder results not found at {orthologues_path_file}. Please run --orthofinder first."
    )
ortho_dir = open(orthologues_path_file).read().strip()

species_list = snakemake.params.species_list

print(f"\n{'='*60}")
print(f"ALG DISCOVERY{'  (shared OG mode)' if snakemake.params.shared_ogs else ''}")
print(f"{'='*60}")
print(f"Species: {', '.join(species_list)}")

# Load species maps
species_maps = []
for idx, species in enumerate(species_list):
    print(f"\nLoading data for {species}...")

    if idx == 0:
        df = df_parsing(
            f"{ortho_dir}/Orthologues_{species_list[0]}/{species_list[0]}__v__{species_list[1]}.tsv",
            species_list[0], species_list[1]
        )
    else:
        df = df_parsing(
            f"{ortho_dir}/Orthologues_{species_list[0]}/{species_list[0]}__v__{species}.tsv",
            species_list[0], species
        )

    sp_map = synteny_map_creator(df, species, snakemake.params.tsv_dir)
    species_maps.append(sp_map)
    print(f"  Loaded {len(sp_map)} chromosomes")

# Run Fisher's multi to get filtered map
print(f"\nRunning Fisher's exact tests (all pairs)...")
if snakemake.params.shared_ogs:
    filtered_map = fishers_multi_shared_ogs(
        species_names=species_list,
        orthogroups_tsv_path=snakemake.params.orthogroups_tsv,
        tsv_dir=snakemake.params.tsv_dir,
        alpha=snakemake.params.alpha,
        min_matches=snakemake.params.min_matches
    )
else:
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

# Save results
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
    f.write(f"ALG DISCOVERY SUMMARY{'  (shared OG mode)' if snakemake.params.shared_ogs else ''}\n")
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
