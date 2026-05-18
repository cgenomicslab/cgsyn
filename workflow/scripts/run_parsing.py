import sys
sys.path.insert(0, "workflow/scripts")
from synteny import protein, proteome_filter
import pandas as pd

protein_map = protein(
    prot_file=snakemake.input.proteome,
    gff_file=snakemake.input.gff,
    output='dic'  # Get dictionary output
)

# Write TSV
df = pd.DataFrame({
    'GeneID': list(protein_map.keys()),
    'ProteinID': [v[0] for v in list(protein_map.values())],
    'Chr/Scaffold': [v[2] for v in list(protein_map.values())],
    'Strand': [v[3] for v in list(protein_map.values())],
    'Start': [v[4] for v in list(protein_map.values())],
    'End': [v[5] for v in list(protein_map.values())]
})
df.to_csv(snakemake.output.tsv, sep='\t', index=False)

# Filter Proteomes
proteome_filter(
    df = df,
    prot_file=snakemake.input.proteome,
    filtered_proteome=snakemake.output.filtered_proteome
)
