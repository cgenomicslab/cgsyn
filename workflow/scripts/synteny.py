# BEFORE ORTHOLOGY PREDICTION

## Importing Modules

import gzip
import re
import pandas as pd # pyright: ignore[reportMissingModuleSource]
import os
import traceback
import subprocess
from pathlib import Path
import itertools
import sys
import shutil
from collections import defaultdict
from scipy.stats import fisher_exact # pyright: ignore[reportMissingImports]
from scipy.cluster.hierarchy import linkage, fcluster # type: ignore
from scipy.spatial.distance import squareform # type: ignore
import numpy as np
import matplotlib.pyplot as plt # pyright: ignore[reportMissingModuleSource]

## Proteome/Annotation File Parser

def protein(prot_file, gff_file, output='tsv',):
    '''
    Parse the proteome file and gff file for each species, filter
    through protein isoforms and keep the longest isoform for each
    gene. Output is a .tsv file containing the GeneID, ProteinID, 
    chr/scaffold where it's located, strand it's on and its coords.
    
    Parameters:
    -----------
    prot_file: .faa file
        Fasta proteome file from NCBI
    gff_file: .gff file
        Text functional annotation file from NCBI
    output: str, default = 'tsv'
        string specifying the type of output (options: 'dic', 'df', 'tsv')

    Returns:
    --------
    If output = dic:
        protein to gene map in dictionary format - DEFAULT inside Snakemake
    
    If output is df:
        position of each gene on chromosomes in pd.Dataframe format
    
    If output is tsv:
        position of each gene on chromosomes in .tsv format (file creation) - DEFAULT if used outside Snakemake
    '''

    with gzip.open(gff_file, 'rt') as f:
        g2p = {}
        for line in f.readlines():
            if line.startswith("#") or not line.strip():
                continue

            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            feature_type = parts[2]
            start, end, strand = parts[3], parts[4], parts[6]
            attributes = parts[8]

            if 'CDS' in feature_type:
                if 'GeneID:' in line or 'locus_tag=' in line:
                    g_match = re.search(r'GeneID:(\d+).*?Name=([\w\.]+)', line)
                    l_match = re.search(r'Name=([\w\.]+).*?locus_tag=([^;]+)', line)
                    if g_match:
                        gene_id, protein_id = g_match.groups()
                        g2p.setdefault(gene_id, []).append(protein_id)
                    else:
                        if l_match:
                            protein_id, gene_id = l_match.groups()
                            g2p.setdefault(gene_id, []).append(protein_id)
            
            
    for k, v in g2p.items():
        g2p[k] = list(set(v))
    
    with gzip.open(prot_file, 'rt') as f:
        seq = ''
        seqs = []
        lines = f.readlines()
        headers = [line.strip('\n') for line in lines if '>' in line]
    
        for j, line in enumerate(lines):
            line = line.strip('\n')
            if '>'not in line:
                try:
                    if '>' not in lines[j+1]:
                        seq+=line
                    else:
                        seq+=line
                        seqs.append(seq)
                        seq = ''
                except IndexError:
                    seq+=line
                    seqs.append(seq)
                    seq = ''
                    break

    p2seq = {i.split()[0][1:]:j 
     for i, j in zip(headers, seqs)}
    
    p2g_map = {}

    for k, v in g2p.items():
        longest_seq = ''
        longest_id = ''
        for i in v:
            temp_seq = p2seq[i]
            if len(temp_seq) > len(longest_seq):
                longest_seq = temp_seq
                longest_id = i
        p2g_map[k] = [longest_id, longest_seq]

    current_chr = None
    with gzip.open(gff_file, 'rt') as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            feature_type = parts[2]
            start, end, strand = parts[3], parts[4], parts[6]
            attributes = parts[8]

            if "region" in feature_type:
                chr_match = re.search(r"chromosome=([A-Za-z0-9]+)", attributes, re.IGNORECASE)
                if chr_match:
                    if chr_match.group(1) == "Unknown":
                        current_chr = parts[0]
                    else:
                        current_chr = chr_match.group(1)
                else:
                    name_match = re.search(r"Name=([^;]+)", attributes)
                    current_chr = name_match.group(1) if name_match else parts[0] #fallback to seqID, usually for contig-level assemblies
                continue

            if "gene" in feature_type and current_chr is not None:
                gene_match = re.search(r"GeneID:(\d+)", attributes)
                id_match = re.search(r"locus_tag=([^;]+)", attributes)
                if gene_match:
                    gene_id = gene_match.group(1)
                    if gene_id in p2g_map:
                        p2g_map[gene_id].extend([current_chr, strand, start, end])
                else:
                    if id_match:
                        gene_id = id_match.group(1)
                        if gene_id in p2g_map:
                            p2g_map[gene_id].extend([current_chr, strand, start, end])

    df = pd.DataFrame({
        'GeneID' : list(p2g_map.keys()),
        'ProteinID' : [i[0] for i in list(p2g_map.values())],
        'Chr/Scaffold' : [i[2] for i in list(p2g_map.values())],
        'Strand' : [i[3] for i in list(p2g_map.values())],
        'Start' : [i[4] for i in list(p2g_map.values())],
        'End' : [i[5] for i in list(p2g_map.values())]
        })
    
    if output == 'dic':
        return p2g_map
    elif output == 'df':
        return df
    elif output == 'tsv':
        df.to_csv(f"../../intermediates/tsv/{''.join(gff_file.split('.')[0])}.tsv", sep="\t", index=False)
        return None
    
# Create Filtered Proteomes Only With The Longest Isoforms For Orthofinder

def proteome_filter(df, prot_file, filtered_proteome):
    '''
    Parse the dataframe and proteome file from previous function and
    create filtered proteome with only the longest isoform of each gene
    
    Parameters:
    -----------
    df: pd.DataFrame
        dataframe output from previous function
    prot_file: .faa file
        Fasta proteome file from NCBI
    filtered_proteome: .faa file
        Fasta filtered proteome output
    '''

    # Get allowed protein IDs
    allowed_ids = set(df['ProteinID'].tolist())
    
    # Determine if input is gzipped
    is_gzipped = prot_file.endswith('.gz')
    open_func = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'
    
    # Single-pass filtering: read gzipped input, write filtered output
    keep = False
    current_header = None
    current_seq = []
    
    with open_func(prot_file, mode) as f_in, open(filtered_proteome, 'w') as f_out:
        for line in f_in:
            line = line.strip()
            
            if line.startswith('>'):
                # Write previous sequence if it was kept
                if keep and current_header:
                    f_out.write(current_header + '\n')
                    f_out.write(''.join(current_seq) + '\n')
                
                # Check if we should keep this sequence
                # Extract protein ID (everything after '>' up to first whitespace)
                protein_id = line[1:].split()[0]
                
                if protein_id in allowed_ids:
                    keep = True
                    current_header = line
                    current_seq = []
                else:
                    keep = False
                    current_header = None
                    current_seq = []
            else:
                # Accumulate sequence if we're keeping it
                if keep:
                    current_seq.append(line)
        
        # Write last sequence if kept
        if keep and current_header:
            f_out.write(current_header + '\n')
            f_out.write(''.join(current_seq) + '\n')

# ORTHOLOGY PREDICTION

## With RBH (Reciprocal Best BLASTP Hits)

def run_alignment(query_fasta, subject_fasta, output_file, num_threads=12, aligner='diamond'):  
    """
    Run BLASTP or DIAMOND between query and subject proteomes.
    
    Parameters:
    -----------
    query_fasta : str
        Path to query protein fasta file
    subject_fasta : str
        Path to subject protein fasta file
    output_file : str
        Path to output alignment results
    num_threads : int
        Number of threads
    aligner : str
        Choice of alignment tool: 'diamond' or 'blast' (default: 'diamond')
    """
    
    if aligner.lower() == 'diamond':
        # Create DIAMOND database for subject
        db_name = f"{subject_fasta}.dmnd"
        
        print(f"Creating DIAMOND database for {subject_fasta}...")
        makedb_cmd = [
            "diamond", "makedb",
            "--in", subject_fasta,
            "--db", db_name.replace('.dmnd', '')  # Diamond adds .dmnd automatically
        ]
        subprocess.run(makedb_cmd, check=True)
        
        # Run DIAMOND
        print(f"Running DIAMOND: {query_fasta} vs {subject_fasta}...")
        diamond_cmd = [
            "diamond", "blastp",
            "--query", query_fasta,
            "--db", db_name.replace('.dmnd', ''),
            "--out", output_file,
            "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch", 
                        "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore",
            "--max-target-seqs", "1",  # Keep only best hit
            "--threads", str(num_threads),
            "--evalue", "1e-5",
            "--sensitive"  # More sensitive mode (slower but better)
        ]
        subprocess.run(diamond_cmd, check=True)
        
    elif aligner.lower() == 'blast':
        # Create BLAST database for subject
        db_name = f"{subject_fasta}.db"
        
        print(f"Creating BLAST database for {subject_fasta}...")
        makedb_cmd = [
            "makeblastdb",
            "-in", subject_fasta,
            "-dbtype", "prot",
            "-out", db_name
        ]
        subprocess.run(makedb_cmd, check=True)
        
        # Run BLASTP
        print(f"Running BLASTP: {query_fasta} vs {subject_fasta}...")
        blastp_cmd = [
            "blastp",
            "-query", query_fasta,
            "-db", db_name,
            "-out", output_file,
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
            "-max_target_seqs", "1",  # Keep only best hit
            "-num_threads", str(num_threads),
            "-evalue", "1e-5"
        ]
        subprocess.run(blastp_cmd, check=True)
    
    else:
        raise ValueError(f"Unknown aligner: {aligner}. Must be 'diamond' or 'blast'")
    
    print(f"Alignment results saved to {output_file}")


def parse_alignment_results(alignment_file):  
    """
    Parse BLAST/DIAMOND output and extract best hits.
    
    Parameters:
    -----------
    alignment_file : str
        Path to BLAST/DIAMOND output file
    
    Returns:
    --------
    dict
        Dictionary mapping query IDs to their best hit subject IDs
    """
    
    best_hits = {}
    
    if not os.path.exists(alignment_file) or os.path.getsize(alignment_file) == 0:
        print(f"Warning: {alignment_file} is empty or doesn't exist")
        return best_hits
    
    with open(alignment_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            query = fields[0]
            subject = fields[1]
            bitscore = float(fields[11])
            
            # Keep best hit (highest bitscore) for each query
            if query not in best_hits:
                best_hits[query] = (subject, bitscore)
            else:
                if bitscore > best_hits[query][1]:
                    best_hits[query] = (subject, bitscore)
    
    # Convert to simple dict (remove bitscore)
    best_hits = {query: subject for query, (subject, score) in best_hits.items()}
    
    return best_hits


def find_reciprocal_best_hits(proteome1, proteome2, sp1, sp2, output_dir, num_threads=12,  aligner='diamond'):   
    """
    Find reciprocal best BLASTP hits between two proteomes.
    
    Parameters:
    -----------
    proteome1: str
        Path to proteome fasta file for species 1
    proteome2: str
        Path to proteome fasta file for species 2
    sp1 : str
        Name of species 1 (for output headers)
    sp2 : str
        Name of species 2 (for output headers)
    output_dir : str
        Directory to save output files
    num_threads : int
        Number of threads for aligner
    aligner : str
        Choice of alignment tool: 'diamond' (default) or 'blast'
        DIAMOND is ~100-1000x faster than BLAST
    
    Returns:
    --------
    pd.DataFrame
        DataFrame with reciprocal best hits
    str
        Path to output TSV file
    """
    
    # Define output file names
    alignment_forward = f"{output_dir}/{sp1}_vs_{sp2}.{aligner}"
    alignment_reverse = f"{output_dir}/{sp2}_vs_{sp1}.{aligner}"
    rbh_output = f"{output_dir}/{sp1}__RBH__{sp2}.tsv"
    
    aligner_name = "DIAMOND" if aligner.lower() == 'diamond' else "BLASTP"
    
    print("="*80)
    print(f"RECIPROCAL BEST HIT ANALYSIS ({aligner_name}): {sp1} <-> {sp2}")
    print("="*80)
    
    # Step 1: Forward alignment (species1 vs species2)
    print(f"\n[1/4] Forward {aligner_name}: {sp1} -> {sp2}")
    run_alignment(proteome1, proteome2, alignment_forward, num_threads, aligner)
    
    # Step 2: Reverse alignment (species2 vs species1)
    print(f"\n[2/4] Reverse {aligner_name}: {sp2} -> {sp1}")
    run_alignment(proteome2, proteome1, alignment_reverse, num_threads, aligner)
    
    # Step 3: Parse alignment results
    print(f"\n[3/4] Parsing {aligner_name} results...")
    forward_hits = parse_alignment_results(alignment_forward)
    reverse_hits = parse_alignment_results(alignment_reverse)
    
    print(f"  Forward hits: {len(forward_hits)}")
    print(f"  Reverse hits: {len(reverse_hits)}")
    
    # Step 4: Find reciprocal best hits
    print(f"\n[4/4] Identifying reciprocal best hits...")
    rbh_pairs = []
    
    for sp1_gene, sp2_gene in forward_hits.items():
        # Check if reverse hit exists and points back to the same gene
        if sp2_gene in reverse_hits:
            if reverse_hits[sp2_gene] == sp1_gene:
                rbh_pairs.append((sp1_gene, sp2_gene))
    
    print(f"  Found {len(rbh_pairs)} reciprocal best hits")
    
    # Step 5: Create DataFrame and save
    rbh_df = pd.DataFrame(rbh_pairs, columns=[sp1, sp2])
    rbh_df.to_csv(rbh_output, sep='\t', index=False)
    
    print(f"\n{'='*80}")
    print(f"RESULTS SAVED TO: {rbh_output}")
    print(f"{'='*80}\n")
    
    return rbh_df, rbh_output

def find_mutual_best_hits_multi(species_list, proteome_dir, output_dir, 
                                  num_threads=12, aligner='diamond'):
    '''
    Find Mutual Best Hits across multiple species and create an OrthoFinder-like
    orthogroups matrix and pairwise TSV files that can be used for fishers_multi,
    ALG discovery and multi-ribbon plots.
    
    Parameters:
    -----------
    species_list : list of str
        List of species IDs
    proteome_dir : str
        Directory containing filtered proteome files ({species}.faa)
    output_dir : str
        Directory to save output files
    num_threads : int
        Number of threads for aligner
    aligner : str
        'diamond' or 'blast'
    
    Returns:
    --------
    orthogroups_df : pd.DataFrame
        OrthoFinder-like matrix with orthogroup IDs as index and species as columns
    orthogroups_path : str
        Path to saved orthogroups TSV file
    pairwise_tsv_paths : dict
        {(sp1, sp2): path_to_pairwise_tsv}
    '''
    
    os.makedirs(output_dir, exist_ok=True)
    
    print("="*80)
    print(f"MULTI-SPECIES MUTUAL BEST HIT ANALYSIS")
    print(f"Species: {', '.join(species_list)}")
    print(f"Aligner: {aligner.upper()}")
    print("="*80)
    
    # ===== STEP 1: Run RBH for ALL pairs =====
    all_rbh = {}
    
    for sp1, sp2 in itertools.combinations(species_list, 2):
        print(f"\nRunning RBH: {sp1} <-> {sp2}")
        
        proteome1 = f"{proteome_dir}/{sp1}.faa"
        proteome2 = f"{proteome_dir}/{sp2}.faa"
        
        rbh_df, rbh_path = find_reciprocal_best_hits(
            proteome1=proteome1,
            proteome2=proteome2,
            sp1=sp1,
            sp2=sp2,
            output_dir=output_dir,
            num_threads=num_threads,
            aligner=aligner
        )
        
        # Store as lookup dicts in both directions
        all_rbh[(sp1, sp2)] = dict(zip(rbh_df[sp1], rbh_df[sp2]))
        all_rbh[(sp2, sp1)] = dict(zip(rbh_df[sp2], rbh_df[sp1]))
        
        print(f"  {sp1} <-> {sp2}: {len(rbh_df)} RBH pairs")
    
    # ===== STEP 2: Build gene families =====
    print(f"\n{'='*80}")
    print("BUILDING GENE FAMILIES FROM MBH")
    print("="*80)
    
    # Collect all genes for first species
    all_genes_first = set(all_rbh.get((species_list[0], species_list[1]), {}).keys())
    
    gene_families = []
    used_genes = {sp: set() for sp in species_list}
    
    sp_first = species_list[0]
    
    for seed_gene in sorted(all_genes_first):
        if seed_gene in used_genes[sp_first]:
            continue
        
        # Try to build a family starting from this seed gene
        family = {sp_first: seed_gene}
        valid = True
        
        # For each other species, find the RBH with the first species
        for sp_other in species_list[1:]:
            if (sp_first, sp_other) not in all_rbh:
                valid = False
                break
            
            if seed_gene not in all_rbh[(sp_first, sp_other)]:
                valid = False
                break
            
            candidate = all_rbh[(sp_first, sp_other)][seed_gene]
            
            if candidate in used_genes[sp_other]:
                valid = False
                break
            
            family[sp_other] = candidate
        
        if not valid or len(family) < 2:
            continue
        
        # Verify mutual best hits between ALL pairs (not just vs first species)
        all_pairs_valid = True
        for sp1, sp2 in itertools.combinations(species_list, 2):
            if sp1 not in family or sp2 not in family:
                continue
            
            gene1 = family[sp1]
            gene2 = family[sp2]
            
            # Check sp1 -> sp2
            if (sp1, sp2) in all_rbh:
                if all_rbh[(sp1, sp2)].get(gene1) != gene2:
                    all_pairs_valid = False
                    break
            
            # Check sp2 -> sp1
            if (sp2, sp1) in all_rbh:
                if all_rbh[(sp2, sp1)].get(gene2) != gene1:
                    all_pairs_valid = False
                    break
        
        if not all_pairs_valid:
            continue
        
        # Valid family - add it and mark genes as used
        gene_families.append(family)
        for sp, gene in family.items():
            used_genes[sp].add(gene)
    
    print(f"\nFound {len(gene_families)} mutual best hit gene families")
    
    # ===== STEP 3: Create OrthoFinder-like orthogroups matrix =====
    print(f"\n{'='*80}")
    print("CREATING ORTHOGROUPS MATRIX")
    print("="*80)
    
    rows = []
    for idx, family in enumerate(gene_families):
        og_id = f"MBH{idx:07d}"
        row = {'Orthogroup': og_id}
        for sp in species_list:
            row[sp] = family.get(sp, '')
        rows.append(row)
    
    orthogroups_df = pd.DataFrame(rows)
    orthogroups_df = orthogroups_df.set_index('Orthogroup')
    
    orthogroups_path = f"{output_dir}/MBH_Orthogroups.tsv"
    orthogroups_df.reset_index().to_csv(orthogroups_path, sep='\t', index=False)
    
    print(f"Orthogroups matrix saved to: {orthogroups_path}")
    print(f"Total gene families: {len(orthogroups_df)}")
    print(f"Shape: {orthogroups_df.shape}")
    
    # ===== STEP 4: Create pairwise 1-to-1 TSV files =====
    # Mirror OrthoFinder's Orthologues_{sp1}/{sp1}__v__{sp2}.tsv format
    # so all downstream functions work without modification
    print(f"\n{'='*80}")
    print("CREATING PAIRWISE ORTHOLOG TSV FILES")
    print("="*80)
    
    pairwise_tsv_paths = {}
    
    for sp1, sp2 in itertools.combinations(species_list, 2):
        pair_dir = f"{output_dir}/Orthologues_{sp1}"
        os.makedirs(pair_dir, exist_ok=True)
        pair_path = f"{pair_dir}/{sp1}__v__{sp2}.tsv"
        
        pair_rows = []
        for idx, family in enumerate(gene_families):
            if sp1 in family and sp2 in family:
                og_id = f"MBH{idx:07d}"
                pair_rows.append({
                    'Orthogroup': og_id,
                    sp1: family[sp1],
                    sp2: family[sp2]
                })
        
        pair_df = pd.DataFrame(pair_rows)
        pair_df.to_csv(pair_path, sep='\t', index=False)
        pairwise_tsv_paths[(sp1, sp2)] = pair_path
        
        print(f"  {sp1} vs {sp2}: {len(pair_rows)} pairs → {pair_path}")
    
    return orthogroups_df, orthogroups_path, pairwise_tsv_paths

## With OrthoFinder

def run_orthofinder(filtered_prot_dir, output_dir, num_threads=12, aligner='diamond', extra_args=None):   
    """
    Run OrthoFinder on a directory of proteome files.
    
    Parameters:
    -----------
    filtered_prot_dir : str
        Path to directory containing protein fasta files (one per species)
        Files should be named like: Species1.fasta, Species2.fasta, etc.
    num_threads : int
        Number of threads/CPUs to use (default: 4)
    aligner : str
        Alignment tool: 'diamond' (default, fast) or 'blast' (slower, sensitive)
    output_dir : str, optional
        Custom output directory. If None, OrthoFinder creates one automatically
    extra_args : list, optional
        Additional OrthoFinder arguments as a list
        Example: ['-M', 'msa', '-T', 'iqtree']
    
    Returns:
    --------
    str
        Path to OrthoFinder results directory
    
    Examples:
    ---------
    # Basic usage with DIAMOND (fast)
    >>> results_dir = run_orthofinder('/path/to/proteomes', num_threads=16)
    
    # Using BLASTP (slower but sensitive)
    >>> results_dir = run_orthofinder('/path/to/proteomes', 
                                      num_threads=16, 
                                      aligner='blast')
    
    # With custom output directory
    >>> results_dir = run_orthofinder('/path/to/proteomes',
                                      num_threads=16,
                                      output_dir='./My_OrthoFinder_Results')
    
    # With additional arguments (gene trees + species tree)
    >>> results_dir = run_orthofinder('/path/to/proteomes',
                                      num_threads=16,
                                      extra_args=['-M', 'msa', '-T', 'iqtree'])
    """
    
    # Validate proteome directory exists
    if not os.path.exists(filtered_prot_dir):
        raise FileNotFoundError(f"Proteome directory not found: {filtered_prot_dir}")
    
    # Check if directory contains fasta files
    fasta_files = list(Path(filtered_prot_dir).glob('*.fa*'))
    if not fasta_files:
        raise ValueError(f"No fasta files found in {filtered_prot_dir}")
    
    print("="*80)
    print("RUNNING ORTHOFINDER")
    print("="*80)
    print(f"Proteome directory: {filtered_prot_dir}")
    print(f"Number of species: {len(fasta_files)}")
    print(f"Aligner: {aligner.upper()}")
    print(f"Threads: {num_threads}")
    print("="*80 + "\n")
    
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    # Build OrthoFinder command
    orthofinder_cmd = [
        "orthofinder",
        "-f", filtered_prot_dir,
        "-t", str(num_threads),
    ]
    
    # Add aligner selection
    if aligner.lower() == 'diamond':
        orthofinder_cmd.extend(["-S", "diamond"])
    elif aligner.lower() == 'blast':
        orthofinder_cmd.extend(["-S", "blast"])
    else:
        raise ValueError(f"Unknown aligner: {aligner}. Must be 'diamond' or 'blast'")
    
    # Add output directory if specified
    if output_dir:
        orthofinder_cmd.extend(["-o", output_dir])
    
    # Add any extra arguments
    if extra_args:
        orthofinder_cmd.extend(extra_args)
    
    print(f"Running command: {' '.join(orthofinder_cmd)}\n")
    
    # Run OrthoFinder
    try:
        result = subprocess.run(
            orthofinder_cmd,
            check=True,
            stdout=sys.stdout,
            stderr=sys.stderr,
            text=True
        )
        
        # Print output
        print(result.stdout)
        
        # Find the results directory
        if output_dir:
            results_base = output_dir
        else:
            results_base = f"{filtered_prot_dir}/OrthoFinder"
        
        # Find the most recent Results_XXX directory
        results_dirs = sorted(Path(results_base).glob("Results_*"))
        if results_dirs:
            latest_results = str(results_dirs[-1])
            print("\n" + "="*80)
            print("ORTHOFINDER COMPLETED SUCCESSFULLY")
            print("="*80)
            print(f"Results directory: {latest_results}")
            print("="*80 + "\n")
            orthologues_dir = f"{latest_results}/Orthologues"
            return orthologues_dir
        else:
            raise RuntimeError(
                "OrthoFinder did not produce any Results_* directory. "
                "Check logs for errors."
            )
            
    except subprocess.CalledProcessError as e:
        print(f"\nERROR: OrthoFinder failed with exit code {e.returncode}")
        print(f"STDERR: {e.stderr}")
        raise


# AFTER ORTHOLOGY PREDICTION

## Orthofinder Results (1-1 Orthologues) Parsing

def df_parsing(o_dir, species1, species2): 
    '''
    Parse the .tsv files generated by Orthofinder listing the 1 to 1
    orthologue matches between species pairs, filter out non-inferred
    orthologous pairs (multiple genes in each "pair") and output a filtered
    pandas dataframe.
    
    Parameters:
    -----------
    o_dir: str
        "Orthologues" directory from Orthofinder results
    species1: str
        species 1 ID
    species2: str
        species 2 ID

    Returns:
    --------
    df_sc: pd.Dataframe
        filtered 1 to 1 orthologue pairs in dataframe format
    '''

    df = pd.read_csv(o_dir, sep="\t")
    
    df_sc = df[
        ~(
            df[species1].str.contains(",", na=False) |
            df[species2].str.contains(",", na=False)
        )
    ]
    
    df_sc["Orthogroup"] = (
        df_sc["Orthogroup"]
        + "."
        + (df_sc.groupby("Orthogroup").cumcount() + 1).astype(str)
    )

    return df_sc

## Synteny Map Creation (with Orthofinder)

def synteny_map_creator(df, species, tsv_dir):
    '''
    Combine the orthology data with the chromosome mapping data in order to
    create a physical map of the chromosome, but annotated with OrthogroupID 
    and ProteinID instead of GeneID. Output is said map in dict format.
    
    Parameters:
    -----------
    df: pd.Dataframe
        filtered 1-1 orthologue pairs
    species: str
        species ID

    Returns:
    --------
    sp_map: dict
        species physical chromosome map with genes in order, identified primarily by their OG ID
    '''

    df_sp = df[['Orthogroup', species]]
    df_sp_core = df_sp[df_sp.apply(lambda x: all((x.str.len() > 0)), axis=1)].reset_index(drop=True)

    sp_tsv = pd.read_csv(f"{tsv_dir}/{species}.tsv", sep = "\t")
    d = {}
    pids = list(sp_tsv['ProteinID'])
    for k in pids:
        for j,i in enumerate(df_sp_core[species]):
            if k in i:
                d[k] = df_sp_core['Orthogroup'][j]
        d.setdefault(k, 'NA')

    ogs = list(d.values())
    sp_tsv['Orthogroups'] = ogs

    sp_map = {}
    for index, row in sp_tsv.iterrows():
        if row['Orthogroups'] != 'NA':
            sp_map.setdefault(row['Chr/Scaffold'], []).append((row['Orthogroups'], row['ProteinID'], row['Start'], row['End']))
    return sp_map

## Pairwise Comparison Map Creation

def create_comparison_map(species1_map, species2_map, species1_name, species2_name):
    '''
    Input the 2 physical chromosome maps for a pair of species and create
    comparison map (dict format).
    
    Parameters:
    -----------
    species1_map: dict
        species 1 physical chr map
    species2_map: dict
        species 2 physical chr map
    species1_name: str
        species 1 ID
    species2_name: str
        species 2 ID
    
    Returns:
    --------
    comparison_map: dict
        species 1 vs species 2 comparison map
    '''

    species1_og_positions = {}
    for chrom, genes in species1_map.items():
        for pos_idx, (og_id, protein_id, start, end) in enumerate(genes):
            species1_og_positions[og_id] = (species1_name, chrom, pos_idx)
    
    species2_og_positions = {}
    for chrom, genes in species2_map.items():
        for pos_idx, (og_id, protein_id, start, end) in enumerate(genes):
            species2_og_positions[og_id] = (species2_name, chrom, pos_idx)

    comparison_map = {}

    shared_ogs = set(species1_og_positions.keys()) & set(species2_og_positions.keys())
    for og_id in shared_ogs:
        comparison_map[og_id] = [
            species1_og_positions[og_id],
            species2_og_positions[og_id]
        ]
    
    return comparison_map

## Custom sort function: numbers first, then X, Y, Z, MT, etc.

def chrom_sort_key(chrom):
    '''
    Sort chromosomes in a biologically meaningful order.
    
    Handles multiple chromosome naming conventions and sorts them in the following order:
    1. Numeric chromosomes (1, 2, 3, ..., 22)
    2. Roman numeral chromosomes (I, II, III, ..., XX)
    3. Letter-number combinations (BMI1, BMI2, ..., BMI10, ...)
    4. Sex chromosomes (X, Y, Z, W)
    5. Mitochondrial chromosomes (MT, M, Mt, mt)
    6. Other chromosomes/scaffolds (alphabetically)
    
    Parameters:
    -----------
    chrom : str or int
        Chromosome name or identifier
    
    Returns:
    --------
    tuple
        A tuple (priority, sort_value) where:
        - priority determines the sorting group (0=autosomes, 1=sex chr, 2=MT, 3=other)
        - sort_value is the numeric or string value for sorting within that group
    '''
    # Handle None
    if chrom is None:
        return (999, 'None')
    
    # Handle NaN (pandas specific)
    if isinstance(chrom, float) and np.isnan(chrom):
        return (999, 'NaN')
    
    # Handle empty strings
    if chrom == '' or str(chrom).strip() == '':
        return (999, 'Empty')
    
    chrom_str = str(chrom)
    
    # Map Roman numerals to numbers
    roman_map = {'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5, 
                 'VI': 6, 'VII': 7, 'VIII': 8, 'IX': 9, 'X': 10,
                 'XI': 11, 'XII': 12, 'XIII': 13, 'XIV': 14, 'XV': 15,
                 'XVI': 16, 'XVII': 17, 'XVIII': 18, 'XIX': 19, 'XX': 20}
    
    # Try numeric chromosomes first
    if chrom_str.isdigit():
        return (0, int(chrom_str))
    
    # Check if it's a Roman numeral
    elif chrom_str in roman_map and (len(chrom_str) > 1 or chrom_str in ('I', 'V')):
        return (0, roman_map[chrom_str])
    
    elif '_sca' in chrom_str.lower():
        # Extract the scaffold number
        match = re.search(r'sca(\d+)', chrom_str.lower())
        if match:
            scaffold_num = int(match.group(1))
            return (4, scaffold_num)  # Category 4 for scaffolds
        else:
            return (4, chrom_str)  # Scaffolds without numbers
    
    elif chrom_str[-1].isdigit():
        # Match: optional letters/underscores, then digits
        match = re.match(r'^([A-Za-z_]+?)(\d+)$', chrom_str)
        if match:
            prefix = match.group(1)
            number = int(match.group(2))
            return (0, number)  # Regular chromosomes
        else:
            return (3, chrom_str)  # Fallback
    
    # Sex chromosomes
    elif chrom_str == 'X':
        return (1, 0)
    elif chrom_str == 'Y':
        return (1, 1)
    elif chrom_str == 'Z':
        return (1, 2)
    elif chrom_str == 'W':
        return (1, 3)
    
    # Mitochondrial
    elif chrom_str in ('MT', 'M', 'Mt', 'mt'):
        return (2, 0)
    
    # Any other chromosomes
    else:
        return (3, chrom_str)

## Fisher's Test

def fishers(comparison_map, alpha=0.01, min_matches=False, gene_filtering=False):
    '''
    Perform Fisher's exact test on all chromosome pairs to identify statistically significant synteny.
    
    For each pair of chromosomes between two species, constructs a 2x2 contingency table and
    tests whether the observed number of shared orthologs is significantly greater than expected
    by chance. Applies Bonferroni correction for multiple testing.
    
    Parameters:
    -----------
    comparison_map : dict
        Dictionary mapping orthogroup IDs to ortholog pairs.
        Format: {OG_ID: [(species1, chromosome, position), (species2, chromosome, position)]}
    alpha : float, default=0.01
        Significance threshold for Bonferroni-corrected p-values
    min_matches : bool, default=False
        If True, ensures every chromosome has at least one significant match by adding
        the best match for chromosomes without any significant pairs
    gene_filtering : bool, default=False
        If True, returns a filtered comparison_map containing only orthologs from
        statistically significant chromosome pairs
    
    Returns:
    --------
    If gene_filtering=False (default):
        counts_df : pd.DataFrame
            Contingency table showing number of shared orthologs between each chromosome pair
        results_df : pd.DataFrame
            Complete results with columns: sp1_chr, sp2_chr, a, b, c, d, odds_ratio,
            p_value, corrected_p, significant
        significant_pairs : list of tuples
            List of significant chromosome pairs: (chr1, chr2, count, p_value, corrected_p)
    
    If gene_filtering=True:
        filtered_comparison_map : dict
            Comparison map containing only orthologs from significant chromosome pairs
    '''

    counts = defaultdict(lambda: defaultdict(int))
    for i in comparison_map.values():
        outer_key = i[0][1]
        inner_key = i[1][1]
        counts[outer_key][inner_key] += 1

    counts_df = (
        pd.DataFrame.from_dict(counts, orient="index")
          .fillna(0)
          .astype(int)
    )
    
    counts_df = counts_df.loc[
        sorted(counts_df.index, key=chrom_sort_key),
        sorted(counts_df.columns, key=chrom_sort_key)
    ]

    # Total number of orthologs
    total_orthologs = len(comparison_map)
    
    # Get all chromosome pairs
    sp1_chrs = counts_df.index.tolist()
    sp2_chrs = counts_df.columns.tolist()
    
    # Calculate total per chromosome
    sp1_totals = counts_df.sum(axis=1)  # Total orthologs per human chromosome
    sp2_totals = counts_df.sum(axis=0)   # Total orthologs per mouse chromosome
    
    results = []
    
    # Run Fisher's exact test for each chromosome pair
    for sp1_chr in sp1_chrs:
        for sp2_chr in sp2_chrs:
            # a = orthologs on BOTH chromosomes
            a = counts_df.loc[sp1_chr, sp2_chr]
            
            # b = orthologs on sp1_chr but NOT sp2_chr
            b = sp1_totals[sp1_chr] - a
            
            # c = orthologs on sp2_chr but NOT sp1_chr
            c = sp2_totals[sp2_chr] - a
            
            # d = orthologs on neither
            d = total_orthologs - a - b - c
            
            # Construct 2x2 contingency table
            table = [[a, b],
                     [c, d]]
            
            # Run Fisher's exact test (one-sided, testing for enrichment)
            odds_ratio, p_value = fisher_exact(table, alternative='greater')
            
            results.append({
                'sp1_chr': sp1_chr,
                'sp2_chr': sp2_chr,
                'a': a,
                'b': b,
                'c': c,
                'd': d,
                'odds_ratio': odds_ratio,
                'p_value': p_value
            })
    
    # Create results dataframe
    results_df = pd.DataFrame(results)
    
    # Apply Bonferroni correction
    n_tests = len(results_df)
    results_df['corrected_p'] = results_df['p_value'] * n_tests
    results_df['corrected_p'] = results_df['corrected_p'].clip(upper=1.0)  # Cap at 1.0
    
    # Mark significant pairs
    results_df['significant'] = results_df['corrected_p'] < alpha
    
    # Get significant pairs
    significant_pairs = [
        (row['sp1_chr'], row['sp2_chr'], row['a'], row['p_value'], row['corrected_p']) 
        for _, row in results_df.iterrows() 
        if row['significant']
    ]
    
    # Sort by corrected p-value
    results_df = results_df.sort_values('corrected_p')

    #If min_matches is True, Run Fisher's exact test with progressive alpha
    #relaxation to ensure every chromosome has at least one significant match.
    if min_matches:
    
        # Check which chromosomes have no matches
        sp1_chroms_with_matches = set(pair[0] for pair in significant_pairs)
        sp2_chroms_with_matches = set(pair[1] for pair in significant_pairs)
    
        all_sp1_chroms = set(counts_df.index)
        all_sp2_chroms = set(counts_df.columns)
    
        sp1_chroms_without_matches = all_sp1_chroms - sp1_chroms_with_matches
        sp2_chroms_without_matches = all_sp2_chroms - sp2_chroms_with_matches
    
        # If all chromosomes have matches, we're done
        if not sp1_chroms_without_matches and not sp2_chroms_without_matches:
            return counts_df, results_df, significant_pairs
    
        print(f"Initial alpha={alpha}: {len(significant_pairs)} pairs")
        print(f"  Sp1 chroms without matches: {sp1_chroms_without_matches}")
        print(f"  Sp2 chroms without matches: {sp2_chroms_without_matches}")
    
        # For chromosomes without matches, find their best match
        # even if it doesn't meet the alpha threshold
        additional_pairs = []
    
        for chrom in sp1_chroms_without_matches:
            # Find best match for this sp1 chromosome
            best_match = results_df[results_df['sp1_chr'] == chrom].iloc[0]
            additional_pairs.append((best_match['sp1_chr'], best_match['sp2_chr'], best_match['a'], best_match['p_value'], best_match['corrected_p']))
            print(f"  Adding best match for {chrom}: {best_match['sp2_chr']} (corrected_p={best_match['corrected_p']:.4f})")

        for chrom in sp2_chroms_without_matches:
            # Find best match for this sp2 chromosome
            best_match = results_df[results_df['sp2_chr'] == chrom].iloc[0]
            # Only add if not already in significant_pairs or additional_pairs
            pair = (best_match['sp1_chr'], best_match['sp2_chr'], best_match['a'], best_match['p_value'], best_match['corrected_p'])
            if pair not in significant_pairs and pair not in additional_pairs:
                additional_pairs.append(pair)
                print(f"  Adding best match for {chrom}: {best_match['sp1_chr']} (corrected_p={best_match['corrected_p']:.4f})")
    
        # Combine original significant pairs with additional pairs
        significant_pairs = significant_pairs + additional_pairs

        print(f"Final: {len(significant_pairs)} pairs total")

    # If gene_filtering is True, filter comparison_map
    if gene_filtering:
        # Create set of significant chromosome pairs for fast lookup
        significant_set = {(sp1_chr, sp2_chr) for sp1_chr, sp2_chr, count, p_val, corrected_p in significant_pairs}
        
        # Filter comparison_map to only keep orthologs in significant chromosome pairs
        filtered_comparison_map = {}
        for og_id, orthologs in comparison_map.items():
            sp1_chrom = orthologs[0][1]  # chromosome from species 1
            sp2_chrom = orthologs[1][1]  # chromosome from species 2
            
            # Only keep if this chromosome pair is significant
            if (sp1_chrom, sp2_chrom) in significant_set:
                filtered_comparison_map[og_id] = orthologs
        
        print(f"Filtered from {len(comparison_map)} to {len(filtered_comparison_map)} orthologs")
        print(f"Kept genes from {len(significant_pairs)} significant chromosome pairs")
        
        return filtered_comparison_map
    
    return counts_df, results_df, significant_pairs

def fishers_multi(species_maps, species_names, alpha=0.01, min_matches=False):
    '''
    Perform pairwise Fisher's exact tests between ALL species pairs
    (not just adjacent) and create a multi-species comparison map.
    
    Parameters:
    -----------
    species_maps : list of dict
        List of physical chromosome maps for each species
    species_names : list of str
        List of species IDs in same order as species_maps
    alpha : float
        Significance threshold
    min_matches : bool
        Ensure every chromosome has at least one match
    
    Returns:
    --------
    multi_filtered_map : dict
        Format: {og_id: {
            'positions': [('sp', 'chr', pos) or None, ...],
            'significant_segments': [(i, j), ...]  # ALL significant pairs now
        }}
    '''
    
    n_species = len(species_names)
    
    print(f"\n{'='*60}")
    print(f"MULTI-SPECIES FISHER'S EXACT TEST (ALL PAIRS)")
    print(f"{'='*60}")
    
    # Test ALL pairs, not just adjacent
    all_pairs = list(itertools.combinations(range(n_species), 2))
    
    pairwise_sig_genes = {}
    
    for i, j in all_pairs:
        sp1_name = species_names[i]
        sp2_name = species_names[j]
        sp1_map = species_maps[i]
        sp2_map = species_maps[j]
        
        print(f"\n--- Testing {sp1_name} vs {sp2_name} ---")
        
        pairwise_comparison_map = create_comparison_map(
            sp1_map, sp2_map, sp1_name, sp2_name
        )
        
        print(f"Total orthogroups in pair: {len(pairwise_comparison_map)}")
        
        filtered_pairwise_map = fishers(
            pairwise_comparison_map,
            alpha=alpha,
            min_matches=min_matches,
            gene_filtering=True
        )
        
        print(f"Significant orthogroups after Fisher's: {len(filtered_pairwise_map)}")
        
        sig_genes = {}
        for og_id, positions in filtered_pairwise_map.items():
            chr1 = positions[0][1]
            chr2 = positions[1][1]
            sig_genes[og_id] = (chr1, chr2)
        
        pairwise_sig_genes[(i, j)] = sig_genes
    
    # Build position lookup for each species
    species_og_positions = []
    for sp_idx, (species_name, species_map) in enumerate(zip(species_names, species_maps)):
        og_positions = {}
        for chrom, genes in species_map.items():
            for pos_idx, (og_id, protein_id, start, end) in enumerate(genes):
                og_positions[og_id] = (species_name, chrom, pos_idx)
        species_og_positions.append(og_positions)
    
    # Collect all significant orthogroups
    all_significant_ogs = set()
    for sig_genes in pairwise_sig_genes.values():
        all_significant_ogs.update(sig_genes.keys())
    
    print(f"\n{'='*60}")
    print(f"BUILDING MULTI-SPECIES MAP")
    print(f"{'='*60}")
    print(f"Total unique orthogroups: {len(all_significant_ogs)}")
    
    # Build multi-species filtered map
    multi_filtered_map = {}
    
    for og_id in all_significant_ogs:
        positions = []
        for sp_idx in range(n_species):
            if og_id in species_og_positions[sp_idx]:
                positions.append(species_og_positions[sp_idx][og_id])
            else:
                positions.append(None)
        
        # Store ALL significant pairs (not just adjacent segments)
        # Format changes: significant_segments now stores (i, j) tuples
        significant_pairs = []
        
        for i, j in all_pairs:
            if og_id in pairwise_sig_genes[(i, j)]:
                if positions[i] is not None and positions[j] is not None:
                    significant_pairs.append((i, j))
        
        if significant_pairs:
            multi_filtered_map[og_id] = {
                'positions': positions,
                'significant_segments': significant_pairs  # Now stores (i,j) tuples
            }
    
    # Print summary
    print(f"\n{'='*60}")
    print(f"FINAL RESULTS")
    print(f"{'='*60}")
    print(f"Total genes in multi-species map: {len(multi_filtered_map)}")
    
    for i, j in all_pairs:
        count = sum(
            1 for data in multi_filtered_map.values()
            if (i, j) in data['significant_segments']
        )
        print(f"Genes significant in {species_names[i]}-{species_names[j]}: {count}")
    
    return multi_filtered_map

## Create Color Mapping for up to 30 Distinct Chromosomes for Dot Plots and Ribbon Plots 

custom_colors = [
    # ===== BATCH 1: Colors 1-25 (maximally distinct) =====
    "#FF0000",  # 1.  Pure Red
    "#0000FF",  # 2.  Pure Blue
    "#00CC00",  # 3.  Pure Green
    "#FF00FF",  # 4.  Magenta
    "#FF8C00",  # 5.  Dark Orange
    "#00FFFF",  # 6.  Cyan
    "#8B008B",  # 7.  Dark Magenta
    "#FFFF00",  # 8.  Yellow
    "#006400",  # 9.  Dark Green
    "#FF1493",  # 10. Deep Pink
    "#00008B",  # 11. Dark Blue
    "#FF6347",  # 12. Tomato
    "#00FA9A",  # 13. Medium Spring Green
    "#9400D3",  # 14. Dark Violet
    "#FFA500",  # 15. Orange
    "#1E90FF",  # 16. Dodger Blue
    "#8B0000",  # 17. Dark Red
    "#00CED1",  # 18. Dark Turquoise
    "#556B2F",  # 19. Dark Olive Green
    "#FF69B4",  # 20. Hot Pink
    "#4B0082",  # 21. Indigo
    "#ADFF2F",  # 22. Green Yellow
    "#A0522D",  # 23. Sienna
    "#40E0D0",  # 24. Turquoise
    "#DC143C",  # 25. Crimson

    # ===== BATCH 2: Colors 26-50 (distinct from each other, ok if close to batch 1) =====
    "#FA8072",  # 26. Salmon          (close to Red/Tomato)
    "#4169E1",  # 27. Royal Blue      (close to Blue)
    "#90EE90",  # 28. Light Green     (close to Green)
    "#DDA0DD",  # 29. Plum            (close to Magenta)
    "#FFD700",  # 30. Gold            (close to Orange)
    "#87CEEB",  # 31. Sky Blue        (close to Cyan)
    "#BA55D3",  # 32. Medium Orchid   (close to Dark Magenta)
    "#F0E68C",  # 33. Khaki           (close to Yellow)
    "#228B22",  # 34. Forest Green    (close to Dark Green)
    "#FFB6C1",  # 35. Light Pink      (close to Deep Pink)
    "#6495ED",  # 36. Cornflower Blue (close to Dark Blue)
    "#E9967A",  # 37. Dark Salmon     (close to Tomato)
    "#98FB98",  # 38. Pale Green      (close to Spring Green)
    "#9370DB",  # 39. Medium Purple   (close to Dark Violet)
    "#FFDAB9",  # 40. Peach Puff      (close to Orange)
    "#B0C4DE",  # 41. Light Steel Blue(close to Dodger Blue)
    "#CD5C5C",  # 42. Indian Red      (close to Dark Red)
    "#AFEEEE",  # 43. Pale Turquoise  (close to Dark Turquoise)
    "#6B8E23",  # 44. Olive Drab      (close to Olive Green)
    "#FFB7C5",  # 45. Cherry Blossom  (close to Hot Pink)
    "#7B68EE",  # 46. Medium Slate Blue (close to Indigo)
    "#CCFF66",  # 47. Lime Fizz       (close to Green Yellow)
    "#C68642",  # 48. Peru            (close to Sienna)
    "#7FFFD4",  # 49. Aquamarine      (close to Turquoise)
    "#F08080",  # 50. Light Coral     (close to Crimson)

    # ===== BATCH 3: Colors 51-70 (distinct from each other, ok if close to batch 1 or 2) =====
    "#B22222",  # 51. Firebrick       (close to Red family)
    "#191970",  # 52. Midnight Blue   (close to Blue family)
    "#32CD32",  # 53. Lime Green      (close to Green family)
    "#FF00AA",  # 54. Rose            (close to Magenta family)
    "#D2691E",  # 55. Chocolate       (close to Orange family)
    "#00B2EE",  # 56. Deep Sky Blue   (close to Cyan family)
    "#6A0DAD",  # 57. Purple          (close to Violet family)
    "#DAA520",  # 58. Goldenrod       (close to Yellow family)
    "#2E8B57",  # 59. Sea Green       (close to Green family)
    "#FF82AB",  # 60. Pink            (close to Pink family)
    "#003153",  # 61. Prussian Blue   (close to Dark Blue family)
    "#FF4500",  # 62. Orange Red      (close to Red/Orange family)
    "#00BFFF",  # 63. Deep Sky Blue 2 (close to Turquoise family)
    "#7CFC00",  # 64. Lawn Green      (close to Green Yellow family)
    "#C71585",  # 65. Medium Violet Red (close to Pink family)
    "#8FBC8F",  # 66. Dark Sea Green  (close to Green family)
    "#483D8B",  # 67. Dark Slate Blue (close to Indigo family)
    "#FF8247",  # 68. Sienna 1        (close to Orange/Brown family)
    "#20B2AA",  # 69. Light Sea Green (close to Teal family)
    "#EE1289",  # 70. DeepPink 2      (close to Pink/Magenta family)
]

## Dot Plot Creation

def plot_synteny_dotplot(comparison_map, sp1_map, sp2_map, significant_pairs, species1, species2, 
                          figsize=(15, 15), dot_size=1, dot_alpha=0.5, default_color='#909090'):
    """
    Create a synteny dot plot showing ortholog positions between two species.
    Chromosomes are scaled by their actual genomic length, but axes show cumulative gene counts.
    
    Parameters:
    -----------
    comparison_map : dict
        Dictionary with OG IDs as keys and list of tuples as values
        Each tuple contains (species_name, chromosome, position)
    sp1_map : dict
        Dictionary with chromosome as key and list of gene tuples as values for species 1
        Each tuple contains (OG_ID, protein_ID, start, end)
    sp2_map : dict
        Dictionary with chromosome as key and list of gene tuples as values for species 2
        Each tuple contains (OG_ID, protein_ID, start, end)
    significant_pairs : list
        List of tuples (sp1_chr, sp2_chr, count) from Fisher's exact test
    species1 : str
        species 1 name
    species2 : str
        species 2 name
    figsize : tuple
        Figure size (default: (15, 15))
    dot_size : float
        Size of dots (default: 1)
    dot_alpha : float
        Transparency of dots (default: 0.5)
    default_color : str
        Color for non-significant dots (default: 'lightgray')
    """
    
    # Extract data for each species
    data = []
    for og_id, orthologs in comparison_map.items():
        sp1_data = None
        sp2_data = None
        
        for species, chrom, pos in orthologs:
            if species == species1:
                sp1_data = (chrom, pos)
            elif species == species2:
                sp2_data = (chrom, pos)
        
        if sp1_data and sp2_data and sp1_data[1] > 0 and sp2_data[1] > 0:
            data.append((sp1_data[0], sp1_data[1], sp2_data[0], sp2_data[1]))
    
    # Get unique chromosomes that appear in the comparison data and sort them
    sp1_chroms = sorted(set(d[0] for d in data), key=chrom_sort_key)
    sp2_chroms = sorted(set(d[2] for d in data), key=chrom_sort_key)

    # Get chromosome lengths from the maps (maximum end coordinate)
    # ONLY for chromosomes that have orthologs
    sp1_chrom_lengths = {}
    for chrom in sp1_chroms:
        if chrom in sp1_map:
            max_end = max(gene[3] for gene in sp1_map[chrom])
            sp1_chrom_lengths[chrom] = max_end

    sp2_chrom_lengths = {}
    for chrom in sp2_chroms:
        if chrom in sp2_map:
            max_end = max(gene[3] for gene in sp2_map[chrom])
            sp2_chrom_lengths[chrom] = max_end

    all_colors = np.array([plt.matplotlib.colors.to_rgba(c) for c in custom_colors])
    
    chrom_to_color = {}
    for idx, chrom in enumerate(sp1_chroms):
        chrom_to_color[chrom] = all_colors[idx % len(custom_colors)]
    
    # Create set of significant pairs for fast lookup
    significant_set = set()
    if significant_pairs:
        significant_set = {(sp1_chr, sp2_chr) for sp1_chr, sp2_chr, count, p_val, corrected_p in significant_pairs}
    
    # Count genes per chromosome in the comparison data
    sp1_gene_counts = {}
    sp2_gene_counts = {}
    
    for chrom in sp1_chroms:
        sp1_gene_counts[chrom] = max(d[1] for d in data if d[0] == chrom)
    
    for chrom in sp2_chroms:
        sp2_gene_counts[chrom] = max(d[3] for d in data if d[2] == chrom)
    
    # Calculate cumulative gene positions for plotting
    sp1_cumulative_genes = {}
    sp2_cumulative_genes = {}
    
    cumsum = 0
    for chrom in sp1_chroms:
        sp1_cumulative_genes[chrom] = cumsum
        cumsum += sp1_gene_counts[chrom]
    sp1_total_genes = cumsum
    
    cumsum = 0
    for chrom in sp2_chroms:
        sp2_cumulative_genes[chrom] = cumsum
        cumsum += sp2_gene_counts[chrom]
    sp2_total_genes = cumsum
    
    # Calculate cumulative lengths (for scaling chromosome boundaries)
    sp1_cumulative_length = {}
    sp2_cumulative_length = {}
    
    cumsum = 0
    for chrom in sp1_chroms:
        sp1_cumulative_length[chrom] = cumsum
        cumsum += sp1_chrom_lengths[chrom]
    sp1_total_length = cumsum
    
    cumsum = 0
    for chrom in sp2_chroms:
        sp2_cumulative_length[chrom] = cumsum
        cumsum += sp2_chrom_lengths[chrom]
    sp2_total_length = cumsum
    
    # Create the plot
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot all orthologs
    x_coords = []
    y_coords = []
    
    # Separate dots into significant and non-significant
    significant_dots = {chrom: {'x': [], 'y': []} for chrom in sp1_chroms}
    non_significant_x = []
    non_significant_y = []
    
    for sp1_chrom, sp1_pos, sp2_chrom, sp2_pos in data:
        # Map position within chromosome to proportional position in length space
        fraction_x = sp1_pos / sp1_gene_counts[sp1_chrom]
        fraction_y = sp2_pos / sp2_gene_counts[sp2_chrom]
        
        x = sp1_cumulative_length[sp1_chrom] + fraction_x * sp1_chrom_lengths[sp1_chrom]
        y = sp2_cumulative_length[sp2_chrom] + fraction_y * sp2_chrom_lengths[sp2_chrom]
        
        # Check if this pair is significant
        if (sp1_chrom, sp2_chrom) in significant_set:
            significant_dots[sp1_chrom]['x'].append(x)
            significant_dots[sp1_chrom]['y'].append(y)
        else:
            non_significant_x.append(x)
            non_significant_y.append(y)
    
    # Plot non-significant dots first (so they're in the background)
    if non_significant_x:
        ax.scatter(non_significant_x, non_significant_y, s=dot_size, alpha=dot_alpha/2, 
                   c=default_color, rasterized=True, label='Non-significant')
    
    # Plot significant dots by chromosome (each with its own color)
    for chrom in sp1_chroms:
        if significant_dots[chrom]['x']:
            ax.scatter(significant_dots[chrom]['x'], significant_dots[chrom]['y'], 
                      s=dot_size, alpha=dot_alpha, c=[chrom_to_color[chrom]], 
                      rasterized=True, label=f'{species1} Chr{chrom}')
    
    # Set up axes limits
    ax.set_xlim(0, sp1_total_length)
    ax.set_ylim(0, sp2_total_length)
    
    # Add chromosome boundaries as grid lines
    for chrom in sp1_chroms:
        ax.axvline(sp1_cumulative_length[chrom], color='black', linewidth=0.5, alpha=0.5)
    for chrom in sp2_chroms:
        ax.axhline(sp2_cumulative_length[chrom], color='black', linewidth=0.5, alpha=0.5)
    
    # Set up tick marks showing cumulative gene counts
    # Map gene count positions to length space for x-axis
    gene_tick_interval = 1000
    x_gene_ticks = list(range(0, sp1_total_genes + 1, gene_tick_interval))
    x_length_ticks = []
    x_labels = []
    
    for gene_count in x_gene_ticks:
        # Find which chromosome this gene count falls into
        current_chrom = None
        for chrom in sp1_chroms:
            chrom_start = sp1_cumulative_genes[chrom]
            chrom_end = chrom_start + sp1_gene_counts[chrom]
            if chrom_start <= gene_count <= chrom_end:
                current_chrom = chrom
                relative_pos = gene_count - chrom_start
                fraction = relative_pos / sp1_gene_counts[chrom]
                length_pos = sp1_cumulative_length[chrom] + fraction * sp1_chrom_lengths[chrom]
                x_length_ticks.append(length_pos)
                x_labels.append(str(gene_count))
                break
    
    # Map gene count positions to length space for y-axis
    y_gene_ticks = list(range(0, sp2_total_genes + 1, gene_tick_interval))
    y_length_ticks = []
    y_labels = []
    
    for gene_count in y_gene_ticks:
        for chrom in sp2_chroms:
            chrom_start = sp2_cumulative_genes[chrom]
            chrom_end = chrom_start + sp2_gene_counts[chrom]
            if chrom_start <= gene_count <= chrom_end:
                relative_pos = gene_count - chrom_start
                fraction = relative_pos / sp2_gene_counts[chrom]
                length_pos = sp2_cumulative_length[chrom] + fraction * sp2_chrom_lengths[chrom]
                y_length_ticks.append(length_pos)
                y_labels.append(str(gene_count))
                break
    
    ax.set_xticks(x_length_ticks)
    ax.set_yticks(y_length_ticks)
    ax.set_xticklabels(x_labels, fontsize=8)
    ax.set_yticklabels(y_labels, fontsize=8)
    
    # Add chromosome labels on secondary axes (top and right)
    ax2 = ax.twiny()
    ax3 = ax.twinx()
    
    # Calculate midpoints for chromosome labels based on length
    sp1_midpoints = [sp1_cumulative_length[chrom] + sp1_chrom_lengths[chrom]/2 for chrom in sp1_chroms]
    sp2_midpoints = [sp2_cumulative_length[chrom] + sp2_chrom_lengths[chrom]/2 for chrom in sp2_chroms]
    
    ax2.set_xlim(0, sp1_total_length)
    ax2.set_xticks(sp1_midpoints)
    ax2.set_xticklabels(sp1_chroms, fontsize=10)
    ax2.set_xlabel(f'{species1}', fontsize=12, fontweight='bold')
    
    ax3.set_ylim(0, sp2_total_length)
    ax3.set_yticks(sp2_midpoints)
    ax3.set_yticklabels(sp2_chroms, fontsize=10)
    ax3.set_ylabel(f'{species2}', fontsize=12, fontweight='bold', rotation=270, labelpad=20)
    
    # Label bottom and left axes
    ax.set_xlabel('Cumulative gene count', fontsize=10)
    ax.set_ylabel('Cumulative gene count', fontsize=10)
    
    # Add light grid (optional)
    #ax.grid(True, alpha=0.2, linewidth=0.5, linestyle='dotted')
    ax.grid(False)
    
    plt.tight_layout()

    # Print Fisher's test results for significant pairs
    if significant_pairs:
        print("\n" + "="*100)
        print(f"SIGNIFICANT SYNTENIC CHROMOSOME PAIRS ({species1} vs {species2})")
        print("="*100)
        print(f"{'Sp1 Chr':<12} {'Sp2 Chr':<12} {'Shared':<12} {'Raw p-value':<20} {'Corrected p':<20}")
        print("-"*100)
        
        for sp1_chr, sp2_chr, count, p_val, corrected_p in significant_pairs:
            print(f"{str(sp1_chr):<12} {str(sp2_chr):<12} {count:<12} {p_val:<20.2e} {corrected_p:<20.2e}")
        
        print("="*100 + "\n")
    
    return fig, ax

## ALG Discovery

def compute_synteny_similarity(species_maps, species_names, fishers_multi_results):
    """
    Compute pairwise synteny similarity between all species.
    Uses direct Fisher's significant pairs (all pairs tested).
    """
    n = len(species_names)
    similarity = np.zeros((n, n))
    
    print("\n" + "="*70)
    print("COMPUTING SYNTENY SIMILARITY MATRIX")
    print("="*70)
    
    for i in range(n):
        for j in range(i+1, n):
            sp1, sp2 = species_names[i], species_names[j]
            
            shared_orthologs = 0
            total_orthologs = 0
            
            for og_id, data in fishers_multi_results.items():
                positions = data['positions']
                significant_segments = data['significant_segments']
                
                if positions[i] is not None and positions[j] is not None:
                    total_orthologs += 1
                    
                    # significant_segments is now a list of (i,j) tuples
                    # Check direct significance between this pair
                    if (i, j) in significant_segments or (j, i) in significant_segments:
                        shared_orthologs += 1
            
            sim = shared_orthologs / total_orthologs if total_orthologs > 0 else 0.0
            similarity[i, j] = sim
            similarity[j, i] = sim
            
            print(f"{sp1} <-> {sp2}: {shared_orthologs}/{total_orthologs} = {sim:.3f}")
    
    np.fill_diagonal(similarity, 1.0)
    return similarity

def cluster_species_by_synteny(similarity_matrix, species_names, 
                                threshold=0.15, method='average'):
    """
    Cluster species into synteny groups using hierarchical clustering
    on a pairwise synteny similarity matrix.
    
    Converts the similarity matrix to a distance matrix (distance = 1 - similarity),
    then applies scipy hierarchical clustering to group species that share more
    synteny signal than the given threshold.
    
    Parameters:
    -----------
    similarity_matrix : np.ndarray
        N x N symmetric matrix of pairwise synteny similarities (output of
        compute_synteny_similarity). Values range from 0 (no shared synteny)
        to 1 (identical synteny).
    species_names : list of str
        List of species IDs in the same order as the rows/columns of
        similarity_matrix.
    threshold : float
        Similarity threshold for cluster separation (default: 0.15).
        Species pairs with similarity below this threshold will be placed
        in different clusters. Increase to merge more species into fewer
        clusters; decrease to split into more clusters.
    method : str
        Linkage method for hierarchical clustering (default: 'average').
        Passed directly to scipy.cluster.hierarchy.linkage.
        Options: 'single', 'complete', 'average', 'ward', etc.
    
    Returns:
    --------
    clusters : dict
        Dictionary mapping cluster ID (int) to list of species names.
        Example: {1: ['HCA', 'BMI'], 2: ['CLA', 'EMU', 'BFL', 'RES']}
    cluster_labels : np.ndarray
        Array of cluster label integers, one per species, in the same
        order as species_names.
        Example: array([2, 2, 2, 2, 1, 1]) for 6 species
    """
    
    print("\n" + "="*70)
    print("HIERARCHICAL CLUSTERING OF SPECIES")
    print("="*70)
    print(f"Threshold: {threshold} (distance = {1-threshold:.3f})")
    
    distance = 1 - similarity_matrix
    condensed_dist = squareform(distance, checks=False)
    Z = linkage(condensed_dist, method=method)
    cluster_labels = fcluster(Z, 1 - threshold, criterion='distance')
    
    clusters = defaultdict(list)
    for idx, label in enumerate(cluster_labels):
        clusters[label].append(species_names[idx])
    
    print(f"\nFound {len(clusters)} synteny cluster(s):")
    for cluster_id, members in sorted(clusters.items()):
        print(f"  Cluster {cluster_id}: {', '.join(members)}")
    
    return dict(clusters), cluster_labels


def precompute_significant_chrom_pairs(cluster_species, species_names,
                                        fishers_multi_results, min_orthologs=10):
    """
    Extract significant chromosome pairs for all species combinations
    within a cluster, directly from fishers_multi results.

    Parameters:
    -----------
    cluster_species : list of str
        Species IDs belonging to this cluster
    species_names : list of str
        Full species list (all species, not just cluster)
    fishers_multi_results : dict
        Output of fishers_multi
    min_orthologs : int
        Minimum orthologs required for a chromosome pair to be retained (default: 10)

    Returns:
    --------
    sig_chrom_pairs : dict
        {(sp1, sp2): {(chr1, chr2): count}}
        For each species pair, a dictionary of significant chromosome pairs
        and the number of orthologs supporting them.
    """
    
    cluster_indices = {sp: species_names.index(sp) for sp in cluster_species}
    pair_counts = defaultdict(lambda: defaultdict(int))
    
    for og_id, data in fishers_multi_results.items():
        positions = data['positions']
        significant_segments = data['significant_segments']  # List of (i, j) tuples
        
        for sp1 in cluster_species:
            for sp2 in cluster_species:
                if sp1 >= sp2:
                    continue
                
                idx1 = cluster_indices[sp1]
                idx2 = cluster_indices[sp2]
                
                if positions[idx1] is None or positions[idx2] is None:
                    continue
                
                # Check if this pair is directly significant
                is_significant = (
                    (idx1, idx2) in significant_segments or
                    (idx2, idx1) in significant_segments
                )
                
                if is_significant:
                    chr1 = positions[idx1][1]
                    chr2 = positions[idx2][1]
                    pair_counts[(sp1, sp2)][(chr1, chr2)] += 1
    
    # Filter by min_orthologs
    sig_chrom_pairs = {}
    for sp_pair, chrom_pair_counts in pair_counts.items():
        sig_chrom_pairs[sp_pair] = {
            chrom_pair: count
            for chrom_pair, count in chrom_pair_counts.items()
            if count >= min_orthologs
        }
    
    for sp1 in cluster_species:
        for sp2 in cluster_species:
            if sp1 >= sp2:
                continue
            sp_pair = (sp1, sp2)
            if sp_pair in sig_chrom_pairs:
                print(f"  {sp1} ↔ {sp2}: {len(sig_chrom_pairs[sp_pair])} significant chromosome pairs")
    
    return sig_chrom_pairs


def verify_alg_chain(chain, cluster_species, sig_chrom_pairs):
    """
    Verify that ALL pairwise chromosome combinations in a candidate ALG
    chain are Fisher's-significant and above the minimum ortholog threshold.

    Parameters:
    -----------
    chain : dict
        {species: chromosome} mapping for this candidate chain
    cluster_species : list of str
        Ordered list of species in this cluster
    sig_chrom_pairs : dict
        Output of precompute_significant_chrom_pairs

    Returns:
    --------
    bool
        True if all pairwise combinations are significant, False otherwise
    """
    
    species_in_chain = [sp for sp in cluster_species if sp in chain]
    
    if len(species_in_chain) < 2:
        return False
    
    for i in range(len(species_in_chain)):
        for j in range(i + 1, len(species_in_chain)):
            sp1 = species_in_chain[i]
            sp2 = species_in_chain[j]
            chr1 = chain[sp1]
            chr2 = chain[sp2]
            
            # Look up in both orderings
            if (sp1, sp2) in sig_chrom_pairs:
                sp_pair = (sp1, sp2)
                chrom_pair = (chr1, chr2)
            elif (sp2, sp1) in sig_chrom_pairs:
                sp_pair = (sp2, sp1)
                chrom_pair = (chr2, chr1)
            else:
                return False
            
            if chrom_pair not in sig_chrom_pairs[sp_pair]:
                return False
    
    return True


def identify_algs_for_cluster(cluster_species, species_maps, species_names, 
                                   fishers_multi_results, min_orthologs=10):
    """
    Identify Ancestral Linkage Groups for a single synteny cluster.

    Builds candidate chromosome chains by following significant adjacent
    pairwise matches, then verifies each chain by checking that ALL
    pairwise combinations of species are Fisher's-significant.

    Parameters:
    -----------
    cluster_species : list of str
        Ordered list of species in this cluster
    species_maps : list of dict
        Physical chromosome maps for all species (full list, not just cluster)
    species_names : list of str
        Full species list in the same order as species_maps
    fishers_multi_results : dict
        Output of fishers_multi
    min_orthologs : int
        Minimum orthologs for a chromosome pair to be considered (default: 10)

    Returns:
    --------
    alg_assignments : dict
        {species: {chromosome: [alg_id, ...]}}
        Each chromosome may belong to one or more ALGs (integer IDs, not yet
        cluster-prefixed). Cluster prefix is added by compute_algs_full_pipeline.
    """
    
    print(f"\n{'─'*70}")
    print(f"IDENTIFYING ALGs FOR CLUSTER: {', '.join(cluster_species)}")
    print(f"{'─'*70}")
    
    if len(cluster_species) < 2:
        print("⚠️  Only 1 species in cluster - cannot define ALGs")
        return {}
    
    cluster_indices = [species_names.index(sp) for sp in cluster_species]
    
    # Step 0: Precompute significant chromosome pairs
    print("\nPrecomputing significant chromosome pairs...")
    sig_chrom_pairs = precompute_significant_chrom_pairs(
        cluster_species, species_names, fishers_multi_results, min_orthologs
    )
    
    # Step 1: Build adjacent pairwise matches for chain construction
    pairwise_matches = {}
    
    for i in range(len(cluster_species) - 1):
        sp1 = cluster_species[i]
        sp2 = cluster_species[i + 1]
        
        if (sp1, sp2) in sig_chrom_pairs:
            sp_pair = (sp1, sp2)
            matches = [(chr1, chr2, count) 
                      for (chr1, chr2), count in sig_chrom_pairs[sp_pair].items()]
        elif (sp2, sp1) in sig_chrom_pairs:
            sp_pair = (sp2, sp1)
            matches = [(chr2, chr1, count) 
                      for (chr1, chr2), count in sig_chrom_pairs[sp_pair].items()]
        else:
            matches = []
        
        matches.sort(key=lambda x: x[2], reverse=True)
        pairwise_matches[(sp1, sp2)] = matches
        
        print(f"\n{sp1} ↔ {sp2}: {len(matches)} significant chromosome pairs")
        for chr1, chr2, count in matches[:5]:
            print(f"  {chr1} ↔ {chr2}: {count} orthologs")
        if len(matches) > 5:
            print(f"  ... and {len(matches) - 5} more")
    
    # Step 2: Build ALG chains recursively from adjacent matches
    def extend_chain(chain, remaining_species):
        if not remaining_species:
            return [chain]
        
        current_sp = remaining_species[0]
        next_remaining = remaining_species[1:]
        
        prev_sp = None
        for sp in reversed(cluster_species):
            if sp in chain:
                prev_sp = sp
                break
        
        if prev_sp is None:
            return [chain]
        
        if (prev_sp, current_sp) in pairwise_matches:
            matches = pairwise_matches[(prev_sp, current_sp)]
        elif (current_sp, prev_sp) in pairwise_matches:
            matches = [(chr2, chr1, cnt) 
                      for chr1, chr2, cnt in pairwise_matches[(current_sp, prev_sp)]]
        else:
            return [chain]
        
        prev_chr = chain[prev_sp]
        possible_extensions = [chr2 for chr1, chr2, cnt in matches if chr1 == prev_chr]
        
        if not possible_extensions:
            return [chain]
        
        all_chains = []
        for chr_current in possible_extensions:
            new_chain = chain.copy()
            new_chain[current_sp] = chr_current
            all_chains.extend(extend_chain(new_chain, next_remaining))
        
        return all_chains
    
    # Build chains starting from first species
    sp1 = cluster_species[0]
    sp1_chroms = sorted(species_maps[cluster_indices[0]].keys(), key=chrom_sort_key)
    
    all_candidate_chains = []
    for chr1 in sp1_chroms:
        initial_chain = {sp1: chr1}
        extended_chains = extend_chain(initial_chain, cluster_species[1:])
        all_candidate_chains.extend(extended_chains)
    
    print(f"\n{'─'*70}")
    print(f"VERIFYING {len(all_candidate_chains)} candidate chains...")
    print(f"{'─'*70}")
    
    # Step 3: Verify all pairs in each complete chain
    valid_algs = []
    for chain in all_candidate_chains:
        if len(chain) == len(cluster_species):
            if verify_alg_chain(chain, cluster_species, sig_chrom_pairs):
                valid_algs.append(chain)
    
    print(f"Valid ALGs: {len(valid_algs)}")
    
    # Step 4: Assign ALG IDs
    alg_assignments = {}
    
    for alg_id, chain in enumerate(valid_algs, start=1):
        print(f"\nALG {alg_id}:")
        for species in cluster_species:
            if species in chain:
                chrom = chain[species]
                print(f"  {species}: {chrom}")
                
                if species not in alg_assignments:
                    alg_assignments[species] = {}
                if chrom not in alg_assignments[species]:
                    alg_assignments[species][chrom] = []
                
                alg_assignments[species][chrom].append(alg_id)
    
    # Report chromosomes in multiple ALGs
    multi_alg_chroms = [
        (species, chrom, alg_list)
        for species, chr_dict in alg_assignments.items()
        for chrom, alg_list in chr_dict.items()
        if len(alg_list) > 1
    ]
    
    if multi_alg_chroms:
        print(f"\n📍 Chromosomes in multiple ALGs:")
        for species, chrom, alg_list in multi_alg_chroms:
            print(f"  {species} {chrom}: ALG{', ALG'.join(map(str, alg_list))}")
    
    return alg_assignments


def compute_algs_full_pipeline(species_maps, species_names, fishers_multi_results,
                                 similarity_threshold=0.15, min_orthologs=10):
    """
    Complete Ancestral Linkage Group (ALG) identification pipeline.
    
    Orchestrates the full ALG discovery workflow:
    1. Computes pairwise synteny similarity between all species
    2. Clusters species by synteny similarity
    3. For each cluster, identifies ALGs by building and verifying
       chromosome chains across all species in the cluster
    
    Parameters:
    -----------
    species_maps : list of dict
        List of physical chromosome maps for each species, in the same
        order as species_names. Each map is the output of
        synteny_map_creator: {chromosome: [(og_id, protein_id, start, end), ...]}
    species_names : list of str
        List of species IDs in the same order as species_maps.
        Example: ['CLA', 'EMU', 'BFL', 'RES', 'HCA', 'BMI']
    fishers_multi_results : dict
        Output of fishers_multi. Format:
        {og_id: {'positions': [...], 'significant_segments': [(i,j), ...]}}
    similarity_threshold : float
        Minimum synteny similarity for two species to be placed in the
        same cluster (default: 0.15). Passed to cluster_species_by_synteny.
    min_orthologs : int
        Minimum number of significant orthologs required for a chromosome
        pair to be considered a valid ALG link (default: 10). Passed to
        precompute_significant_chrom_pairs.
    
    Returns:
    --------
    dict with the following keys:
        'similarity_matrix' : np.ndarray
            N x N pairwise synteny similarity matrix
        'clusters' : dict
            {cluster_id: [species_names]} grouping of species
        'cluster_labels' : np.ndarray
            Per-species cluster assignment array
        'alg_assignments' : dict
            Nested dict of ALG assignments per species per chromosome.
            Format: {species: {chromosome: [alg_id, ...]}}
            Example: {'CLA': {'CLA1': ['C2_ALG1'], 'CLA5': ['C2_ALG3', 'C2_ALG6']}}
            ALG IDs are prefixed with cluster ID: C1_ALG1 = cluster 1, ALG 1
        'clusters_info' : dict
            {cluster_id: [species_names]} - identical to 'clusters',
            kept separately for use by downstream plotting functions
    """
    
    print("\n" + "="*70)
    print("ANCESTRAL LINKAGE GROUP (ALG) IDENTIFICATION")
    print("="*70)
    
    similarity = compute_synteny_similarity(
        species_maps, species_names, fishers_multi_results
    )
    
    clusters, cluster_labels = cluster_species_by_synteny(
        similarity, species_names, threshold=similarity_threshold
    )
    
    all_alg_assignments = {}
    all_clusters_info = {}
    
    for cluster_id, cluster_species in sorted(clusters.items()):
        alg_assignments = identify_algs_for_cluster(
            cluster_species, species_maps, species_names,
            fishers_multi_results, min_orthologs=min_orthologs
        )
        
        for species, chrom_to_algs in alg_assignments.items():
            if species not in all_alg_assignments:
                all_alg_assignments[species] = {}
            
            for chrom, alg_list in chrom_to_algs.items():
                prefixed_algs = [f"C{cluster_id}_ALG{alg}" for alg in alg_list]
                all_alg_assignments[species][chrom] = prefixed_algs
        
        all_clusters_info[cluster_id] = cluster_species
    
    return {
        'similarity_matrix': similarity,
        'clusters': clusters,
        'alg_assignments': all_alg_assignments,
        'cluster_labels': cluster_labels,
        'clusters_info': all_clusters_info
    }

## Ribbon Plot Creation

def plot_synteny_ribbons(comparison_map, sp1_map, sp2_map, 
                         species1='Species1', species2='Species2',
                         figsize=(20, 12), ribbon_alpha=0.3, curve_style='bezier'):
    """
    Create a synteny ribbon plot showing chromosome relationships between two species.
    Each ortholog pair is shown as an individual ribbon.
    Species 2 chromosomes are rearranged to maximize visual alignment with species 1.
    
    Parameters:
    -----------
    comparison_map : dict
        Dictionary with OG IDs as keys and list of tuples as values
    sp1_map : dict
        Dictionary with chromosome as key and list of gene tuples for species 1
    sp2_map : dict
        Dictionary with chromosome as key and list of gene tuples for species 2
    species1 : str
        Name of first species
    species2 : str
        Name of second species
    figsize : tuple
        Figure size
    ribbon_alpha : float
        Transparency of ribbons (default: 0.3)
    curve_style : str
        'bezier' for smooth curves, 'straight' for straight lines (default: 'bezier')
    """
    
    # Bezier curve function
    def get_bezier_curve(x0, y0, x1, y1, n_points=100):
        """Generate Bezier curve between two points."""
        if curve_style == 'straight':
            return [x0, x1], [y0, y1]
        
        # Control points at 1/3 and 2/3 of vertical distance
        cy0 = y0 - (y0 - y1) * 0.33
        cy1 = y0 - (y0 - y1) * 0.67
        
        t = np.linspace(0, 1, n_points)
        
        # Cubic Bezier formula
        x = (1-t)**3 * x0 + 3*(1-t)**2*t * x0 + 3*(1-t)*t**2 * x1 + t**3 * x1
        y = (1-t)**3 * y0 + 3*(1-t)**2*t * cy0 + 3*(1-t)*t**2 * cy1 + t**3 * y1
        
        return x, y
    
    # Get chromosome lengths
    sp1_chrom_lengths = {}
    for chrom, genes in sp1_map.items():
        max_end = max(gene[3] for gene in genes)
        sp1_chrom_lengths[chrom] = max_end
    
    sp2_chrom_lengths = {}
    for chrom, genes in sp2_map.items():
        max_end = max(gene[3] for gene in genes)
        sp2_chrom_lengths[chrom] = max_end
    
    # Get sorted chromosomes that have orthologs
    sp1_chroms_with_orthologs = sorted(set(v[0][1] for v in comparison_map.values()), key=chrom_sort_key)
    sp2_chroms_with_orthologs = sorted(set(v[1][1] for v in comparison_map.values()), key=chrom_sort_key)
    
    # Filter to only chromosomes with orthologs
    sp1_chroms = [c for c in sp1_chroms_with_orthologs if c in sp1_chrom_lengths]
    sp2_chroms_unsorted = [c for c in sp2_chroms_with_orthologs if c in sp2_chrom_lengths]
    
    # Count orthologs between each chromosome pair
    chrom_pair_counts = defaultdict(int)
    for og_id, orthologs in comparison_map.items():
        sp1_chrom = orthologs[0][1]
        sp2_chrom = orthologs[1][1]
        if sp1_chrom in sp1_chroms and sp2_chrom in sp2_chroms_unsorted:
            chrom_pair_counts[(sp1_chrom, sp2_chrom)] += 1
    
    # Rearrange sp2 chromosomes for better alignment
    sp2_chroms = []
    used_sp2_chroms = set()
    
    for sp1_chrom in sp1_chroms:
        matches = [(sp2_chrom, chrom_pair_counts[(sp1_chrom, sp2_chrom)]) 
                   for sp2_chrom in sp2_chroms_unsorted 
                   if sp2_chrom not in used_sp2_chroms]
        matches.sort(key=lambda x: x[1], reverse=True)
        
        for sp2_chrom, count in matches:
            if sp2_chrom not in used_sp2_chroms:
                sp2_chroms.append(sp2_chrom)
                used_sp2_chroms.add(sp2_chrom)
                print(f"Aligned {species1} {sp1_chrom} with {species2} {sp2_chrom} ({count} orthologs)")
                break
    
    # Add any remaining sp2 chromosomes that weren't matched
    for sp2_chrom in sp2_chroms_unsorted:
        if sp2_chrom not in used_sp2_chroms:
            sp2_chroms.append(sp2_chrom)
            print(f"Added unmatched {species2} {sp2_chrom} at end")
    
    # Map chromosomes to colors (for ribbons only)
    sp1_chrom_to_color = {}
    for idx, chrom in enumerate(sp1_chroms):
        sp1_chrom_to_color[chrom] = custom_colors[idx % len(custom_colors)]
    
    # Calculate total genome lengths
    sp1_total_length = sum(sp1_chrom_lengths[c] for c in sp1_chroms)
    sp2_total_length = sum(sp2_chrom_lengths[c] for c in sp2_chroms)

    # Use a standard horizontal length
    standard_length = 1000000
    gap_size = standard_length * 0.005  # Small gap (0.5% of total length)

    # Calculate cumulative positions for species 1 (top) with gaps
    sp1_cumulative = {}
    cumsum = 0
    for chrom in sp1_chroms:
        sp1_cumulative[chrom] = cumsum
        cumsum += (sp1_chrom_lengths[chrom] / sp1_total_length) * standard_length
        cumsum += gap_size  # Add gap after each chromosome

    # Calculate cumulative positions for species 2 (bottom) with gaps
    sp2_cumulative = {}
    cumsum = 0
    for chrom in sp2_chroms:
        sp2_cumulative[chrom] = cumsum
        cumsum += (sp2_chrom_lengths[chrom] / sp2_total_length) * standard_length
        cumsum += gap_size  # Add gap after each chromosome
    
    # Get gene counts per chromosome for positioning
    sp1_gene_counts = {}
    sp2_gene_counts = {}
    
    for chrom in sp1_chroms:
        sp1_gene_counts[chrom] = len(sp1_map[chrom])
    
    for chrom in sp2_chroms:
        sp2_gene_counts[chrom] = len(sp2_map[chrom])
    
    # Create plot
    fig, ax = plt.subplots(figsize=figsize)
    
    # Draw individual gene ribbons with Bezier curves
    for og_id, orthologs in comparison_map.items():
        sp1_chrom = orthologs[0][1]
        sp1_pos = orthologs[0][2]
        sp2_chrom = orthologs[1][1]
        sp2_pos = orthologs[1][2]
        
        if sp1_chrom in sp1_chroms and sp2_chrom in sp2_chroms and sp1_pos > 0 and sp2_pos > 0:
            # Calculate normalized positions based on ordinal position
            sp1_chrom_width = (sp1_chrom_lengths[sp1_chrom] / sp1_total_length) * standard_length
            sp2_chrom_width = (sp2_chrom_lengths[sp2_chrom] / sp2_total_length) * standard_length
            
            # Position within chromosome as fraction
            sp1_fraction = sp1_pos / sp1_gene_counts[sp1_chrom]
            sp2_fraction = sp2_pos / sp2_gene_counts[sp2_chrom]
            
            # Calculate absolute position
            sp1_x = sp1_cumulative[sp1_chrom] + sp1_fraction * sp1_chrom_width
            sp2_x = sp2_cumulative[sp2_chrom] + sp2_fraction * sp2_chrom_width
            
            # Get Bezier curve points
            curve_x, curve_y = get_bezier_curve(sp1_x, 1, sp2_x, 0)
            
            color = sp1_chrom_to_color[sp1_chrom]
            ax.plot(curve_x, curve_y, 
                   color=color, alpha=ribbon_alpha, linewidth=1.5, 
                   solid_capstyle='round')
    
    # Draw chromosome lines for species 1 (top) - thin black lines
    for chrom in sp1_chroms:
        start = sp1_cumulative[chrom]
        width = (sp1_chrom_lengths[chrom] / sp1_total_length) * standard_length
        
        # Draw as simple black line
        ax.plot([start, start + width], [1, 1], 
               color='black', linewidth=3, solid_capstyle='butt')
        
        # Add chromosome labels
        mid = start + width/2
        ax.text(mid, 1 + 0.03, str(chrom), 
               ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    # Draw chromosome lines for species 2 (bottom) - thin black lines
    for chrom in sp2_chroms:
        start = sp2_cumulative[chrom]
        width = (sp2_chrom_lengths[chrom] / sp2_total_length) * standard_length
        
        # Draw as simple black line
        ax.plot([start, start + width], [0, 0], 
               color='black', linewidth=3, solid_capstyle='butt')
        
        # Add chromosome labels
        mid = start + width/2
        ax.text(mid, -0.03, str(chrom),
               ha='center', va='top', fontsize=10, fontweight='bold')
    
    # Add species labels
    ax.text(standard_length/2, 1.15, species1, ha='center', va='bottom', 
           fontsize=14, fontweight='bold')
    ax.text(standard_length/2, -0.15, species2, ha='center', va='top',
           fontsize=14, fontweight='bold')
    
    # Set limits and clean up
    max_x = max(
        max(
            sp1_cumulative[chrom] +
            (sp1_chrom_lengths[chrom] / sp1_total_length) * standard_length
            for chrom in sp1_chroms
        ),
        max(
            sp2_cumulative[chrom] +
            (sp2_chrom_lengths[chrom] / sp2_total_length) * standard_length
            for chrom in sp2_chroms
        )
    )

    ax.set_xlim(-standard_length * 0.05, max_x * 1.02)
    ax.set_ylim(-0.2, 1.2)
    ax.axis('off')
    
    plt.tight_layout()
    return fig, ax

def plot_synteny_ribbons_multi(filtered_map, species_maps, species_names,
                                alg_results=None, ribbon_alpha=0.3,
                                figsize=(25, 15), curve_style='bezier'):

    """Create a multi-species synteny ribbon plot for N species.

    Draws Bezier curve ribbons connecting orthologous genes across all
    species. Species are arranged vertically, with chromosomes displayed
    as horizontal black lines. If ALG results are provided, ribbons are
    colored by ALG identity and cross-cluster ribbons are suppressed.
    Otherwise, ribbons are colored by the first species' chromosomes.

    Chromosome ordering: the first species in each cluster is sorted
    nominally; subsequent species are ordered to minimize ribbon crossings
    with the species above them.

    Parameters:
    -----------
    filtered_map : dict
        Output of fishers_multi.
        Format: {og_id: {'positions': [...], 'significant_segments': [(i,j), ...]}}
    species_maps : list of dict
        Physical chromosome maps for each species, in the same order as
        species_names. Output of synteny_map_creator.
    species_names : list of str
        Species IDs in display order (top to bottom in the plot)
    alg_results : dict or None
        Output of compute_algs_full_pipeline. If provided, enables ALG-based
        coloring and per-cluster chromosome ordering. If None, falls back to
        first-species chromosome coloring (default: None)
    ribbon_alpha : float
        Transparency of ribbon lines (default: 0.3)
    figsize : tuple
        Figure dimensions as (width, height) in inches (default: (25, 15))
    curve_style : str
        'bezier' for smooth cubic Bezier curves, 'straight' for straight
        lines (default: 'bezier')

    Returns:
    --------
    fig : matplotlib.figure.Figure
    ax : matplotlib.axes.Axes
    """
    
    n_species = len(species_names)

    # ===== ALG SETUP =====
    if alg_results is not None:
        clusters_info = alg_results.get('clusters_info', {})
        alg_assignments = alg_results.get('alg_assignments', {})

        all_algs = set()
        for sp_algs in alg_assignments.values():
            for alg_list in sp_algs.values():
                if isinstance(alg_list, list):
                    all_algs.update(alg_list)
                else:
                    all_algs.add(alg_list)

        alg_colors = {}
        if all_algs:
            for idx, alg in enumerate(sorted(all_algs)):
                alg_colors[alg] = custom_colors[idx % len(custom_colors)]

        def get_alg_for_segment(og_id, sp1_name, chrom1, sp2_name, chrom2):
            if sp1_name not in alg_assignments or sp2_name not in alg_assignments:
                return None
            if chrom1 not in alg_assignments[sp1_name]:
                return None
            if chrom2 not in alg_assignments[sp2_name]:
                return None
            alg_list1 = alg_assignments[sp1_name][chrom1]
            alg_list2 = alg_assignments[sp2_name][chrom2]
            if not isinstance(alg_list1, list):
                alg_list1 = [alg_list1]
            if not isinstance(alg_list2, list):
                alg_list2 = [alg_list2]
            common = set(alg_list1) & set(alg_list2)
            if common:
                return sorted(common)[0]
            return None

        same_cluster_pairs = set()
        for cluster_species in clusters_info.values():
            for i in range(len(cluster_species)):
                for j in range(len(cluster_species)):
                    if i != j:
                        same_cluster_pairs.add((cluster_species[i], cluster_species[j]))

        # Build mapping: species_name -> cluster_id
        species_to_cluster = {}
        for cluster_id, cluster_species in clusters_info.items():
            for sp in cluster_species:
                species_to_cluster[sp] = cluster_id
    else:
        clusters_info = {}
        alg_assignments = {}
        alg_colors = {}
        same_cluster_pairs = None
        species_to_cluster = {}

    # ===== CHROMOSOME ORDERING (PER CLUSTER) =====
    # First, determine which species belong to which cluster and in what order
    # Species order within each cluster follows their order in species_names

    # Build cluster-aware species ordering
    if clusters_info:
        # Get ordered species per cluster (preserving species_names order)
        cluster_ordered_species = {}
        for cluster_id, cluster_species in clusters_info.items():
            cluster_ordered_species[cluster_id] = [
                sp for sp in species_names if sp in cluster_species
            ]
    else:
        # No clusters - treat all as one cluster
        cluster_ordered_species = {1: species_names}

    all_species_chroms = [None] * n_species

    for cluster_id, ordered_species in cluster_ordered_species.items():
        # First species in this cluster: sorted order
        first_sp = ordered_species[0]
        first_sp_idx = species_names.index(first_sp)
        first_chroms = sorted(species_maps[first_sp_idx].keys(), key=chrom_sort_key)
        all_species_chroms[first_sp_idx] = first_chroms

        # Subsequent species in cluster: ordered by ortholog counts with previous
        for i in range(1, len(ordered_species)):
            prev_sp = ordered_species[i - 1]
            curr_sp = ordered_species[i]
            prev_sp_idx = species_names.index(prev_sp)
            curr_sp_idx = species_names.index(curr_sp)

            prev_chroms = all_species_chroms[prev_sp_idx]
            current_chroms_unsorted = sorted(species_maps[curr_sp_idx].keys(), key=chrom_sort_key)

            chrom_pair_counts = defaultdict(int)
            for og_id, data in filtered_map.items():
                positions = data['positions']
                if positions[prev_sp_idx] is not None and positions[curr_sp_idx] is not None:
                    prev_chrom = positions[prev_sp_idx][1]
                    curr_chrom = positions[curr_sp_idx][1]
                    if prev_chrom in prev_chroms and curr_chrom in current_chroms_unsorted:
                        chrom_pair_counts[(prev_chrom, curr_chrom)] += 1

            ordered_chroms = []
            used_chroms = set()

            for prev_chrom in prev_chroms:
                matches = [
                    (curr_chrom, chrom_pair_counts[(prev_chrom, curr_chrom)])
                    for curr_chrom in current_chroms_unsorted
                    if curr_chrom not in used_chroms
                ]
                matches.sort(key=lambda x: x[1], reverse=True)

                for curr_chrom, count in matches:
                    if curr_chrom not in used_chroms:
                        ordered_chroms.append(curr_chrom)
                        used_chroms.add(curr_chrom)
                        break

            for curr_chrom in current_chroms_unsorted:
                if curr_chrom not in used_chroms:
                    ordered_chroms.append(curr_chrom)

            all_species_chroms[curr_sp_idx] = ordered_chroms

    # ===== CHROMOSOME LENGTHS =====
    species_chr_lengths = []
    for sp_idx, species_map in enumerate(species_maps):
        chr_lengths = {}
        for chrom in all_species_chroms[sp_idx]:
            if chrom in species_map:
                max_end = max(gene[3] for gene in species_map[chrom])
                chr_lengths[chrom] = max_end
        species_chr_lengths.append(chr_lengths)

    # ===== COLOR MAPPING FOR NON-ALG MODE =====
    sp1_chrom_to_color = {}
    for idx, chrom in enumerate(all_species_chroms[0]):
        sp1_chrom_to_color[chrom] = custom_colors[idx % len(custom_colors)]

    # ===== NORMALIZED CUMULATIVE POSITIONS WITH GAPS =====
    standard_length = 1000000
    gap_size = standard_length * 0.005

    species_cumulative = []
    species_total_lengths = []

    for sp_idx in range(n_species):
        total_length = sum(species_chr_lengths[sp_idx].values())
        species_total_lengths.append(total_length)

        cumulative = {}
        cumsum = 0
        for chrom in all_species_chroms[sp_idx]:
            cumulative[chrom] = cumsum
            cumsum += (species_chr_lengths[sp_idx][chrom] / total_length) * standard_length
            cumsum += gap_size

        species_cumulative.append(cumulative)

    # ===== GENE COUNTS =====
    species_gene_counts = []
    for sp_idx in range(n_species):
        gene_counts = {}
        for chrom in all_species_chroms[sp_idx]:
            if chrom in species_maps[sp_idx]:
                gene_counts[chrom] = len(species_maps[sp_idx][chrom])
        species_gene_counts.append(gene_counts)

    # ===== BEZIER CURVE FUNCTION =====
    def get_bezier_curve(points, n_points=100):
        if curve_style == 'straight' or len(points) < 2:
            xs, ys = zip(*points)
            return list(xs), list(ys)

        all_x = []
        all_y = []

        for i in range(len(points) - 1):
            x0, y0 = points[i]
            x1, y1 = points[i + 1]

            cy0 = y0 - (y0 - y1) * 0.33
            cy1 = y0 - (y0 - y1) * 0.67

            t = np.linspace(0, 1, n_points // max(1, len(points) - 1))

            x = (1-t)**3 * x0 + 3*(1-t)**2*t * x0 + 3*(1-t)*t**2 * x1 + t**3 * x1
            y = (1-t)**3 * y0 + 3*(1-t)**2*t * cy0 + 3*(1-t)*t**2 * cy1 + t**3 * y1

            all_x.extend(x)
            all_y.extend(y)

        return all_x, all_y

    # ===== CREATE PLOT =====
    fig, ax = plt.subplots(figsize=figsize)
    y_positions = np.linspace(1, 0, n_species)

    # ===== PRE-COMPUTE CLUSTER BOUNDARY Y POSITIONS =====
    # Store ALL boundary y positions (one per cluster boundary)
    cluster_boundary_ys = []
    if alg_results is not None and clusters_info:
        for sp_idx in range(1, n_species):
            sp_prev = species_names[sp_idx - 1]
            sp_curr = species_names[sp_idx]
            different_cluster = all(
                not (sp_prev in cluster_species and sp_curr in cluster_species)
                for cluster_species in clusters_info.values()
            )
            if different_cluster:
                cluster_boundary_ys.append(
                    (y_positions[sp_idx - 1] + y_positions[sp_idx]) / 2
                )

    # ===== DRAW RIBBONS =====
    ribbons_drawn = 0

    for og_id, data in filtered_map.items():
        positions = data['positions']
        significant_segments = data['significant_segments']

        adjacent_segments = [
            (i, j) for (i, j) in significant_segments
            if j == i + 1
        ]

        for (seg_i, seg_j) in adjacent_segments:
            pos1 = positions[seg_i]
            pos2 = positions[seg_j]

            if pos1 is None or pos2 is None:
                continue

            species1_name, chrom1, pos_idx1 = pos1
            species2_name, chrom2, pos_idx2 = pos2

            if same_cluster_pairs is not None:
                if (species1_name, species2_name) not in same_cluster_pairs:
                    continue

            chr_width1 = (species_chr_lengths[seg_i][chrom1] / species_total_lengths[seg_i]) * standard_length
            fraction1 = pos_idx1 / species_gene_counts[seg_i][chrom1]
            x1 = species_cumulative[seg_i][chrom1] + fraction1 * chr_width1
            y1 = y_positions[seg_i]

            chr_width2 = (species_chr_lengths[seg_j][chrom2] / species_total_lengths[seg_j]) * standard_length
            fraction2 = pos_idx2 / species_gene_counts[seg_j][chrom2]
            x2 = species_cumulative[seg_j][chrom2] + fraction2 * chr_width2
            y2 = y_positions[seg_j]

            if alg_results is not None and alg_colors:
                alg_id = get_alg_for_segment(og_id, species1_name, chrom1, species2_name, chrom2)
                if alg_id is None or alg_id not in alg_colors:
                    continue
                ribbon_color = alg_colors[alg_id]
            else:
                first_valid_pos = next((p for p in positions if p is not None), None)
                if first_valid_pos and first_valid_pos[1] in sp1_chrom_to_color:
                    ribbon_color = sp1_chrom_to_color[first_valid_pos[1]]
                else:
                    ribbon_color = 'gray'

            segment_points = [(x1, y1), (x2, y2)]
            curve_x, curve_y = get_bezier_curve(segment_points)

            ax.plot(
                curve_x, curve_y,
                color=ribbon_color,
                alpha=ribbon_alpha,
                linewidth=1.5,
                solid_capstyle='round',
                zorder=1
            )
            ribbons_drawn += 1

    print(f"\nTotal ribbon segments drawn: {ribbons_drawn}")

    # ===== DRAW CHROMOSOME LINES =====
    for sp_idx in range(n_species):
        y_pos = y_positions[sp_idx]

        for chrom in all_species_chroms[sp_idx]:
            start = species_cumulative[sp_idx][chrom]
            width = (species_chr_lengths[sp_idx][chrom] / species_total_lengths[sp_idx]) * standard_length

            ax.plot([start, start + width], [y_pos, y_pos],
                   color='black', linewidth=3, solid_capstyle='butt', zorder=2)

            mid = start + width / 2
            if sp_idx == 0:
                ax.text(mid, y_pos + 0.03, str(chrom),
                       ha='center', va='bottom', fontsize=10, fontweight='bold')
            else:
                ax.text(mid, y_pos - 0.03, str(chrom),
                       ha='center', va='top', fontsize=10, fontweight='bold')

    # ===== SPECIES LABELS =====
    for sp_idx in range(n_species):
        y_pos = y_positions[sp_idx]
        ax.text(-standard_length * 0.05, y_pos, species_names[sp_idx],
               ha='right', va='center', fontsize=14, fontweight='bold')

    # ===== CLUSTER BOUNDARY MARKERS =====
    for boundary_y in cluster_boundary_ys:
        ax.plot([0, standard_length], [boundary_y, boundary_y],
               color='red', linestyle='--', linewidth=2, alpha=0.7, zorder=3)
        ax.text(standard_length * 1.01, boundary_y, 'Cluster Boundary',
               va='center', fontsize=10, color='red', style='italic')

    # ===== ALG LEGEND =====
    if alg_results is not None and alg_colors:

        cluster_bottom_species = {}
        for cluster_id, cluster_species in clusters_info.items():
            species_y = {sp: y_positions[species_names.index(sp)] for sp in cluster_species}
            bottom_sp = min(species_y, key=species_y.get)
            bottom_sp_idx = species_names.index(bottom_sp)
            cluster_bottom_species[cluster_id] = (bottom_sp, bottom_sp_idx)

        rect_height = 0.012
        rect_gap = 0.005
        label_fontsize = 9

        for cluster_id, (bottom_sp, bottom_sp_idx) in cluster_bottom_species.items():
            y_pos = y_positions[bottom_sp_idx]

            # Find the boundary BELOW this cluster (if any)
            # The boundary below this cluster is the one with y value just below y_pos
            boundary_below = None
            for boundary_y in cluster_boundary_ys:
                if boundary_y < y_pos:
                    # This boundary is below this cluster's bottom species
                    if boundary_below is None or boundary_y > boundary_below:
                        boundary_below = boundary_y

            # Calculate max stack depth
            max_stack = 1
            if bottom_sp in alg_assignments:
                for chrom in all_species_chroms[bottom_sp_idx]:
                    if chrom in alg_assignments[bottom_sp]:
                        chrom_algs = alg_assignments[bottom_sp][chrom]
                        stack = len(chrom_algs) if isinstance(chrom_algs, list) else 1
                        max_stack = max(max_stack, stack)

            # Start drawing below chromosome labels
            legend_y_start = y_pos - 0.07

            # Clamp above boundary below this cluster
            if boundary_below is not None:
                deepest_point = legend_y_start - (max_stack - 1) * (rect_height + rect_gap) - rect_height
                if deepest_point < boundary_below + 0.01:
                    shift = (boundary_below + 0.01) - deepest_point
                    legend_y_start = legend_y_start + shift

            for chrom in all_species_chroms[bottom_sp_idx]:
                start = species_cumulative[bottom_sp_idx][chrom]
                width = (species_chr_lengths[bottom_sp_idx][chrom] / species_total_lengths[bottom_sp_idx]) * standard_length

                if bottom_sp in alg_assignments and chrom in alg_assignments[bottom_sp]:
                    chrom_algs = alg_assignments[bottom_sp][chrom]
                    if not isinstance(chrom_algs, list):
                        chrom_algs = [chrom_algs]

                    chrom_algs = sorted(chrom_algs)

                    for alg_rank, alg_id in enumerate(chrom_algs):
                        if alg_id not in alg_colors:
                            continue

                        rect_y = legend_y_start - alg_rank * (rect_height + rect_gap)
                        rect_width = width * 0.425
                        rect_x = start + (width - rect_width) / 2

                        rect = plt.Rectangle(
                            (rect_x, rect_y - rect_height),
                            rect_width,
                            rect_height,
                            facecolor=alg_colors[alg_id],
                            edgecolor='none',
                            zorder=3
                        )
                        ax.add_patch(rect)

                        alg_label = alg_id.split('_')[-1] if '_' in alg_id else str(alg_id)

                        ax.text(
                            rect_x + rect_width + width * 0.05,
                            rect_y - rect_height / 2,
                            alg_label,
                            ha='left',
                            va='center',
                            fontsize=label_fontsize,
                            fontstyle='italic',
                            color='black',
                            zorder=4
                        )

    # ===== SET AXIS LIMITS =====
    max_x = max(
        species_cumulative[sp_idx][chrom] +
        (species_chr_lengths[sp_idx][chrom] / species_total_lengths[sp_idx]) * standard_length
        for sp_idx in range(n_species)
        for chrom in all_species_chroms[sp_idx]
    )

    max_algs_per_chrom = 1
    if alg_results is not None:
        for sp_algs in alg_assignments.values():
            for alg_list in sp_algs.values():
                if isinstance(alg_list, list):
                    max_algs_per_chrom = max(max_algs_per_chrom, len(alg_list))

    legend_depth = 0.07 + max_algs_per_chrom * (0.012 + 0.005)

    ax.set_xlim(-standard_length * 0.08, max_x * 1.05)
    ax.set_ylim(-legend_depth - 0.05, 1.2)
    ax.axis('off')

    plt.tight_layout()

    return fig, ax

## EXTRAS

### Simakov's gene families analysis

def simakov_gene_families(triangle_species, rest_species, rbh_dir, tsv_dir, mbh=False):

    """
    Recreate Simakov's method for identifying conserved gene families across species.
    
    This function implements a stringent orthology inference approach:
    - For triangle mode (mbh=False): Starts with 3 core species forming mutual best hits (MBH),
      then adds additional species by requiring they have MBH with ALL core species
    - For full MBH mode (mbh=True): Requires mutual best hits across ALL species simultaneously
    
    Parameters:
    -----------
    triangle_species : list
        List of 3 species IDs for the core triangle (e.g., ['Bfloridae', 'Pmaximus', 'Emuelleri'])
        Only used when mbh=False
    rest_species : list
        List of additional species IDs (e.g., ['Hvulgaris', 'Resculentum'])
    rbh_dir : str
        Directory containing RBH TSV files
        Files should be named: 'Species1__RBH__Species2.tsv'
    tsv_dir : str
        Directory containing chromosome position TSV files
        Files should be named: '{species}.tsv'
    mbh : bool, default=False
        If False: Use triangle-based approach (3 core species, then add others)
        If True: Require MBH across ALL species simultaneously (more stringent)
    
    Returns:
    --------
    gene_families : list of tuples
        Each tuple contains chromosomes where orthologous genes are found
        Format: (sp1_chr, sp2_chr, sp3_chr, ...) in order of species list
    gene_family_details : list of dicts
        Detailed information for each gene family including protein IDs
    
    Examples:
    ---------
    # Triangle mode (Simakov's original approach)
    >>> families, details = simakov_gene_families(
            triangle_species=['Bfloridae', 'Pmaximus', 'Emuelleri'],
            rest_species=['Hvulgaris', 'Resculentum'],
            rbh_dir
            tsv_dir
            mbh=False
        )
    
    # Full MBH mode (all species must be MBH with each other)
    >>> families, details = simakov_gene_families(
            triangle_species=[],  # Not used in MBH mode
            rest_species=['Bfloridae', 'Pmaximus', 'Emuelleri', 'Hvulgaris', 'Resculentum'],
            rbh_dir
            tsv_dir
            mbh=True
        )
    """
    
    def load_rbh(sp1, sp2, rbh_dir):
        """Load RBH file and return dict mapping sp1 proteins to sp2 proteins."""
        # Try both possible file name orders
        file1 = Path(rbh_dir) / f"{sp1}__RBH__{sp2}.tsv"
        file2 = Path(rbh_dir) / f"{sp2}__RBH__{sp1}.tsv"
        
        if file1.exists():
            df = pd.read_csv(file1, sep='\t')
            return dict(zip(df[sp1], df[sp2]))
        elif file2.exists():
            df = pd.read_csv(file2, sep='\t')
            return dict(zip(df[sp1], df[sp2]))
        else:
            raise FileNotFoundError(f"RBH file not found for {sp1} vs {sp2}")
    
    def load_chromosome_map(species, tsv_dir):
        """Load chromosome positions for a species."""
        tsv_file = Path(tsv_dir) / f"{species}.tsv"
        if not tsv_file.exists():
            raise FileNotFoundError(f"Chromosome file not found: {tsv_file}")
        
        df = pd.read_csv(tsv_file, sep='\t')
        # Create mapping: protein_id -> chromosome
        return dict(zip(df['ProteinID'], df['Chr/Scaffold']))
    
    print("="*100)
    if mbh:
        all_species = rest_species
        print(f"FULL MBH MODE: Finding {len(all_species)}-way mutual best hits")
        print(f"Species: {', '.join(all_species)}")
    else:
        all_species = triangle_species + rest_species
        print(f"TRIANGLE MODE: 3-way core + {len(rest_species)} additional species")
        print(f"Core triangle: {', '.join(triangle_species)}")
        print(f"Additional species: {', '.join(rest_species)}")
    print("="*100 + "\n")
    
    # Load chromosome maps for all species
    print("Loading chromosome maps...")
    chrom_maps = {}
    for species in all_species:
        chrom_maps[species] = load_chromosome_map(species, tsv_dir)
        print(f"  {species}: {len(chrom_maps[species])} proteins")
    
    if mbh:
        # ===== FULL MBH MODE =====
        print(f"\nLoading all pairwise RBH files...")
        
        # Load all pairwise RBH relationships
        rbh_data = {}
        for i, sp1 in enumerate(all_species):
            for sp2 in all_species[i+1:]:
                rbh_data[(sp1, sp2)] = load_rbh(sp1, sp2, rbh_dir)
                print(f"  {sp1} <-> {sp2}: {len(rbh_data[(sp1, sp2)])} RBH pairs")
        
        print(f"\nFinding {len(all_species)}-way mutual best hits...")
        
        # Start with proteins from first species
        gene_families = []
        gene_family_details = []
        
        sp0 = all_species[0]
        
        for protein0 in chrom_maps[sp0].keys():
            # Track proteins across all species
            proteins_in_family = {sp0: protein0}
            
            # Try to extend to all other species
            valid_family = True
            
            # For each subsequent species, find the MBH
            for i, sp_current in enumerate(all_species[1:], 1):
                # Get the protein from previous species
                sp_prev = all_species[i-1]
                protein_prev = proteins_in_family.get(sp_prev)
                
                if not protein_prev:
                    valid_family = False
                    break
                
                # Find RBH between previous species and current species
                pair_key = (sp_prev, sp_current) if (sp_prev, sp_current) in rbh_data else (sp_current, sp_prev)
                
                if pair_key not in rbh_data:
                    valid_family = False
                    break
                
                rbh_map = rbh_data[pair_key]
                
                if pair_key[0] == sp_prev:
                    protein_current = rbh_map.get(protein_prev)
                else:
                    # Reverse mapping
                    reverse_map = {v: k for k, v in rbh_map.items()}
                    protein_current = reverse_map.get(protein_prev)
                
                if not protein_current:
                    valid_family = False
                    break
                
                proteins_in_family[sp_current] = protein_current
            
            # Verify this is a complete MBH across ALL species
            if valid_family and len(proteins_in_family) == len(all_species):
                # Additional check: verify ALL pairwise RBHs exist
                all_pairs_valid = True
                for i, sp1 in enumerate(all_species):
                    for sp2 in all_species[i+1:]:
                        pair_key = (sp1, sp2) if (sp1, sp2) in rbh_data else (sp2, sp1)
                        rbh_map = rbh_data[pair_key]
                        
                        prot1 = proteins_in_family[sp1]
                        prot2 = proteins_in_family[sp2]
                        
                        # Check if this pair exists in RBH
                        if pair_key[0] == sp1:
                            if rbh_map.get(prot1) != prot2:
                                all_pairs_valid = False
                                break
                        else:
                            reverse_map = {v: k for k, v in rbh_map.items()}
                            if reverse_map.get(prot1) != prot2:
                                all_pairs_valid = False
                                break
                    if not all_pairs_valid:
                        break
                
                if all_pairs_valid:
                    # Get chromosomes for all proteins
                    try:
                        chromosomes = tuple(
                            chrom_maps[sp][proteins_in_family[sp]] 
                            for sp in all_species
                        )
                    except KeyError as e:
                        print(f'Missing key: {e}')
                        continue
                        
                    gene_families.append(chromosomes)
                    gene_family_details.append({
                        'proteins': proteins_in_family,
                        'chromosomes': chromosomes
                    })
        
        print(f"\nFound {len(gene_families)} complete {len(all_species)}-way MBH families")
    
    else:
        # ===== TRIANGLE MODE =====
        if len(triangle_species) != 3:
            raise ValueError("Triangle mode requires exactly 3 core species")
        
        sp1, sp2, sp3 = triangle_species
        
        # Step 1: Build 3-way MBH core
        print(f"\nStep 1: Building 3-way MBH core ({sp1}, {sp2}, {sp3})...")
        
        # Load RBH files for triangle
        rbh_12 = load_rbh(sp1, sp2, rbh_dir)
        rbh_23 = load_rbh(sp2, sp3, rbh_dir)
        rbh_13 = load_rbh(sp1, sp3, rbh_dir)
        
        print(f"  {sp1} <-> {sp2}: {len(rbh_12)} RBH pairs")
        print(f"  {sp2} <-> {sp3}: {len(rbh_23)} RBH pairs")
        print(f"  {sp1} <-> {sp3}: {len(rbh_13)} RBH pairs")
        
        # Find 3-way MBH
        core_families = []
        
        for prot1, prot2 in rbh_12.items():
            # Check if prot2 has RBH with sp3
            prot3 = rbh_23.get(prot2)
            if prot3:
                # Check if prot1 and prot3 are RBH
                if rbh_13.get(prot1) == prot3:
                    core_families.append({
                        sp1: prot1,
                        sp2: prot2,
                        sp3: prot3
                    })
        
        print(f"  Found {len(core_families)} 3-way MBH families\n")
        
        # Step 2: Add additional species
        extended_families = core_families.copy()
        
        def get_rbh_protein(rbh_map, query_protein, sp_query, sp_target):
                if query_protein in rbh_map:
                    return rbh_map.get(query_protein)
    
                reverse_map = {v: k for k, v in rbh_map.items()}
                return reverse_map.get(query_protein)
        
        for add_sp in rest_species:
            print(f"Step 2.{rest_species.index(add_sp)+1}: Adding {add_sp}...")
            
            # Load RBH files between new species and each core species
            rbh_add_sp1 = load_rbh(add_sp, sp1, rbh_dir)
            rbh_add_sp2 = load_rbh(add_sp, sp2, rbh_dir)
            rbh_add_sp3 = load_rbh(add_sp, sp3, rbh_dir)
            
            print(f"  {add_sp} <-> {sp1}: {len(rbh_add_sp1)} RBH pairs")
            print(f"  {add_sp} <-> {sp2}: {len(rbh_add_sp2)} RBH pairs")
            print(f"  {add_sp} <-> {sp3}: {len(rbh_add_sp3)} RBH pairs")
            
            new_extended = []
            for family in extended_families:
                # Try to find protein in new species that is RBH with all 3 core proteins
                prot1 = family[sp1]
                prot2 = family[sp2]
                prot3 = family[sp3]
                
                # Get RBH proteins from new species
                prot_add_1 = get_rbh_protein(rbh_add_sp1, prot1, sp1, add_sp)
                prot_add_2 = get_rbh_protein(rbh_add_sp2, prot2, sp2, add_sp)
                prot_add_3 = get_rbh_protein(rbh_add_sp3, prot3, sp3, add_sp)
                
                # Check if all three point to the SAME protein in new species
                if prot_add_1 and prot_add_1 == prot_add_2 == prot_add_3:
                    # Add this protein to the family
                    new_family = family.copy()
                    new_family[add_sp] = prot_add_1
                    new_extended.append(new_family)
            
            extended_families = new_extended
            print(f"  Families remaining: {len(extended_families)}\n")
        
        # Convert to output format
        gene_families = []
        gene_family_details = []
        
        for family in extended_families:
            # Get chromosomes in order of all_species
            try:    
                chromosomes = tuple(
                    chrom_maps[sp][family[sp]]
                    for sp in all_species
                )
            except KeyError as e:
                print(f'Missing key: {e}')
                continue
            
            gene_families.append(chromosomes)
            gene_family_details.append({
                'proteins': family,
                'chromosomes': chromosomes
            })
        
        print(f"Final: {len(gene_families)} gene families across all {len(all_species)} species")
    
    print("\n" + "="*100)
    print("SUMMARY")
    print("="*100)
    print(f"Total gene families: {len(gene_families)}")
    print(f"Species order: {', '.join(all_species)}")
    print("\nFirst 10 families (chromosomes):")
    for i, family in enumerate(gene_families[:10], 1):
        print(f"  Family {i}: {family}")
    print("="*100 + "\n")
    
    return gene_families, gene_family_details