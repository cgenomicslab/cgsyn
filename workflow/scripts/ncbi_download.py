"""
NCBI Genome Downloader for Synteny Pipeline
Intelligent batch download with assembly quality fallbacks
"""

import subprocess
import json
import sys
import os
import gzip
import shutil
from pathlib import Path
import re

def check_datasets_installed():
    """Check if NCBI datasets is installed."""
    try:
        subprocess.run(['datasets', '--version'], 
                      capture_output=True, check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False

def is_taxid(query):
    """Check if query is a taxonomy ID (pure digits)."""
    return query.strip().isdigit()

def search_genome(species_query, assembly_level='complete', 
                  refseq_only=True, verbose=True):
    """
    Search for optimal genome assembly with smart fallbacks.
    
    Parameters:
    -----------
    species_query : str
        Species name or taxonomy ID
    assembly_level : str
        'complete', 'chromosome', 'scaffold', or 'contig'
    refseq_only : bool
        Only search RefSeq
    verbose : bool
        Print search details
    
    Returns:
    --------
    dict : Assembly metadata or None
    """
    
    if verbose:
        source_str = "RefSeq" if refseq_only else "RefSeq+GenBank"
        print(f"   Searching {source_str} for {assembly_level}-level assemblies...")
    
    # Build search command
    cmd = [
        'datasets', 'summary', 'genome', 'taxon', species_query,
        '--assembly-level', assembly_level,
        '--exclude-atypical',
        '--annotated'
    ]
    
    if refseq_only:
        cmd.extend(['--assembly-source', 'RefSeq'])
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        data = json.loads(result.stdout)
        
        if 'reports' not in data or len(data['reports']) == 0:
            return None
        
        assemblies = data['reports']
        
        # Sort by priority
        def assembly_priority(asm):
            # Reference genome = highest priority
            ref_priority = 0 if asm.get('assembly_info', {}).get('refseq_category') == 'reference genome' else 10
            
            # Complete > Chromosome > Scaffold > Contig
            level_priority = {
                'Complete Genome': 0,
                'Chromosome': 1,
                'Scaffold': 2,
                'Contig': 3
            }.get(asm.get('assembly_info', {}).get('assembly_level', ''), 999)
            
            # Latest assembly date (parse as negative for sorting)
            try:
                date_str = asm.get('assembly_info', {}).get('release_date', '1900-01-01')
                date_priority = -int(date_str.replace('-', ''))
            except:
                date_priority = 0
            
            return (ref_priority, level_priority, date_priority)
        
        # Get best assembly
        best = sorted(assemblies, key=assembly_priority)[0]
        
        return best
        
    except subprocess.CalledProcessError:
        return None

def get_assembly_info_string(assembly):
    """Format assembly info for display."""
    accession = assembly['accession']
    name = assembly['assembly_info']['assembly_name']
    level = assembly['assembly_info']['assembly_level']
    organism = assembly['organism']['organism_name']
    
    # Get scaffold/contig count if applicable
    stats = assembly.get('assembly_stats', {})
    scaffold_count = stats.get('number_of_scaffolds', 0)
    contig_count = stats.get('number_of_contigs', 0)
    
    info = f"{organism}\n"
    info += f"      Assembly: {name} ({accession})\n"
    info += f"      Level: {level}"
    
    if level == 'Scaffold' and scaffold_count > 0:
        info += f" ({scaffold_count:,} scaffolds)"
    elif level == 'Contig' and contig_count > 0:
        info += f" ({contig_count:,} contigs)"
    
    return info

def find_optimal_assembly(species_query, interactive=True):
    """
    Find optimal assembly with intelligent fallbacks.
    
    Search order:
    1. Complete genome in RefSeq
    2. Chromosome-level in RefSeq
    3. Complete genome in GenBank (with user confirmation)
    4. Chromosome-level in GenBank (with user confirmation)
    5. Scaffold-level in RefSeq (with user confirmation + warning)
    6. Scaffold-level in GenBank (with user confirmation + warning)
    7. Contig-level in RefSeq (with user confirmation + warning)
    8. Contig-level in GenBank (with user confirmation + warning)
    
    Returns:
    --------
    dict : Best assembly or None
    """
    
    query_type = "TaxID" if is_taxid(species_query) else "species name"
    print(f"\n🔍 Searching NCBI for {query_type}: {species_query}")
    
    # Try 1: Complete genome in RefSeq
    assembly = search_genome(species_query, 'complete', refseq_only=True)
    if assembly:
        print(f"✅ Found complete RefSeq genome:")
        print(f"   {get_assembly_info_string(assembly)}")
        return assembly
    
    # Try 2: Chromosome-level in RefSeq
    assembly = search_genome(species_query, 'chromosome', refseq_only=True)
    if assembly:
        print(f"✅ Found chromosome-level RefSeq genome:")
        print(f"   {get_assembly_info_string(assembly)}")
        return assembly
    
    # Try 3: Complete genome in GenBank
    print("   ⚠️  No RefSeq complete/chromosome assemblies found")
    assembly = search_genome(species_query, 'complete', refseq_only=False)
    if assembly:
        print(f"⚠️  Found complete GenBank genome (not RefSeq-curated):")
        print(f"   {get_assembly_info_string(assembly)}")
        
        if interactive:
            response = input("   Use this assembly? [y/N]: ").strip().lower()
            if response == 'y':
                return assembly
        else:
            print("   Using GenBank assembly (non-interactive mode)")
            return assembly
    
    # Try 4: Chromosome-level in GenBank
    assembly = search_genome(species_query, 'chromosome', refseq_only=False)
    if assembly:
        print(f"⚠️  Found chromosome-level GenBank genome (not RefSeq-curated):")
        print(f"   {get_assembly_info_string(assembly)}")
        
        if interactive:
            response = input("   Use this assembly? [y/N]: ").strip().lower()
            if response == 'y':
                return assembly
        else:
            print("   Using GenBank assembly (non-interactive mode)")
            return assembly
    
    # Try 5: Scaffold-level in RefSeq
    assembly = search_genome(species_query, 'scaffold', refseq_only=True, verbose=False)
    if assembly:
        print(f"⚠️  Only scaffold-level RefSeq assembly found:")
        print(f"   {get_assembly_info_string(assembly)}")
        
        stats = assembly.get('assembly_stats', {})
        scaffold_count = stats.get('number_of_scaffolds', 0)
        
        print(f"\n   ⚠️  WARNING: Scaffold-level assemblies may not be suitable for synteny analysis!")
        print(f"   This assembly has {scaffold_count:,} scaffolds (not chromosomes)")
        
        if interactive:
            response = input("   Proceed with scaffold-level assembly? [y/N]: ").strip().lower()
            if response == 'y':
                return assembly
        else:
            print("   Skipping scaffold-level assembly (non-interactive mode)")
            return None
    
    # Try 6: Scaffold-level in GenBank
    assembly = search_genome(species_query, 'scaffold', refseq_only=False, verbose=False)
    if assembly:
        print(f"⚠️  Only scaffold-level GenBank assembly found:")
        print(f"   {get_assembly_info_string(assembly)}")
        
        stats = assembly.get('assembly_stats', {})
        scaffold_count = stats.get('number_of_scaffolds', 0)
        
        print(f"\n   ⚠️  WARNING: Scaffold-level assemblies may not be suitable for synteny analysis!")
        print(f"   This assembly has {scaffold_count:,} scaffolds (not chromosomes)")
        
        if interactive:
            response = input("   Proceed with scaffold-level GenBank assembly? [y/N]: ").strip().lower()
            if response == 'y':
                return assembly
        else:
            print("   Skipping scaffold-level assembly (non-interactive mode)")
            return None
    
    # Try 7: Contig-level in RefSeq
    assembly = search_genome(species_query, 'contig', refseq_only=True, verbose=False)
    if assembly:
        print(f"⚠️  Only contig-level RefSeq assembly found:")
        print(f"   {get_assembly_info_string(assembly)}")
        
        stats = assembly.get('assembly_stats', {})
        contig_count = stats.get('number_of_contigs', 0)
        
        print(f"\n   🛑 WARNING: Contig-level assemblies are NOT suitable for synteny analysis!")
        print(f"   This assembly has {contig_count:,} contigs (highly fragmented)")
        print(f"   Synteny analysis requires chromosome or scaffold-level assemblies.")
        
        if interactive:
            response = input("   Proceed anyway? (NOT recommended) [y/N]: ").strip().lower()
            if response == 'y':
                return assembly
        else:
            print("   Skipping contig-level assembly (non-interactive mode)")
            return None
    
    # Try 8: Contig-level in GenBank
    assembly = search_genome(species_query, 'contig', refseq_only=False, verbose=False)
    if assembly:
        print(f"⚠️  Only contig-level GenBank assembly found:")
        print(f"   {get_assembly_info_string(assembly)}")
        
        stats = assembly.get('assembly_stats', {})
        contig_count = stats.get('number_of_contigs', 0)
        
        print(f"\n   🛑 WARNING: Contig-level assemblies are NOT suitable for synteny analysis!")
        print(f"   This assembly has {contig_count:,} contigs (highly fragmented)")
        print(f"   Synteny analysis requires chromosome or scaffold-level assemblies.")
        
        if interactive:
            response = input("   Proceed anyway? (NOT recommended) [y/N]: ").strip().lower()
            if response == 'y':
                return assembly
        else:
            print("   Skipping contig-level assembly (non-interactive mode)")
            return None
    
    # Nothing found at all
    print(f"❌ No assemblies with functional annotation found for: {species_query}")
    print(f"   This species may have genome assemblies available in NCBI, but without functional annotation files (.gff), which are necessary for synteny analysis.")
    return None

def download_genome(accession, output_dir, species_label):
    """
    Download genome files using NCBI datasets.
    
    Parameters:
    -----------
    accession : str
        Assembly accession
    output_dir : str
        Output directory
    species_label : str
        User-friendly label
    
    Returns:
    --------
    bool : Success status
    """
    
    print(f"   ⬇️  Downloading {species_label}...")
    
    # Create temp directory
    temp_dir = Path(output_dir) / f'temp_download_{species_label}'
    temp_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # Download
        cmd = [
            'datasets', 'download', 'genome', 'accession', accession,
            '--include', 'gff3,protein',
            '--filename', str(temp_dir / 'genome.zip')
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"   ❌ Download command failed: {result.stderr}")
            return False
        
        # Unzip
        subprocess.run(['unzip', '-q', str(temp_dir / 'genome.zip'), 
                       '-d', str(temp_dir)], 
                      capture_output=True, check=True)
        
        # Find files
        data_dir = temp_dir / 'ncbi_dataset' / 'data' / accession
        
        gff_files = list(data_dir.glob('genomic.gff')) or list(data_dir.glob('*.gff'))
        protein_files = list(data_dir.glob('protein.faa'))
        
        if not gff_files or not protein_files:
            print(f"   ❌ Could not find GFF or protein files in download")
            return False
        
        # Create output directories
        gff_out_dir = Path(output_dir) / 'gff'
        prot_out_dir = Path(output_dir) / 'proteomes'
        gff_out_dir.mkdir(parents=True, exist_ok=True)
        prot_out_dir.mkdir(parents=True, exist_ok=True)
        
        # Copy and compress
        gff_output = gff_out_dir / f'{species_label}.gff.gz'
        prot_output = prot_out_dir / f'{species_label}.faa.gz'
        
        with open(gff_files[0], 'rb') as f_in:
            with gzip.open(gff_output, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        with open(protein_files[0], 'rb') as f_in:
            with gzip.open(prot_output, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        print(f"   ✅ Saved to resources/gff/{species_label}.gff.gz")
        print(f"   ✅ Saved to resources/proteomes/{species_label}.faa.gz")
        
        # Cleanup
        shutil.rmtree(temp_dir)
        return True
        
    except Exception as e:
        print(f"   ❌ Download failed: {e}")
        shutil.rmtree(temp_dir, ignore_errors=True)
        return False

def main():
    """Main batch download function."""
    
    if not check_datasets_installed():
        print("❌ Error: NCBI datasets tool not found!")
        print("\nInstall with:")
        print("  conda install -c conda-forge ncbi-datasets-cli")
        print("  OR download from: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/")
        sys.exit(1)
    
    import argparse
    parser = argparse.ArgumentParser(
        description='Download genomes from NCBI',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single species
  %(prog)s --species "Homo sapiens" --labels Hsapiens
  
  # Multiple species (same number of labels required)
  %(prog)s --species "Homo sapiens,Mus musculus,7227" --labels Hsapiens,Mmusculus,Dmelanogaster
  
  # Mix taxids and names
  %(prog)s --species "9606,Mus musculus,Drosophila" --labels human,mouse,fly
        """
    )
    
    parser.add_argument('--species', required=True,
                       help='Comma-separated list of species names or taxids')
    parser.add_argument('--labels', required=True,
                       help='Comma-separated list of output labels (must match species count)')
    parser.add_argument('--output-dir', default='resources',
                       help='Output directory (default: resources)')
    parser.add_argument('--non-interactive', action='store_true',
                       help='Skip user prompts (auto-accept RefSeq/GenBank, skip scaffolds)')
    
    args = parser.parse_args()
    
    # Parse lists
    species_list = [s.strip() for s in args.species.split(',')]
    label_list = [l.strip() for l in args.labels.split(',')]
    
    if len(species_list) != len(label_list):
        print(f"❌ Error: Number of species ({len(species_list)}) must match number of labels ({len(label_list)})")
        sys.exit(1)
    
    print(f"\n{'='*70}")
    print(f"NCBI GENOME BATCH DOWNLOAD")
    print(f"{'='*70}")
    print(f"Species to download: {len(species_list)}")
    
    results = []
    
    for species_query, label in zip(species_list, label_list):
        print(f"\n{'─'*70}")
        print(f"Processing: {species_query} → {label}")
        print(f"{'─'*70}")
        
        # Find optimal assembly
        assembly = find_optimal_assembly(
            species_query, 
            interactive=not args.non_interactive
        )
        
        if not assembly:
            print(f"⏭️  Skipping {species_query}")
            results.append((species_query, label, False))
            continue
        
        # Download
        accession = assembly['accession']
        success = download_genome(accession, args.output_dir, label)
        results.append((species_query, label, success))
    
    # Summary
    print(f"\n{'='*70}")
    print(f"DOWNLOAD SUMMARY")
    print(f"{'='*70}")
    
    successful = sum(1 for _, _, success in results if success)
    failed = len(results) - successful
    
    for species_query, label, success in results:
        status = "✅" if success else "❌"
        print(f"{status} {species_query:30s} → {label}")
    
    print(f"\nTotal: {successful}/{len(results)} successful")
    
    if failed > 0:
        print(f"\n⚠️  {failed} download(s) failed")
        sys.exit(1)
    else:
        print(f"\n🎉 All genomes downloaded successfully!")
        sys.exit(0)

if __name__ == '__main__':
    main()
