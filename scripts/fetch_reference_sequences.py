#!/usr/bin/env python3
"""
Fetch COI reference sequences from GenBank based on Hoque et al. 2022
Focus on Southern California mosquito species

Citation:
Hoque MM, Valentine MJ, Kelly PJ, Barua S, Barrantes Murillo DF, Wang C. (2022)
Modification of the Folmer primers for the cytochrome c oxidase gene facilitates
identification of mosquitoes. Parasites & Vectors 15:437.
https://doi.org/10.1186/s13071-022-05494-2
"""

import os
import sys
from Bio import Entrez, SeqIO
from time import sleep
from collections import defaultdict

# Set your email for NCBI Entrez
Entrez.email = "your.email@example.com"  # CHANGE THIS!

# GenBank accessions extracted from Hoque et al. 2022 Tables 2 and 3
# Organized by species with focus on Southern California genera
GENBANK_ACCESSIONS = {
    # Aedes species (Southern California relevant)
    'Aedes_aegypti': [
        'MN299002.1', 'MK300221.1', 'MN298992.1', 'MN298993.1',
        'MK300224.1'
    ],
    'Aedes_albopictus': ['KF211505', 'MN513368'],
    'Aedes_taeniorhynchus': ['MN626442.1'],
    'Aedes_busckii': ['MN626443.1'],
    'Aedes_japonicus': ['KF211494.1'],
    'Aedes_tortilis': ['JX259682.1'],
    'Aedes_vexans': ['MH032639', 'KP954638.1'],
    'Aedes_triseriatus': ['MG242523.1'],

    # Anopheles species (Southern California relevant)
    'Anopheles_quadrimaculatus': ['L04272', 'L04272.1'],
    'Anopheles_punctipennis': ['KR653634.100', 'KR666470.1'],
    'Anopheles_crucians': ['MT040812.1'],
    'Anopheles_albimanus': ['KC354824'],
    'Anopheles_funestus': ['MT917175'],
    'Anopheles_gambiae': ['MG753769'],
    'Anopheles_pseudopunctipennis': ['KC354820'],
    'Anopheles_stephensi': ['KT899888'],

    # Culex species (Southern California relevant)
    'Culex_pipiens': [
        'KP293422.1', 'MK714012.1', 'KP293425.1',
        'MK714001.1'
    ],
    'Culex_quinquefasciatus': [
        'MW509603', 'MH463059.1', 'MN389462.1', 'MN005046.1'
    ],
    'Culex_tarsalis': ['NC_036006.1', 'AF425847'],
    'Culex_erraticus': ['MH129001', 'MH128999.1', 'MN389459'],
    'Culex_nigripalpus': ['NC_037823.1'],
    'Culex_usquatissimus': ['NC_036007.1'],
    'Culex_coronator': ['NC_036006.1'],
    'Culex_sitiens': ['NC_054318'],
    'Culex_tritaeniorhynchus': ['KT852976'],

    # Other genera for reference/comparison
    'Deinocerites_magnus': ['MH376751.1'],
    'Psorophora_howardii': ['MG242538.1'],
    'Psorophora_ferox': ['KY782650', 'MK575485'],
    'Psorophora_pygmaea': ['JX260116.1'],
    'Psorophora_cingulata': ['KM592989.1'],
    'Psorophora_confinis': ['KY859921.1'],
    'Uranotaenia_sapphirina': ['GU908127.1'],
    'Culiseta_annulata': ['MN626442'],
    'Mansonia_annulata': ['HQ341635'],
    'Ochlerotatus_taeniorhynchus': ['MN626442'],
}

# Southern California priority species (will be fetched first)
SOCAL_PRIORITY = [
    'Aedes_aegypti',
    'Aedes_albopictus',
    'Aedes_taeniorhynchus',
    'Culex_pipiens',
    'Culex_quinquefasciatus',
    'Culex_tarsalis',
    'Culex_erraticus',
    'Anopheles_quadrimaculatus',
    'Anopheles_punctipennis',
]


def fetch_sequence(accession, retries=3):
    """Fetch a single sequence from GenBank with retry logic"""
    for attempt in range(retries):
        try:
            print(f"  Fetching {accession} (attempt {attempt + 1}/{retries})...")
            handle = Entrez.efetch(db="nucleotide", id=accession,
                                   rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            sleep(0.34)  # NCBI rate limit: max 3 requests/second
            return record
        except Exception as e:
            print(f"  Error fetching {accession}: {e}")
            if attempt < retries - 1:
                sleep(2)
            else:
                return None
    return None


def extract_coi_region(record):
    """Extract COI gene region if annotated, otherwise return full sequence"""
    # Look for COI gene annotation
    for feature in record.features:
        if feature.type == "CDS" or feature.type == "gene":
            if 'gene' in feature.qualifiers:
                gene_name = feature.qualifiers['gene'][0].upper()
                if 'COI' in gene_name or 'COX1' in gene_name or 'COXI' in gene_name:
                    return feature.extract(record.seq)

    # If no COI annotation found, return the whole sequence
    # (most records in the paper are already COI fragments)
    return record.seq


def create_fasta_header(species_name, accession, record):
    """Create informative FASTA header"""
    # Format: >Species_name_Accession [organism info] [location if available]
    organism = record.annotations.get('organism', species_name.replace('_', ' '))

    # Try to get country/location info
    location = ""
    for feature in record.features:
        if feature.type == "source":
            if 'country' in feature.qualifiers:
                location = f" [{feature.qualifiers['country'][0]}]"
            break

    header = f">{species_name}_{accession} {organism}{location}"
    return header


def main():
    """Main function to fetch and save sequences"""
    output_file = "/Users/lucianocosme/Projects/dna-barcoding-analysis/data/reference/socal_mosquitoes.fasta"

    print("=" * 80)
    print("Fetching COI reference sequences from Hoque et al. 2022")
    print("=" * 80)
    print(f"\nOutput file: {output_file}")
    print(f"Total species: {len(GENBANK_ACCESSIONS)}")
    print(f"Total accessions to fetch: {sum(len(accs) for accs in GENBANK_ACCESSIONS.values())}\n")

    sequences = []
    stats = defaultdict(int)
    failed = []

    # Process priority species first
    print("\n" + "=" * 80)
    print("PRIORITY: Southern California Species")
    print("=" * 80)

    for species in SOCAL_PRIORITY:
        if species not in GENBANK_ACCESSIONS:
            continue

        print(f"\n{species.replace('_', ' ')}:")
        accessions = GENBANK_ACCESSIONS[species]

        for accession in accessions:
            record = fetch_sequence(accession)

            if record:
                coi_seq = extract_coi_region(record)
                header = create_fasta_header(species, accession, record)
                sequences.append((header, str(coi_seq)))
                stats['success'] += 1
                stats[species] += 1
                print(f"    ✓ {accession}: {len(coi_seq)} bp")
            else:
                failed.append((species, accession))
                stats['failed'] += 1
                print(f"    ✗ {accession}: FAILED")

    # Process remaining species
    print("\n" + "=" * 80)
    print("Additional Reference Species")
    print("=" * 80)

    for species, accessions in GENBANK_ACCESSIONS.items():
        if species in SOCAL_PRIORITY:
            continue  # Already processed

        print(f"\n{species.replace('_', ' ')}:")

        for accession in accessions:
            record = fetch_sequence(accession)

            if record:
                coi_seq = extract_coi_region(record)
                header = create_fasta_header(species, accession, record)
                sequences.append((header, str(coi_seq)))
                stats['success'] += 1
                stats[species] += 1
                print(f"    ✓ {accession}: {len(coi_seq)} bp")
            else:
                failed.append((species, accession))
                stats['failed'] += 1
                print(f"    ✗ {accession}: FAILED")

    # Write sequences to file
    print(f"\n\nWriting {len(sequences)} sequences to {output_file}...")

    with open(output_file, 'w') as f:
        for header, seq in sequences:
            f.write(f"{header}\n")
            # Write sequence in 80-character lines
            for i in range(0, len(seq), 80):
                f.write(f"{seq[i:i+80]}\n")

    # Print summary statistics
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"\nTotal sequences fetched: {stats['success']}")
    print(f"Failed accessions: {stats['failed']}")
    print(f"Success rate: {stats['success']/(stats['success']+stats['failed'])*100:.1f}%")

    print("\n\nSequences by species:")
    for species in sorted(GENBANK_ACCESSIONS.keys()):
        if stats[species] > 0:
            priority = " [SoCal Priority]" if species in SOCAL_PRIORITY else ""
            print(f"  {species.replace('_', ' ')}: {stats[species]} sequences{priority}")

    if failed:
        print(f"\n\nFailed accessions ({len(failed)}):")
        for species, accession in failed:
            print(f"  {species}: {accession}")

    print(f"\n✓ Complete! Sequences saved to:\n  {output_file}\n")


if __name__ == "__main__":
    # Check if BioPython is installed
    try:
        import Bio
        print(f"Using BioPython version {Bio.__version__}\n")
    except ImportError:
        print("ERROR: BioPython is not installed!")
        print("Please install it with: pip install biopython")
        sys.exit(1)

    main()
