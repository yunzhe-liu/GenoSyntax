"""Generate synthetic Panaroo/Pharokka-shaped inputs using the standard library."""

import argparse
import csv
from pathlib import Path


def make_demo(destination):
    destination = Path(destination)
    fasta_dir = destination / 'fastas'
    annotation_dir = destination / 'pharokka'
    fasta_dir.mkdir(parents=True, exist_ok=True)
    annotation_dir.mkdir(parents=True, exist_ok=True)
    genomes = [f'demo_{i:02d}' for i in range(12)]
    families = ['core_a', 'core_b', 'accessory_x', 'accessory_y',
                'accessory_z', 'core_c', 'core_d']
    metadata = ['Gene', 'Non-unique Gene name', 'Annotation', 'No. isolates',
                'No. sequences', 'Avg sequences per isolate', 'Genome Fragment',
                'Order within Fragment', 'Accessory Fragment',
                'Accessory Order with Fragment', 'QC', 'Min group size nuc',
                'Max group size nuc', 'Avg group size nuc']
    loci = {}
    for i, genome in enumerate(genomes):
        genome_dir = annotation_dir / genome
        genome_dir.mkdir(exist_ok=True)
        selected = families if i < 6 else [f for f in families if f.startswith('core_')]
        sequence = 'ATG' * 250
        fasta = f'>{genome}\n{sequence}\n'
        (fasta_dir / f'{genome}.fasta').write_text(fasta)
        (genome_dir / f'{genome}.fna').write_text(fasta)
        gff, faa = ['##gff-version 3'], []
        for j, family in enumerate(selected):
            locus = f'{genome}_{family}'
            loci[genome, family] = locus
            start = 51 + j * 70
            gff.append(f'{genome}\tdemo\tCDS\t{start}\t{start + 29}\t.\t+\t0\tID={locus}')
            faa.append(f'>{locus}\nMMMMMMMMMM')
        (genome_dir / f'{genome}.gff').write_text('\n'.join(gff) + '\n')
        (genome_dir / 'phanotate.faa').write_text('\n'.join(faa) + '\n')
    with (destination / 'gene_presence_absence_roary.csv').open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(metadata + genomes)
        for family in families:
            count = sum((g, family) in loci for g in genomes)
            writer.writerow([family, '', 'hypothetical protein', count, count, 1,
                             '', '', '', '', '', 30, 30, 30] +
                            [loci.get((g, family), '') for g in genomes])
    print(f'Synthetic inputs written to {destination}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('destination', type=Path)
    make_demo(parser.parse_args().destination)
