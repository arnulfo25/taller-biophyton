#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Solución del Ejercicio de Filogenia (Versión Robusta con BLAST)

Objetivo:
Construir un árbol filogenético a partir de secuencias de genes extraídas
mediante la localización de primers con BLASTn, un método más preciso.

Pasos:
1.  Para cada genoma, crear una base de datos BLAST temporal.
2.  Usar BLASTn para encontrar la ubicación de los primers forward y reverse.
3.  Extraer la secuencia del amplicón delimitada por los primers.
4.  Guardar los amplicones en un archivo FASTA.
5.  Alinear los amplicones con MUSCLE y construir el árbol.
"""

import os
import subprocess
import tempfile
from Bio.Seq import Seq
from Bio import SeqIO, AlignIO, Phylo
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

# --- Configuración ---
GENOMES_DIR = "genomes"
GENE_FASTA_OUTPUT = "gene_sequences_from_blast.fasta"
ALIGNED_GENES_OUTPUT = "aligned_genes_from_blast.aln"

ALL_PRIMERS = {
    "PED-N_Forward": "CGCAAAGACTGAACCCACTAAC",
    "PED-N_Reverse": "TTGCCTCTGTTGTTACTTGGAG",
    "TTSuV1F": "GTAGATGCCACAGTCGATGAAC",
    "TTSuV1R": "CCTACGTCAAGTGGAACACAGG",
    "TTSuV2F": "ATGGCCAGCCAGACAATAAAGA",
    "TTSuV2R": "GCAACACCCAGTAGGTTGACTT",
}

ASSAYS_FOR_EXTRACTION = {
    "PEDv_N_Gene": {
        "forward_primer": "PED-N_Forward",
        "reverse_primer": "PED-N_Reverse",
        "genome_file": "PEDv_genome_new.fasta",
    },
    "TTSuV1_Region": {
        "forward_primer": "TTSuV1F",
        "reverse_primer": "TTSuV1R",
        "genome_file": "TTSuV1_genome.fasta",
    },
    "TTSuV2_Region": {
        "forward_primer": "TTSuV2F",
        "reverse_primer": "TTSuV2R",
        "genome_file": "TTSuV2_genome.fasta",
    },
}

def run_command(cmd):
    """Ejecuta un comando y devuelve su salida."""
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        # print(f"Error ejecutando `{' '.join(cmd)}`: {e.stderr}")
        return None

def get_primer_hit_coords(primer_seq: str, db_path: str) -> tuple:
    """Usa BLASTn para encontrar la mejor posición de un primer."""
    query_fasta = f">primer\n{primer_seq}"
    cmd = [
        "blastn", "-task", "blastn-short", "-query", "-", "-db", db_path,
        "-outfmt", "6 sstart send", "-max_target_seqs", "1"
    ]
    process = subprocess.run(cmd, input=query_fasta, capture_output=True, text=True)
    if process.returncode == 0 and process.stdout:
        parts = process.stdout.strip().split()
        return int(parts[0]), int(parts[1])
    return None, None

def extract_and_save_gene_sequences():
    print(f"--- PASO 1/2: Extrayendo secuencias de genes usando BLAST y guardando en '{GENE_FASTA_OUTPUT}' ---")
    gene_sequences = []

    with tempfile.TemporaryDirectory() as tempdir:
        for gene_name, assay in ASSAYS_FOR_EXTRACTION.items():
            genome_file = os.path.join(GENOMES_DIR, assay["genome_file"])
            db_path = os.path.join(tempdir, assay["genome_file"])
            
            run_command(["makeblastdb", "-in", genome_file, "-dbtype", "nucl", "-out", db_path])

            fw_primer = ALL_PRIMERS[assay["forward_primer"]]
            rv_primer = ALL_PRIMERS[assay["reverse_primer"]]

            f_start, f_end = get_primer_hit_coords(fw_primer, db_path)
            r_start, r_end = get_primer_hit_coords(reverse_complement(rv_primer), db_path)

            if f_start and r_start:
                start = min(f_start, r_start)
                end = max(f_end, r_end)
                
                genome_record = SeqIO.read(genome_file, "fasta")
                amplicon_seq = genome_record.seq[start-1:end]
                
                record = SeqRecord(amplicon_seq, id=gene_name, description=f"Amplicón de {assay['genome_file']}")
                gene_sequences.append(record)
                print(f"  - Extraído: {gene_name} ({len(amplicon_seq)} pb) de {assay['genome_file']}")
            else:
                print(f"  - ADVERTENCIA: No se pudo encontrar el amplicón para '{gene_name}' en {assay['genome_file']}")

    SeqIO.write(gene_sequences, GENE_FASTA_OUTPUT, "fasta")
    print("¡Secuencias de genes guardadas con éxito!\n")

def run_muscle_alignment(input_fasta: str, output_alignment: str):
    print(f"--- PASO 3: Alineando secuencias de genes con MUSCLE ---")
    muscle_cmd = ["muscle", "-align", input_fasta, "-output", output_alignment]
    run_command(muscle_cmd)
    print(f"Alineamiento guardado en '{output_alignment}'\n")

def build_and_visualize_tree(aligned_fasta: str):
    print(f"--- PASO 4: Construyendo y visualizando el árbol ---")
    try:
        if os.path.getsize(aligned_fasta) > 0:
            aln = AlignIO.read(aligned_fasta, "fasta")
            calculator = DistanceCalculator('identity')
            dist_matrix = calculator.get_distance(aln)
            constructor = DistanceTreeConstructor()
            tree = constructor.nj(dist_matrix)
            if tree:
                tree.root_at_midpoint()
                Phylo.draw_ascii(tree)
    except (ValueError, FileNotFoundError):
        print("No se pudo construir el árbol. El archivo de alineamiento podría estar vacío.")

def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())

def main():
    print("### Inicio del Ejercicio de Filogenia (versión BLAST) ###\n")
    extract_and_save_gene_sequences()
    run_muscle_alignment(GENE_FASTA_OUTPUT, ALIGNED_GENES_OUTPUT)
    build_and_visualize_tree(ALIGNED_GENES_OUTPUT)
    print("\n### Fin del Ejercicio ###")

if __name__ == "__main__":
    main()
