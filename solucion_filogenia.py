#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Solución del Ejercicio de Filogenia: Construcción de un Árbol Filogenético

Este script es la solución completa al ejercicio `ejercicio_filogenia.py`.
"""

import os
import subprocess
from Bio import SeqIO, AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt

# ==============================================================================
# --- CONFIGURACIÓN ---
# ==============================================================================

GENOMES_DIR = "genomes"
GENOME_FILES = [
    "PEDv_genome.fasta",
    "TTSuV1_genome.fasta",
    "TTSuV2_genome.fasta",
]
COMBINED_FASTA = "combined_genomes.fasta"
ALIGNED_FASTA = "aligned_genomes.aln"
SVG_OUTPUT_FILE = "phylogenetic_tree.svg"
SVG_OUTPUT_FILE = "phylogenetic_tree.svg"
SVG_OUTPUT_FILE = "phylogenetic_tree.svg"
SVG_OUTPUT_FILE = "phylogenetic_tree.svg"
SVG_OUTPUT_FILE = "phylogenetic_tree.svg"

# ==============================================================================
# --- PASO 1: Combinar Archivos FASTA ---
# ==============================================================================

def combine_fasta_files(file_list: list, output_file: str):
    """
    Combina múltiples archivos FASTA en uno solo.
    """
    print(f"--- PASO 1: Combinando {len(file_list)} archivos FASTA en '{output_file}' ---")
    
    # ===== SOLUCIÓN TAREA 1 =====
    all_sequences = []
    for file_path in file_list:
        # Usamos extend para añadir todas las secuencias de un archivo a la lista
        all_sequences.extend(list(SeqIO.parse(file_path, "fasta")))

    # Escribimos la lista de secuencias al archivo de salida
    SeqIO.write(all_sequences, output_file, "fasta")
    # ===== FIN SOLUCIÓN TAREA 1 =====

    print("¡Archivos combinados con éxito!\n")

# ==============================================================================
# --- PASO 2: Realizar Alineamiento Múltiple con MUSCLE ---
# ==============================================================================

def run_muscle_alignment(input_fasta: str, output_alignment: str):
    """
    Ejecuta MUSCLE para realizar un alineamiento múltiple de secuencias.
    """
    print(f"--- PASO 2: Alineando secuencias de '{input_fasta}' con MUSCLE ---")
    # Usamos el algoritmo -super5 que es más rápido para secuencias largas
    muscle_cmd = ["muscle", "-super5", input_fasta, "-output", output_alignment]
    
    # ===== SOLUCIÓN TAREA 2 =====
    try:
        # Ejecutamos el comando, verificando que termine correctamente (check=True)
        subprocess.run(muscle_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(f"Alineamiento guardado en '{output_alignment}'\n")
    except FileNotFoundError:
        print("Error: El comando 'muscle' no se encontró.")
        print("Por favor, asegúrate de que MUSCLE esté instalado y en el PATH del sistema.")
        exit()
    except subprocess.CalledProcessError as e:
        print(f"Ocurrió un error al ejecutar MUSCLE: {e}")
        print(f"Stderr: {e.stderr.decode()}")
        exit()
    # ===== FIN SOLUCIÓN TAREA 2 =====

# ==============================================================================
# --- PASO 3: Construir el Árbol Filogenético ---
# ==============================================================================

def build_phylogenetic_tree(aligned_fasta: str):
    """
    Construye y devuelve un árbol filogenético a partir de un alineamiento.
    """
    print(f"--- PASO 3: Construyendo árbol desde '{aligned_fasta}' ---")
    
    # ===== SOLUCIÓN TAREA 3 =====
    # 1. Lee el alineamiento en formato fasta, que es el default de muscle v5
    aln = AlignIO.read(aligned_fasta, "fasta")

    # 2. Calcula la matriz de distancias usando el modelo de identidad
    calculator = DistanceCalculator('identity')
    dist_matrix = calculator.get_distance(aln)

    # 3. Construye el árbol usando el algoritmo Neighbor-Joining
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dist_matrix)
    # ===== FIN SOLUCIÓN TAREA 3 =====

    print("¡Árbol construido con éxito!\n")
    return tree

# ==============================================================================
# --- PASO 4: Visualizar el Árbol ---
# ==============================================================================

def visualize_tree(tree):
    """
    Visualiza un árbol filogenético en la consola en formato ASCII.
    """
    print("--- PASO 4: Visualización del Árbol Filogenético ---")
    
    # ===== SOLUCIÓN TAREA 4 =====
    # Hacemos el árbol un poco más legible estableciendo el root
    tree.root_at_midpoint()
    Phylo.draw_ascii(tree)
    # ===== FIN SOLUCIÓN TAREA 4 =====

# ==============================================================================
# --- BLOQUE PRINCIPAL DE EJECUCIÓN ---
# ==============================================================================

def main():
    """Función principal que orquesta todos los pasos del ejercicio."""
    print("### Inicio del Ejercicio de Análisis Filogenético ###\n")

    genome_paths = [os.path.join(GENOMES_DIR, fname) for fname in GENOME_FILES]

    combine_fasta_files(genome_paths, COMBINED_FASTA)
    run_muscle_alignment(COMBINED_FASTA, ALIGNED_FASTA)
    phylo_tree = build_phylogenetic_tree(ALIGNED_FASTA)

    if phylo_tree:
        visualize_tree(phylo_tree)

    print("\n### Fin del Ejercicio ###")

if __name__ == "__main__":
    main()