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
import matplotlib
matplotlib.use('Agg')  # Usar el backend Agg para no requerir GUI
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
# --- CONFIGURACIÓN ---
# ==============================================================================

# Usaremos un solo archivo de entrada con secuencias de genes, que es más pequeño
INPUT_SEQUENCES = "final_gene_sequences.fasta" # Usamos el archivo que sí tiene contenido
ALIGNED_FASTA = "aligned_genes.aln" # Nuevo archivo de salida para el alineamiento
OUTPUT_IMAGE = "arbol_genes_filogenetico.png" # Nuevo archivo de imagen


# ==============================================================================
# --- PASO 1: Combinar Archivos FASTA (NO NECESARIO) ---
# ==============================================================================

# Ya no necesitamos esta función porque partimos de un solo archivo.

# ==============================================================================
# --- PASO 2: Realizar Alineamiento Múltiple con MUSCLE ---
# ==============================================================================

def run_muscle_alignment(input_fasta: str, output_alignment: str):
    """
    Ejecuta MUSCLE para realizar un alineamiento múltiple de secuencias.
    """
    print(f"--- PASO 2: Alineando secuencias de '{input_fasta}' con MUSCLE ---")
    # Usamos el algoritmo -super5 que es rápido y adecuado para genes
    muscle_cmd = ["muscle", "-super5", input_fasta, "-output", output_alignment]
    
    try:
        # Ejecutamos el comando, pero no verificamos el código de salida (check=False)
        # Ejecutamos el comando y capturamos su salida para depuración.
        result = subprocess.run(muscle_cmd, check=False, capture_output=True, text=True)
        print("SALIDA DE MUSCLE (stdout):")
        print(result.stdout)
        print("ERRORES DE MUSCLE (stderr):")
        print(result.stderr)
        if result.returncode != 0:
            print(f"ADVERTENCIA: MUSCLE terminó con código de salida {result.returncode}")
        print(f"Alineamiento guardado en '{output_alignment}'\n")
    except FileNotFoundError:
        print("Error: El comando 'muscle' no se encontró.")
        print("Por favor, asegúrate de que MUSCLE esté instalado y en el PATH del sistema.")
        exit()

# ==============================================================================
# --- PASO 3: Construir el Árbol Filogenético ---
# ==============================================================================

def build_phylogenetic_tree(aligned_fasta: str):
    """
    Construye y devuelve un árbol filogenético a partir de un alineamiento.
    """
    print(f"--- PASO 3: Construyendo árbol desde '{aligned_fasta}' ---")
    
    try:
        aln = AlignIO.read(aligned_fasta, "fasta")
    except ValueError:
        print(f"Error: El archivo de alineamiento '{aligned_fasta}' está vacío o dañado.")
        exit()

    # Calcula la matriz de distancias
    calculator = DistanceCalculator('identity')
    dist_matrix = calculator.get_distance(aln)

    # Construye el árbol
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dist_matrix)

    print("¡Árbol construido con éxito!\n")
    return tree

# ==============================================================================
# --- PASO 4: Guardar el Árbol como Imagen ---
# ==============================================================================

def save_tree_image(tree, output_file: str):
    """
    Guarda un árbol filogenético como una imagen PNG.
    """
    print(f"--- PASO 4: Guardando el árbol en '{output_file}' ---")
    
    tree.root_at_midpoint()
    tree.ladderize()

    fig = plt.figure(figsize=(12, 8), dpi=100)
    axes = fig.add_subplot(1, 1, 1)
    
    Phylo.draw(tree, axes=axes, do_show=False)
    
    plt.title("Árbol Filogenético de Genes Virales")
    plt.xlabel("Distancia Evolutiva")
    plt.ylabel("Secuencias")
    fig.tight_layout()
    
    plt.savefig(output_file)
    print(f"Gráfico del árbol guardado en '{output_file}'.")

# ==============================================================================
# --- BLOQUE PRINCIPAL DE EJECUCIÓN ---
# ==============================================================================

def main():
    """Función principal que orquesta todos los pasos del ejercicio."""
    print("### Inicio del Ejercicio de Análisis Filogenético con Genes ###\n")

    # Como usamos un solo archivo de entrada, no necesitamos combinar nada.
    # Directamente pasamos al alineamiento.
    run_muscle_alignment(INPUT_SEQUENCES, ALIGNED_FASTA)
    
    phylo_tree = build_phylogenetic_tree(ALIGNED_FASTA)

    if phylo_tree:
        save_tree_image(phylo_tree, OUTPUT_IMAGE)

    print("\n### Fin del Ejercicio ###")

if __name__ == "__main__":
    main()

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
# --- PASO 4: Guardar el Árbol como Imagen ---
# ==============================================================================

def save_tree_image(tree, output_file: str):
    """
    Guarda un árbol filogenético como una imagen PNG.
    """
    print(f"--- PASO 4: Guardando el árbol en '{output_file}' ---")
    
    tree.root_at_midpoint()
    tree.ladderize()

    fig = plt.figure(figsize=(10, 6), dpi=100)
    axes = fig.add_subplot(1, 1, 1)
    
    Phylo.draw(tree, axes=axes)
    
    plt.title("Árbol Filogenético de Virus Porcinos")
    plt.xlabel("Distancia Evolutiva")
    plt.ylabel("Cepas Virales")
    fig.tight_layout()
    
    plt.savefig(output_file)
    print(f"Gráfico del árbol guardado en '{output_file}'.")

# ==============================================================================
# --- BLOQUE PRINCIPAL DE EJECUCIÓN ---
# ==============================================================================

def main():
    """Función principal que orquesta todos los pasos del ejercicio."""
    print("### Inicio del Ejercicio de Análisis Filogenético ###\n")

    OUTPUT_IMAGE = "arbol_filogenetico.png"
    genome_paths = [os.path.join(GENOMES_DIR, fname) for fname in GENOME_FILES]

    combine_fasta_files(genome_paths, COMBINED_FASTA)
    run_muscle_alignment(COMBINED_FASTA, ALIGNED_FASTA)
    phylo_tree = build_phylogenetic_tree(ALIGNED_FASTA)

    if phylo_tree:
        save_tree_image(phylo_tree, OUTPUT_IMAGE)

    print("\n### Fin del Ejercicio ###")

if __name__ == "__main__":
    main()

