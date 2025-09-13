#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ejercicio de Filogenia: Construcción de un Árbol Filogenético de Virus Porcinos

Objetivo:
Construir un árbol filogenético para visualizar las relaciones evolutivas entre 
diferentes cepas de virus porcinos a partir de sus secuencias genómicas.

Pasos a seguir:
1.  Reunir las secuencias de los genomas virales en un solo archivo FASTA.
2.  Realizar un Alineamiento Múltiple de Secuencias (MSA) con MUSCLE.
3.  Construir un árbol filogenético utilizando el método Neighbor-Joining.
4.  Visualizar el árbol en la consola.

Herramientas de Biopython utilizadas:
- SeqIO: Para leer y escribir archivos de secuencias.
- AlignIO: Para leer el resultado del alineamiento.
- Phylo: Para construir y visualizar el árbol filogenético.
- subprocess: Para ejecutar comandos externos como MUSCLE.
"""

import os
import subprocess
from Bio import SeqIO, AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

# ==============================================================================
# --- CONFIGURACIÓN ---
# ==============================================================================

# Directorio que contiene los genomas individuales
GENOMES_DIR = "genomes"

# Lista de los archivos de genomas que se usarán para el árbol
# Asegúrate de que estos archivos existan en el directorio GENOMES_DIR
GENOME_FILES = [
    "PEDv_genome.fasta",
    "TTSuV1_genome.fasta",
    "TTSuV2_genome.fasta",
]

# Nombres de los archivos de salida que se generarán
COMBINED_FASTA = "combined_genomes.fasta"  # Archivo con todos los genomas juntos
ALIGNED_FASTA = "aligned_genomes.aln"      # Archivo de alineamiento (formato clustal)

# ==============================================================================
# --- PASO 1: Combinar Archivos FASTA ---
# ==============================================================================

def combine_fasta_files(file_list: list, output_file: str):
    """
    Combina múltiples archivos FASTA en uno solo.
    
    Args:
        file_list (list): Lista de rutas a los archivos FASTA de entrada.
        output_file (str): Ruta al archivo FASTA de salida.
    """
    print(f"--- PASO 1: Combinando {len(file_list)} archivos FASTA en '{output_file}' ---")
    
    # TAREA DEL ALUMNO: Completa el código para combinar los archivos.
    # Deberás leer las secuencias de cada archivo en `file_list` y 
    # escribirlas en `output_file`.
    
    # Pista: Puedes usar un generador para leer las secuencias de todos los
    # archivos y luego usar SeqIO.write() para escribirlas todas a la vez.

    # ===== INICIO DE LA TAREA 1 =====
    all_sequences = []
    for file_path in file_list:
        # Lee las secuencias del archivo y añádelas a la lista `all_sequences`
        pass # Reemplaza 'pass' con tu código

    # Escribe todas las secuencias combinadas en el archivo de salida
    # SeqIO.write(..., ..., "fasta")
    # ===== FIN DE LA TAREA 1 =====

    print("¡Archivos combinados con éxito!\n")

# ==============================================================================
# --- PASO 2: Realizar Alineamiento Múltiple con MUSCLE ---
# ==============================================================================

def run_muscle_alignment(input_fasta: str, output_alignment: str):
    """
    Ejecuta MUSCLE para realizar un alineamiento múltiple de secuencias.

    Args:
        input_fasta (str): Archivo FASTA con las secuencias a alinear.
        output_alignment (str): Archivo de salida para el alineamiento (formato clustal).
    """
    print(f"--- PASO 2: Alineando secuencias de '{input_fasta}' con MUSCLE ---")
    
    # Para secuencias largas, es mejor usar el algoritmo -super5 de MUSCLE.
    # El comando es: muscle -super5 <archivo_entrada> -output <archivo_salida>
    muscle_cmd = ["muscle", "-super5", input_fasta, "-output", output_alignment]
    
    # TAREA DEL ALUMNO: Ejecuta el comando de MUSCLE.
    # Utiliza el módulo `subprocess` para correr el comando definido en `muscle_cmd`.
    # Asegúrate de capturar la salida y verificar si hay errores.

    # ===== INICIO DE LA TAREA 2 =====
    try:
        # subprocess.run(...)
        print(f"Alineamiento guardado en '{output_alignment}'\n")
    except FileNotFoundError:
        print("Error: El comando 'muscle' no se encontró.")
        print("Por favor, asegúrate de que MUSCLE esté instalado y en el PATH del sistema.")
        exit()
    except Exception as e:
        print(f"Ocurrió un error al ejecutar MUSCLE: {e}")
        exit()
    # ===== FIN DE LA TAREA 2 =====

# ==============================================================================
# --- PASO 3: Construir el Árbol Filogenético ---
# ==============================================================================

def build_phylogenetic_tree(aligned_fasta: str):
    """
    Construye y devuelve un árbol filogenético a partir de un alineamiento.

    Args:
        aligned_fasta (str): Archivo de alineamiento en formato clustal.
    
    Returns:
        Bio.Phylo.BaseTree.Tree: El árbol filogenético construido.
    """
    print(f"--- PASO 3: Construyendo árbol desde '{aligned_fasta}' ---")
    
    # TAREA DEL ALUMNO: Lee el alineamiento y construye el árbol.
    # 1. Lee el archivo de alineamiento con AlignIO en formato "fasta".
    # 2. Calcula la matriz de distancias usando `DistanceCalculator`.
    # 3. Construye el árbol con `DistanceTreeConstructor` y el método Neighbor-Joining ('nj').

    # ===== INICIO DE LA TAREA 3 =====
    # 1. Lee el alineamiento
    aln = None # Usa AlignIO.read(...)

    # 2. Calcula la matriz de distancias (modelo 'identity')
    calculator = DistanceCalculator('identity')
    dist_matrix = None # Usa calculator.get_distance(aln)

    # 3. Construye el árbol
    constructor = DistanceTreeConstructor()
    tree = None # Usa constructor.nj(dist_matrix)
    # ===== FIN DE LA TAREA 3 =====

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
    
    # TAREA DEL ALUMNO: Muestra el árbol en la consola.
    # La librería Phylo puede "dibujar" el árbol en formato de texto.
    
    # ===== INICIO DE LA TAREA 4 =====
    # Phylo.draw_ascii(...)
    # ===== FIN DE LA TAREA 4 =====

# ==============================================================================
# --- BLOQUE PRINCIPAL DE EJECUCIÓN ---
# ==============================================================================

def main():
    """Función principal que orquesta todos los pasos del ejercicio."""
    print("### Inicio del Ejercicio de Análisis Filogenético ###\n")

    # Construye la lista completa de rutas a los archivos de genomas
    genome_paths = [os.path.join(GENOMES_DIR, fname) for fname in GENOME_FILES]

    # --- Ejecutar Pasos ---
    # 1. Combinar archivos FASTA
    combine_fasta_files(genome_paths, COMBINED_FASTA)

    # 2. Realizar alineamiento
    run_muscle_alignment(COMBINED_FASTA, ALIGNED_FASTA)

    # 3. Construir el árbol
    phylo_tree = build_phylogenetic_tree(ALIGNED_FASTA)

    # 4. Visualizar el árbol
    if phylo_tree:
        visualize_tree(phylo_tree)

    print("\n### Fin del Ejercicio ###")

if __name__ == "__main__":
    main()
