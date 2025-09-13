#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Análisis In Silico de Especificidad y Amplificación de Cebadores para PCR.

Este script realiza dos funciones principales:
1.  Verifica la especificidad de cebadores individuales contra una base de datos 
    de genomas virales usando BLASTn.
2.  Simula la amplificación por PCR para pares de cebadores en genomas objetivo y 
    no objetivo para predecir el tamaño del amplicón y detectar posibles 
    amplificaciones inespecíficas.
3.  Genera un gráfico de barras para visualizar y comparar los resultados de la 
    simulación de PCR.

Autor: [Tu Nombre]
Fecha: [Fecha de Hoy]
"""

import os
import subprocess
from Bio.Seq import Seq
from Bio import SeqIO
from collections import defaultdict
import logging
from typing import List, Dict, Set, Tuple, Optional
import matplotlib.pyplot as plt
import numpy as np

# ==============================================================================
# --- CONFIGURACIÓN Y DATOS DE CEBADORES ---
# ==============================================================================

# --- Constantes de Configuración ---
BLAST_DB = "pig_viruses_db"  # Nombre de la base de datos BLAST a utilizar
GENOMES_DIR = "genomes"      # Directorio donde están los archivos FASTA de los genomas
MAX_MISMATCHES = 1           # Máximo de mismatches permitidos en la simulación de PCR
BLAST_PERC_IDENTITY = 85     # Mínimo porcentaje de identidad para hits en BLAST
AMPLICON_SIZE_TOLERANCE = 0.1 # Tolerancia para el tamaño del amplicón (e.g., 0.1 para +/- 10%)
BLAST_OUTFMT = "6 sseqid pident length qstart qend sstart send evalue" # Formato de salida para BLAST

def reverse_complement(seq: str) -> str:
    """Devuelve el complemento reverso de una secuencia de ADN."""
    return str(Seq(seq).reverse_complement())

# --- Definición de Cebadores y Sondas ---
ALL_PRIMERS = {
    # Porcine Epidemic Diarrhea Virus (PEDv)
    "PED-N_Forward": "CGCAAAGACTGAACCCACTAAC",
    "PED-N_Reverse": "TTGCCTCTGTTGTTACTTGGAG",
    "PED-N_Probe": "AGCTTTGAGCCTGCATCTATTGGTCTC",
    "PED-S_Forward": "GCTATTGGTGGCATTACTCCTG",
    "PED-S_Reverse": "CCATCAACATCTACAGGTCGTT",
    "PED-S_Probe": "TGGCTCTAACTGCTCATGCTATGC",
    # Torque Teno Sus Virus (TTSuV)
    "TTSuV1F": "GTAGATGCCACAGTCGATGAAC",
    "TTSuV1R": "CCTACGTCAAGTGGAACACAGG",
    "TTSuV2F": "ATGGCCAGCCAGACAATAAAGA",
    "TTSuV2R": "GCAACACCCAGTAGGTTGACTT",
}

# Estructura de datos unificada para definir los ensayos de PCR
ASSAYS = {
    "PED-N": {
        "forward_primer": "PED-N_Forward",
        "reverse_primer": "PED-N_Reverse",
        "target_genomes": ["PEDv_genome.fasta"],
        "expected_size": 150,
    },
    "PED-S": {
        "forward_primer": "PED-S_Forward",
        "reverse_primer": "PED-S_Reverse",
        "target_genomes": ["PEDv_genome.fasta"],
        "expected_size": 130,
    },
    "TTSuV1": {
        "forward_primer": "TTSuV1F",
        "reverse_primer": "TTSuV1R",
        "target_genomes": ["TTSuV1_genome.fasta"],
        "expected_size": 210,
    },
    "TTSuV2": {
        "forward_primer": "TTSuV2F",
        "reverse_primer": "TTSuV2R",
        "target_genomes": ["TTSuV2_genome.fasta"],
        "expected_size": 180,
    },
}



# ==============================================================================
# --- FUNCIONES DE ANÁLISIS (BLAST y PCR) ---
# ==============================================================================

def run_blastn(query_seq: str, db_name: str) -> List[str]:
    """Ejecuta blastn con la secuencia de consulta contra la base de datos."""
    query_fasta = f">query\n{query_seq}"
    
    cmd = [
        "blastn",
        "-task", "blastn-short",
        "-query", "-",
        "-db", db_name,
        "-outfmt", BLAST_OUTFMT,
        "-word_size", "7",
        "-dust", "no",
        "-soft_masking", "false",
        "-perc_identity", str(BLAST_PERC_IDENTITY)
    ]
    
    try:
        process = subprocess.run(
            cmd, input=query_fasta, capture_output=True, check=True, text=True
        )
        return process.stdout.strip().split('\n')
    except FileNotFoundError:
        logging.error(f"\nError: El comando '{cmd[0]}' no se encontró.")
        logging.error("Por favor, asegúrate de que NCBI BLAST+ esté instalado y en el PATH del sistema.")
        logging.error("En sistemas Debian/Ubuntu, puedes instalarlo con: sudo apt-get install ncbi-blast+")
        exit()
    except subprocess.CalledProcessError as e:
        logging.error(f"Error al ejecutar BLAST para la secuencia '{query_seq[:10]}...': {e}")
        logging.error(f"Stderr: {e.stderr}")
        if "BLAST Database error: No alias or index file found for nucleotide database" in e.stderr:
            logging.warning(f"Sugerencia: La base de datos BLAST '{db_name}' no existe. Asegúrate de crearla con 'makeblastdb'.")
        return []

def analyze_blast_hits(primer_name: str, primer_seq: str, blast_results: List[str], seq_id_to_genome_file: Dict[str, str], all_target_genomes: Set[str]) -> None:
    """Analiza e imprime los resultados de BLAST para un cebador."""
    logging.info(f"\n--- Análisis de BLAST para {primer_name} ({primer_seq}) ---")
    
    if not blast_results or blast_results == ['']:
        logging.info("  No se encontraron hits significativos.")
        return

    hits_by_genome = defaultdict(list)
    for line in blast_results:
        if not line: continue
        parts = line.split('\t')
        sseqid, pident, length, qstart, qend, sstart, send, evalue = parts
        
        genome_file = seq_id_to_genome_file.get(sseqid.split()[0])
        if not genome_file:
            logging.warning(f"  Advertencia: No se pudo mapear el hit '{sseqid}' a un archivo de genoma.")
            continue
        
        hits_by_genome[genome_file].append({
            "sseqid": sseqid, "pident": float(pident), "length": int(length),
            "qstart": int(qstart), "qend": int(qend), "sstart": int(sstart),
            "send": int(send), "evalue": float(evalue)
        })

    for genome_file, hits in hits_by_genome.items():
        is_target = genome_file in all_target_genomes
        logging.info(f"  Resultados en {genome_file} (Genoma Objetivo: {is_target}):")
        for hit in hits:
            strand = "forward" if hit['sstart'] < hit['send'] else "reverse"
            logging.info(f"    - ID: {hit['sseqid']}, Identidad: {hit['pident']}%, Longitud: {hit['length']}, Hebra: {strand}, E-value: {hit['evalue']:.2e}")
            if hit['pident'] < 100:
                logging.warning(f"      **Advertencia**: Posible unión no específica o con mismatches.")

def find_with_mismatches(sequence: str, primer: str, max_mismatches: int = 0) -> List[int]:
    """Encuentra los sitios de unión de un cebador permitiendo mismatches."""
    hits = []
    len_primer = len(primer)
    for i in range(len(sequence) - len_primer + 1):
        window = sequence[i:i+len_primer]
        mismatches = sum(1 for j in range(len_primer) if window[j] != primer[j])
        if mismatches <= max_mismatches:
            hits.append(i)
    return hits

def simulate_pcr_amplification(forward_primer_seq: str, reverse_primer_rc_seq: str, genome_record: SeqIO.SeqRecord, expected_size: Optional[int], is_target: bool, max_mismatches: int = 0) -> List[Tuple[int, int, int]]:
    """Simula la amplificación en una secuencia de genoma específica y devuelve los amplicones."""
    genome_seq = genome_record.seq.upper()
    
    f_hits = find_with_mismatches(genome_seq, forward_primer_seq, max_mismatches)
    r_hits = find_with_mismatches(genome_seq, reverse_primer_rc_seq, max_mismatches)

    predicted_amplicons = []
    for f_start in f_hits:
        for r_start in r_hits:
            if f_start < r_start:
                amplicon_end = r_start + len(reverse_primer_rc_seq)
                predicted_size = amplicon_end - f_start
                
                if expected_size:
                    min_expected = expected_size * (1 - AMPICON_SIZE_TOLERANCE)
                    max_expected = expected_size * (1 + AMPICON_SIZE_TOLERANCE)
                    if min_expected <= predicted_size <= max_expected:
                        predicted_amplicons.append((f_start, amplicon_end, predicted_size))
                else:
                    predicted_amplicons.append((f_start, amplicon_end, predicted_size))

    if predicted_amplicons:
        if not is_target:
            logging.warning(f"  **ADVERTENCIA: Amplificación INESPECÍFICA PREDICHA en {genome_record.id} (No-Objetivo):")
        else:
            logging.info(f"  Amplificación PREDICHA en {genome_record.id}:")
        for amp_start, amp_end, amp_size in predicted_amplicons:
            logging.info(f"    - Rango: {amp_start}-{amp_end}, Tamaño: {amp_size} pb. (Esperado: {expected_size or 'N/A'} pb)")
            if expected_size and not (expected_size * (1 - AMPICON_SIZE_TOLERANCE) <= amp_size <= expected_size * (1 + AMPICON_SIZE_TOLERANCE)):
                 logging.warning(f"      **Advertencia**: El tamaño predicho ({amp_size} pb) difiere del esperado ({expected_size} pb).")
    else:
        if is_target:
            logging.info(f"  No se predice amplificación en {genome_record.id}.")

    return predicted_amplicons

# ==============================================================================
# --- FUNCIÓN DE VISUALIZACIÓN ---
# ==============================================================================

def generate_pcr_results_graph(results: List[Dict], output_filename: str = "pcr_simulation_results.png"):
    """Genera un gráfico de barras y puntos con los resultados de la simulación de PCR."""
    logging.info(f"\n--- Generando gráfico de resultados de la simulación de PCR ---")

    assay_names = sorted(list(set(r['assay'] for r in results)))
    n_assays = len(assay_names)
    
    if not n_assays:
        logging.warning("No hay resultados para graficar.")
        return

    fig, ax = plt.subplots(figsize=(max(10, n_assays * 2), 8))
    
    assay_map = {name: i for i, name in enumerate(assay_names)}
    
    # 1. Plot expected sizes as bars
    expected_sizes = [next(r['expected'] for r in results if r['assay'] == name) for name in assay_names]
    ax.bar(assay_names, expected_sizes, color='skyblue', alpha=0.7, label='Tamaño Esperado')

    # 2. Plot predicted sizes as scatter points
    for res in results:
        if res['predicted']:
            assay_idx = assay_map[res['assay']]
            # Add jitter for visibility if multiple predictions exist for the same assay
            jitter = np.random.normal(0, 0.05)
            for pred_size in res['predicted']:
                if res['is_target']:
                    ax.plot(assay_idx + jitter, pred_size, 'go', markersize=10, alpha=0.8, label='Predicción (Objetivo)')
                else:
                    ax.plot(assay_idx + jitter, pred_size, 'rX', markersize=10, markeredgewidth=2, label='Predicción (No-Objetivo)')

    ax.set_xlabel("Ensayo de PCR", fontsize=12)
    ax.set_ylabel("Tamaño del Amplicón (pb)", fontsize=12)
    ax.set_title("Comparación de Tamaños de Amplicones: Esperado vs. Predicho", fontsize=14, weight='bold')
    ax.set_xticks(range(n_assays))
    ax.set_xticklabels(assay_names, rotation=45, ha="right")
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Create a clean legend without duplicate labels
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc='upper right')
    
    plt.tight_layout()
    try:
        plt.savefig(output_filename)
        logging.info(f"--- Gráfico de resultados guardado en: {output_filename} ---")
    except Exception as e:
        logging.error(f"Error al guardar el gráfico: {e}")


# ==============================================================================
# --- BLOQUE PRINCIPAL DE EJECUCIÓN ---
# ==============================================================================

def load_and_prepare_data(genomes_dir: str) -> Tuple[Dict[str, List[SeqIO.SeqRecord]], Dict[str, str]]:
    """Carga los genomas desde los archivos FASTA y crea mapeos de IDs."""
    genome_files_map = {}
    seq_id_to_genome_file = {}
    for filename in os.listdir(genomes_dir):
        if filename.endswith((".fasta", ".fa", ".fna")):
            filepath = os.path.join(genomes_dir, filename)
            try:
                records = list(SeqIO.parse(filepath, "fasta"))
                genome_files_map[filename] = records
                for record in records:
                    seq_id_to_genome_file[record.id.split()[0]] = filename
            except Exception as e:
                logging.error(f"Error al cargar o procesar {filepath}: {e}")
    return genome_files_map, seq_id_to_genome_file

def perform_blast_analysis(all_primers: Dict[str, str], seq_id_to_genome_file: Dict[str, str], all_target_genomes: Set[str]) -> None:
    """Ejecuta el análisis de especificidad de BLAST para todos los cebadores."""
    logging.info("\n" + "="*60)
    logging.info("--- Verificación de Especificidad Individual de Cebadores con BLAST ---")
    logging.info("="*60)
    for primer_name, primer_seq in all_primers.items():
        if "Probe" in primer_name: continue # Omitir sondas
        
        # Verificar secuencia original
        blast_results = run_blastn(primer_seq, BLAST_DB)
        analyze_blast_hits(primer_name, primer_seq, blast_results, seq_id_to_genome_file, all_target_genomes)
        
        # Verificar complemento reverso para uniones en la otra hebra
        rev_comp_seq = reverse_complement(primer_seq)
        blast_results_rc = run_blastn(rev_comp_seq, BLAST_DB)
        analyze_blast_hits(f"{primer_name}_RC", rev_comp_seq, blast_results_rc, seq_id_to_genome_file, all_target_genomes)

def perform_pcr_simulation(assays: Dict, all_primers: Dict[str, str], genome_files_map: Dict[str, List[SeqIO.SeqRecord]]) -> List[Dict]:
    """Ejecuta la simulación de amplificación por PCR y recopila los resultados."""
    logging.info("\n" + "="*60)
    logging.info("--- Simulación de Amplificación por Pares de Cebadores ---")
    logging.info("="*60)

    all_assay_results = []

    for assay_name, assay_data in assays.items():
        fw_primer_name = assay_data["forward_primer"]
        rv_primer_name = assay_data["reverse_primer"]
        fw_seq = all_primers[fw_primer_name]
        rv_seq = all_primers[rv_primer_name]
        rv_rc_seq = reverse_complement(rv_seq)
        
        logging.info(f"\n--- Ensayo de Simulación para: {assay_name} ---")
        
        for genome_file, records in genome_files_map.items():
            is_target = genome_file in assay_data["target_genomes"]
            expected_size = assay_data["expected_size"] if is_target else None
            
            logging.info(f"\nSimulando en '{genome_file}' (Objetivo: {is_target})")
            
            predicted_amplicons_for_genome = []
            for record in records:
                amplicons = simulate_pcr_amplification(fw_seq, rv_rc_seq, record, expected_size, is_target, MAX_MISMATCHES)
                predicted_amplicons_for_genome.extend(amplicons)

            all_assay_results.append({
                "assay": assay_name,
                "genome": genome_file,
                "is_target": is_target,
                "expected": assay_data["expected_size"],
                "predicted": [amp[2] for amp in predicted_amplicons_for_genome]
            })
            
    return all_assay_results

def main() -> None:
    """Función principal que orquesta el análisis."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(message)s',
    )
    logging.info("--- Cargando genomas y preparando datos ---")
    if not os.path.exists(GENOMES_DIR):
        logging.error(f"Error: El directorio de genomas '{GENOMES_DIR}' no existe.")
        return

    genome_files_map, seq_id_to_genome_file = load_and_prepare_data(GENOMES_DIR)
    if not genome_files_map:
        logging.error(f"Error: No se encontraron archivos de genoma válidos en el directorio '{GENOMES_DIR}'.")
        logging.error("Por favor, asegúrate de que la carpeta contenga archivos .fasta, .fa o .fna.")
        return

    all_target_genomes = {g for assay_data in ASSAYS.values() for g in assay_data["target_genomes"]}
    
    perform_blast_analysis(ALL_PRIMERS, seq_id_to_genome_file, all_target_genomes)
    
    pcr_results = perform_pcr_simulation(ASSAYS, ALL_PRIMERS, genome_files_map)
    
    generate_pcr_results_graph(pcr_results)
    
    logging.info("\n--- Análisis completado. ---")

if __name__ == "__main__":
    main()
