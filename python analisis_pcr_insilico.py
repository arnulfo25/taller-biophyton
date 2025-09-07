# -*- coding: utf-8 -*-
"""
Análisis In Silico de Especificidad y Amplificación de Cebadores para PCR.

Este script realiza dos funciones principales:
1.  Verifica la especificidad de cebadores individuales contra una base de datos 
    de genomas virales usando BLASTn.
2.  Simula la amplificación por PCR para pares de cebadores en genomas objetivo y 
    no objetivo para predecir el tamaño del amplicón y detectar posibles 
    amplificaciones inespecíficas.

Autor: [Tu Nombre]
Fecha: [Fecha de Hoy]
"""

import os
import subprocess
from Bio.Seq import Seq
from Bio import SeqIO

# ==============================================================================
# --- CONFIGURACIÓN Y DATOS DE CEBADORES ---
# ==============================================================================

# --- Constantes de Configuración ---
BLAST_DB = "pig_viruses_db"  # Nombre de la base de datos BLAST a utilizar
GENOMES_DIR = "genomes"      # Directorio donde están los archivos FASTA de los genomas

def reverse_complement(seq):
    """Devuelve el complemento reverso de una secuencia de ADN."""
    return str(Seq(seq).reverse_complement())

# --- Definición de Cebadores, Sondas y Ensayos ---

# Porcine Epidemic Diarrhea Virus (PEDv)
# !! IMPORTANTE: Reemplazar "..." con las secuencias reales.
PEDv_primers_probes = {
    "PED-N_Forward": "...",
    "PED-N_Reverse": "...",
    "PED-N_Probe": "...",
    "PED-S_Forward": "...",
    "PED-S_Reverse": "...",
    "PED-S_Probe": "...",
}

# Torque Teno Sus Virus (TTSuV)
# !! IMPORTANTE: Reemplazar "..." con las secuencias reales.
TTSuV_primers = {
    "TTSuV1F": "...",
    "TTSuV1R": "...",
    "TTSuV2F": "...",
    "TTSuV2R": "...",
}

# Tamaños esperados de los amplicones (en pares de bases)
# !! IMPORTANTE: Reemplazar 0 con los tamaños reales.
expected_amplicons = {
    "PED-N": 150,   # Ejemplo, ajustar al valor real
    "PED-S": 130,   # Ejemplo, ajustar al valor real
    "TTSuV1": 210,  # Ejemplo, ajustar al valor real
    "TTSuV2": 180,  # Ejemplo, ajustar al valor real
}

# Genomas objetivo para cada ensayo
# !! IMPORTANTE: Usar los nombres de archivo correctos que están en la carpeta GENOMES_DIR.
target_genomes = {
    "PED-N": ["PEDv_genome.fasta"],
    "PED-S": ["PEDv_genome.fasta"],
    "TTSuV1": ["TTSuV1_genome.fasta"],
    "TTSuV2": ["TTSuV2_genome.fasta"],
}


# ==============================================================================
# --- FUNCIONES DE ANÁLISIS (BLAST y PCR) ---
# ==============================================================================

def run_blastn(query_seq, db_name):
    """Ejecuta blastn con la secuencia de consulta contra la base de datos."""
    query_fasta = f">query\n{query_seq}"
    
    cmd = [
        "blastn",
        "-task", "blastn-short",
        "-query", "-",
        "-db", db_name,
        "-outfmt", "6 sseqid pident length qstart qend sstart send evalue",
        "-word_size", "7",
        "-dust", "no",
        "-soft_masking", "false",
        "-perc_identity", "85"
    ]
    
    try:
        process = subprocess.run(cmd, input=query_fasta.encode(), capture_output=True, check=True, text=True)
        return process.stdout.strip().split('\n')
    except subprocess.CalledProcessError as e:
        print(f"Error al ejecutar BLAST para la secuencia '{query_seq[:10]}...': {e}")
        print(f"Stderr: {e.stderr}")
        return []

def analyze_blast_hits(primer_name, primer_seq, blast_results, seq_id_to_genome_file, all_target_genomes):
    """Analiza e imprime los resultados de BLAST para un cebador."""
    print(f"\n--- Análisis de BLAST para {primer_name} ({primer_seq}) ---")
    
    if not blast_results or blast_results == ['']:
        print("  No se encontraron hits significativos.")
        return

    hits_by_genome = {}
    for line in blast_results:
        if not line: continue
        parts = line.split('\t')
        sseqid, pident, length, qstart, qend, sstart, send, evalue = parts
        
        genome_file = seq_id_to_genome_file.get(sseqid.split()[0])
        if not genome_file:
            print(f"  Advertencia: No se pudo mapear el hit '{sseqid}' a un archivo de genoma.")
            continue

        if genome_file not in hits_by_genome:
            hits_by_genome[genome_file] = []
        
        hits_by_genome[genome_file].append({
            "sseqid": sseqid, "pident": float(pident), "length": int(length),
            "qstart": int(qstart), "qend": int(qend), "sstart": int(sstart),
            "send": int(send), "evalue": float(evalue)
        })

    for genome_file, hits in hits_by_genome.items():
        is_target = genome_file in all_target_genomes
        print(f"  Resultados en {genome_file} (Genoma Objetivo: {is_target}):")
        for hit in hits:
            print(f"    - ID: {hit['sseqid']}, Identidad: {hit['pident']}%, Longitud: {hit['length']}, E-value: {hit['evalue']:.2e}")
            if hit['pident'] < 100:
                print(f"      **Advertencia**: Posible unión no específica o con mismatches.")

def simulate_pcr_amplification(forward_primer_seq, reverse_primer_seq, genome_record, expected_size):
    """Simula la amplificación en una secuencia de genoma específica."""
    genome_seq = genome_record.seq.upper()
    
    f_hits = [i for i in range(len(genome_seq) - len(forward_primer_seq) + 1)
              if genome_seq[i:i+len(forward_primer_seq)] == forward_primer_seq]
    
    r_rc_seq = reverse_complement(reverse_primer_seq)
    r_hits = [i for i in range(len(genome_seq) - len(r_rc_seq) + 1)
              if genome_seq[i:i+len(r_rc_seq)] == r_rc_seq]

    predicted_amplicons = []
    for f_start in f_hits:
        for r_start in r_hits:
            if f_start < r_start:
                amplicon_end = r_start + len(r_rc_seq)
                predicted_size = amplicon_end - f_start
                
                if expected_size:
                    min_expected = expected_size * 0.9
                    max_expected = expected_size * 1.1
                    if min_expected <= predicted_size <= max_expected:
                        predicted_amplicons.append((f_start, amplicon_end, predicted_size))
                else:
                    predicted_amplicons.append((f_start, amplicon_end, predicted_size))

    if predicted_amplicons:
        print(f"  Amplificación PREDICHA en {genome_record.id}:")
        for amp_start, amp_end, amp_size in predicted_amplicons:
            print(f"    - Rango: {amp_start}-{amp_end}, Tamaño: {amp_size} pb. (Esperado: {expected_size or 'N/A'} pb)")
            if expected_size and not (expected_size * 0.9 <= amp_size <= expected_size * 1.1):
                 print(f"      **Advertencia**: El tamaño predicho ({amp_size} pb) difiere del esperado ({expected_size} pb).")
    else:
        print(f"  No se predice amplificación en {genome_record.id}.")


# ==============================================================================
# --- BLOQUE PRINCIPAL DE EJECUCIÓN ---
# ==============================================================================

if __name__ == "__main__":
    # --- 1. Carga y preparación de datos ---
    print("--- Cargando genomas y preparando datos ---")
    if not os.path.exists(GENOMES_DIR):
        print(f"Error: El directorio de genomas '{GENOMES_DIR}' no existe.")
        exit()

    genome_files_map = {}
    seq_id_to_genome_file = {}
    for filename in os.listdir(GENOMES_DIR):
        if filename.endswith((".fasta", ".fa", ".fna")):
            filepath = os.path.join(GENOMES_DIR, filename)
            try:
                records = list(SeqIO.parse(filepath, "fasta"))
                genome_files_map[filename] = records
                for record in records:
                    seq_id_to_genome_file[record.id.split()[0]] = filename
            except Exception as e:
                print(f"Error al cargar o procesar {filepath}: {e}")

    all_primers = {**PEDv_primers_probes, **TTSuV_primers}
    all_target_genomes = {g for g_list in target_genomes.values() for g in g_list}

    # --- 2. Verificación de Especificidad con BLAST ---
    print("\n--- Verificación de Especificidad Individual de Cebadores con BLAST ---")
    for primer_name, primer_seq in all_primers.items():
        if "Probe" in primer_name: continue # Omitir sondas
        
        # Verificar secuencia original
        blast_results = run_blastn(primer_seq, BLAST_DB)
        analyze_blast_hits(primer_name, primer_seq, blast_results, seq_id_to_genome_file, all_target_genomes)
        
        # Verificar complemento reverso para uniones en la otra hebra
        rev_comp_seq = reverse_complement(primer_seq)
        blast_results_rc = run_blastn(rev_comp_seq, BLAST_DB)
        analyze_blast_hits(f"{primer_name}_RC", rev_comp_seq, blast_results_rc, seq_id_to_genome_file, all_target_genomes)

    # --- 3. Simulación de Amplificación por Pares ---
    print("\n--- Simulación de Amplificación por Pares de Cebadores ---")
    assay_defs = [
        {"name": "PED-N", "fw": "PED-N_Forward", "rv": "PED-N_Reverse", "db": PEDv_primers_probes},
        {"name": "PED-S", "fw": "PED-S_Forward", "rv": "PED-S_Reverse", "db": PEDv_primers_probes},
        {"name": "TTSuV1", "fw": "TTSuV1F", "rv": "TTSuV1R", "db": TTSuV_primers},
        {"name": "TTSuV2", "fw": "TTSuV2F", "rv": "TTSuV2R", "db": TTSuV_primers},
    ]

    for assay in assay_defs:
        assay_name = assay["name"]
        fw_seq = assay["db"][assay["fw"]]
        rv_seq = assay["db"][assay["rv"]]
        
        print(f"\n--- Ensayo de Simulación para: {assay_name} ---")
        
        # Iterar sobre todos los genomas cargados
        for genome_file, records in genome_files_map.items():
            is_target = genome_file in target_genomes.get(assay_name, [])
            expected_size = expected_amplicons.get(assay_name) if is_target else None
            
            print(f"\nSimulando en '{genome_file}' (Objetivo: {is_target})")
            for record in records:
                simulate_pcr_amplification(fw_seq, rv_seq, record, expected_size)

    print("\n--- Análisis completado. ---")
