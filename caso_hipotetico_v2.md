# Caso Hipotético Revisado: Brote de Virus Porcinos en una Granja Mexicana

## Escenario Inicial

Octubre de 2025, Jalisco, México. Eres bioinformático en el INIFAP. Una granja en Tepatitlán reporta diarrea severa en lechones (mortalidad 25%, ~500 animales), letargo en adultos y pérdidas de $75,000 USD. Se recolectan 30 muestras (heces, sangre) para NGS. Datos simulados en `genomes/`: `PEDv_genome_new.fasta`, referencias en `all_viruses.fasta`.

## Preparación de Datos

Ensamblaje con SPAdes, anotación con Prokka. Bases locales: `pedv_db`, `pig_viruses_db`. Archivos: `combined_genomes.fasta`, `genomes_to_align.fasta`.

## Análisis

### 1. BLAST (`solucion_filogenia_blast.py`)
- PEDv: 97.5% identidad en S gene (lechones).
- TTSuV2: 92% en adultos, coinfección 40%.
- Output: `gene_sequences_from_blast.fasta`, `aligned_genes.aln`.
- Interpretación: PEDv principal, TTSuV inmunosupresor.

### 2. Filogenia (`ejercicio_filogenia.py`)
- Alineamiento: `aligned_genomes.aln`.
- Árbol: PEDv clado con MX/2022 (bootstrap 95%), TTSuV2 asiático.
- Interpretación: Linaje emergente, escape vacunal.

### 3. PCR In Silico (`analisis_pcr_insilico.py`)
- PEDv: 650 bp en 80% lechones (`pcr_simulation_results.png`).
- TTSuV: 400 bp en 60%.
- Interpretación: Confirma presencia, guía qPCR.

### 4. ML Riesgos (`ml_risk_predictor.py`)
- Riesgo PEDv: 0.92 alto (`pedv_risk_predictor_model.joblib`).
- Regresión: +40% mortalidad (`linear_regression_prediction.py`, `pedv_factores_riesgo.png`).
- Interpretación: Brote de alto impacto.

## Decisiones

1. Inmediatas: Cuarentena, eutanasia lechones, notificación SENASICA.
2. Control: Vacunación, RT-qPCR masiva.
3. Vigilancia: Muestreo adyacentes, rastreo origen.
4. Largo plazo: Bioseguridad, actualizar modelos.

Este caso integra los scripts del taller para simular respuesta real. Revisa `caso_hipotetico_v2.md` y dime si necesitas más ajustes.