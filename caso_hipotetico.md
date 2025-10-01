# Caso Hipotético: Brote de Virus Porcinos en una Granja Mexicana (Versión Revisada)

## Escenario Inicial

Es octubre de 2025, y trabajas como bioinformático en el Instituto Nacional de Investigaciones Forestales, Agrícolas y Pecuarias (INIFAP) en México. Una# Caso Hipotético: Brote de Virus Porcinos en una Granja Mexicana (Versión Revisada)

## Escenario Inicial

Es octubre de 2025, y trabajas como especialista en bioinformática en el Instituto Nacional de Investigaciones Forestales, Agrícolas y Pecuarias (INIFAP) en México. Una granja porcina en el municipio de Tep

En el laboratorio, se extrae el ADN/ARN viral de las muestras utilizando kits comerciales. Se realiza secuenciación de nueva generación (NGS) para obtener lecturas genómicas crudas. Para este taller, asumimos que las secuencias procesadas están disponibles en archivos como `PEDv_genome_new.fasta`, `TTSuV1_genome_new.fasta` y `TTSuV2_genome_new.fasta` en la carpeta `genomes/`, representando genomas virales ensamblados de las muestras.

Se combinan con bases de datos conocidas: `all_viruses.fasta` incluye genomas de referencia de PEDv (virus de la diarrea epidémica porcina), TTSuV1 y TTSuV2 (torque teno sus virus), así como otros virus porcinos comunes.

## Análisis Realizados

### 1. Identificación de Secuencias vía BLAST
Utilizando el script implícito en `solucion_filogenia_blast.py`, se ejecuta BLAST contra bases de datos locales como `pedv_db` y `pig_viruses_db`. Resultados simulados:
- En muestras de lechones: Alto match (98% identidad) con el genoma de PEDv en la región S (proteína de espiga), indicando un aislado virulento similar a brotes de 2013-2014 en EE.UU.
- En adultos: Detección de TTSuV2 (95% identidad) en genes no estructurales, y co-infección con TTSuV1 en el 30% de las muestras.
- Archivo de salida: `gene_sequences_from_blast.fasta` contiene las secuencias de genes identificados, alineadas en `aligned_genes.aln`.

**Interpretación**: Confirma presencia de PEDv como agente principal de la diarrea, con TTSuV como coinfectantes potencialmente inmunosupresores, aumentando la severidad.

### 2. Alineamiento Múltiple y Análisis Filogenético
Se alinean los genomas nuevos con referencias usando herramientas como MAFFT, generando `aligned_genomes.aln` y `aligned_viruses.afa`. Luego, se construye un árbol filogenético con `ejercicio_filogenia.py` o `solucion_filogenia.py` (usando PhyML o similar).

Resultados simulados:
- El aislado PEDv de la granja forma un clado con variantes US/IN/2014, con soporte bootstrap >90%, sugiriendo introducción reciente vía importaciones o aves silvestres.
- TTSuV2 se agrupa con linajes asiáticos, indicando circulación endémica.
- Visualización: El árbol muestra que el PEDv es más cercano a cepas epidémicas que a atenuadas vacunales.

**Interpretación**: El PEDv es un linaje emergente, no controlado por vacunas estándar, requiriendo respuesta rápida.

### 3. Simulación de PCR In Silico
Con `analisis_pcr_insilico.py`, se simula PCR para primers específicos de PEDv (S gene) y TTSuV (ORF1/2). Resultados en `pcr_simulation_results.png`:
- Amplicones positivos para PEDv en todas las muestras de lechones (tamaño ~800 bp), confirmando carga viral alta.
- TTSuV detectado en 70% de muestras adultas, con bandas dobles indicando genotipos mixtos.
- No se detectan otros enterovirus comunes.

**Interpretación**: Valida la presencia viral y sugiere coinfecciones que complican el cuadro clínico. En un laboratorio real, esto guiaría pruebas confirmatorias.

### 4. Predicción de Riesgos con Machine Learning
Usando `ml_risk_predictor.py` y el modelo `pedv_risk_predictor_model.joblib`, se predicen riesgos basados en features como identidad genómica, carga viral simulada y factores clínicos (edad, mortalidad).
- Predicción para PEDv: Riesgo alto (0.85/1.0), con regresión lineal en `linear_regression_prediction.py` proyectando un aumento del 30% en mortalidad si no se interviene.
- Para TTSuV2: Carga viral media (ver `ttsuv2_carga_viral_analisis.png`), pero riesgo de inmunosupresión (0.65/1.0).
- Gráficos: `pedv_factores_riesgo.png` muestra que la edad <7 días y coinfección elevan el riesgo.

**Interpretación**: El modelo indica un brote de alto impacto, con potencial para propagación a granjas vecinas.

## Toma de Decisiones Basadas en Resultados

Integrando todos los análisis:

1. **Acciones Inmediatas**:
   - Aislamiento de la granja: Cuarentena total, desinfección y restricción de movimientos. (Basado en match BLAST y PCR positivo).
   - Sacrificio humanitario de lechones afectados para reducir sufrimiento y reservorio viral.

2. **Medidas de Control**:
   - Vacunación de emergencia: Usar vacunas atenuadas contra PEDv, aunque el linaje filogenético sugiere posible escape; monitorear eficacia.
   - Pruebas adicionales: qPCR en tiempo real para cuantificar carga viral en todos los animales.

3. **Vigilancia y Prevención**:
   - Reporte a SENASICA (Servicio Nacional de Sanidad, Inocuidad y Calidad Agroalimentaria): Notificar brote para alertas nacionales.
   - Muestreo en granjas adyacentes: Usar el modelo ML para priorizar basadas en factores de riesgo (e.g., densidad porcina).
   - Investigación epidemiológica: Rastrear origen (e.g., piensos contaminados o vectores).

4. **Recomendaciones a Largo Plazo**:
   - Mejora bioseguridad: Implementar protocolos de all-in/all-out y ventilación.
   - Actualización de modelos: Incorporar este brote al dataset de ML para refinar predicciones futuras.
   - Colaboración: Compartir secuencias en GenBank para vigilancia global.

Este caso hipotético simula un brote real, destacando cómo los análisis bioinformáticos guían decisiones críticas en salud animal, previniendo pérdidas mayores y protegiendo la industria porcina. En el taller, los estudiantes pueden replicar estos pasos con los datos proporcionados para "resolver" el caso.