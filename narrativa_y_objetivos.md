# Narrativa y Objetivos del Taller de Bioinformática en Virus Porcinos

## Narrativa del Proyecto

Este taller de bioinformática se centra en el análisis de secuencias genómicas de virus que afectan al ganado porcino, con énfasis en el virus de la diarrea epidémica porcina (PEDv) y los virus Torque teno sus (TTSuV1 y TTSuV2). En el contexto de la salud animal y la producción porcina, la detección temprana y el análisis de riesgos virales son cruciales para prevenir brotes epidémicos que pueden causar pérdidas económicas significativas.

El proyecto simula un escenario real donde los estudiantes actúan como investigadores en un laboratorio veterinario. Utilizando herramientas bioinformáticas, se exploran secuencias genómicas de virus aislados de muestras porcinas. A través de alineamientos de secuencias, análisis filogenéticos, simulaciones de PCR in silico y modelos de machine learning para predicción de riesgos, los participantes aprenden a interpretar datos genéticos para apoyar decisiones en salud animal.

Los datos incluyen genomas virales alineados, secuencias de genes identificados vía BLAST, y resultados de predicciones de carga viral y factores de riesgo. Este enfoque práctico integra conceptos de biología molecular, computación y estadística, preparando a los estudiantes para aplicaciones reales en biotecnología y veterinaria.

## Objetivos del Taller

1. **Entender los fundamentos de bioinformática en virología porcina**: Familiarizarse con el manejo de secuencias genómicas, alineamientos múltiples y herramientas como BLAST para identificación de genes virales.

2. **Realizar análisis filogenéticos**: Construir árboles filogenéticos a partir de alineamientos para determinar relaciones evolutivas entre virus porcinos y aislamientos conocidos.

3. **Simular experimentos moleculares**: Ejecutar PCR in silico para detectar la presencia de secuencias virales específicas en muestras genómicas.

4. **Aplicar machine learning en predicción de riesgos**: Desarrollar y utilizar modelos para predecir la carga viral y factores de riesgo basados en datos genómicos y clínicos.

5. **Interpretar resultados y tomar decisiones informadas**: Aprender a analizar outputs de los análisis para recomendar acciones como cuarentenas, pruebas adicionales o estrategias de control epidémico.

## Cómo Interpretar Resultados y Tomar Decisiones

Los estudiantes deben seguir un flujo lógico para interpretar los resultados de los análisis y derivar decisiones basadas en evidencia. A continuación, se detalla el proceso paso a paso:

### 1. **Análisis de Alineamientos y BLAST**
   - **Resultados clave**: Archivos como `aligned_genes.aln`, `gene_sequences_from_blast.fasta`.
   - **Interpretación**: Verifica la similitud entre secuencias de genes virales identificados y bases de datos conocidas (e.g., PEDv, TTSuV). Un alto porcentaje de identidad (>95%) indica un match fuerte con patógenos conocidos.
   - **Decisión**: Si se detectan genes virales patogénicos, prioriza muestras para confirmación experimental. Usa visualizaciones (e.g., `analisis_grafico.py`) para identificar regiones conservadas o variables que indiquen mutaciones.

### 2. **Análisis Filogenético**
   - **Resultados clave**: Árboles generados en `ejercicio_filogenia.py` o `solucion_filogenia.py`, archivos como `aligned_genomes.aln`.
   - **Interpretación**: Observa la posición de tus secuencias en el árbol respecto a clados conocidos. Secuencias cercanas a linajes epidémicos sugieren riesgo alto de transmisión.
   - **Decisión**: Si el aislamiento forma un clado con virus virulentos recientes, recomienda aislamiento de animales afectados y vigilancia epidemiológica. Usa métricas como bootstrap para validar la robustez del árbol.

### 3. **Simulación de PCR In Silico**
   - **Resultados clave**: Outputs de `analisis_pcr_insilico.py`, e.g., `pcr_simulation_results.png`.
   - **Interpretación**: Bandas de PCR positivas indican amplificación de targets virales. Múltiples bandas o tamaños inesperados sugieren coinfecciones o variantes.
   - **Decisión**: Resultados positivos confirman presencia viral; diseña primers específicos para subtipos si hay ambigüedad. Negativos no descartan infección baja; sugiere métodos más sensibles como qPCR.

### 4. **Predicción de Riesgos con Machine Learning**
   - **Resultados clave**: Modelos en `ml_risk_predictor.py`, `pedv_risk_predictor_model.joblib`, gráficos como `pedv_factores_riesgo.png`, `ttsuv2_carga_viral_analisis.png`.
   - **Interpretación**: Predicciones de riesgo (e.g., alta/baja carga viral) basadas en features genómicas. Usa regresión lineal (`linear_regression_prediction.py`) para tendencias cuantitativas.
   - **Decisión**: Riesgo alto (> umbral, e.g., 0.7) implica acciones inmediatas como vacunación o sacrificio. Evalúa factores como edad animal o condiciones de granja para contextualizar. Valida el modelo con datos nuevos para evitar falsos positivos.

### Flujo General de Toma de Decisiones
- **Integración de resultados**: Combina evidencia de todos los análisis. Por ejemplo, un match BLAST + clado epidémico + PCR positivo + riesgo ML alto = alta prioridad para intervención.
- **Umbrales y criterios**: Define umbrales basados en literatura (e.g., identidad >90% para preocupación). Documenta incertidumbres (e.g., baja cobertura de secuencia).
- **Recomendaciones prácticas**: Basado en hallazgos, propone pasos como:
  - Monitoreo continuo si riesgo medio.
  - Reporte a autoridades sanitarias si riesgo alto.
  - Estudios adicionales (e.g., secuenciación completa) para confirmación.
- **Herramientas para visualización**: Usa scripts como `analisis_grafico.py` para gráficos que apoyen decisiones, asegurando reproducibilidad.

Este enfoque fomenta el pensamiento crítico, donde los resultados no son absolutos sino guías para acciones éticas y científicas en salud animal.