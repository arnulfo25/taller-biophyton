# Taller Biophyton - Estructura de Ejercicios

## analisis_pcr_insilico.py
- **Objetivo**: Simular y validar virtualmente (in silico) la efectividad y especificidad de cebadores de PCR.
  - Verifica la unión de los cebadores contra una base de datos de genomas (BLAST).
  - Simula la reacción de PCR para predecir el tamaño de los productos (amplicones).
  - Genera un gráfico para comparar los resultados esperados vs. los simulados.

## analisis_grafico.py
- **Objetivo**: Construir y visualizar un árbol filogenético para entender las relaciones evolutivas entre secuencias de genes virales.
  - Combina secuencias de diferentes archivos.
  - Realiza un alineamiento múltiple de secuencias (con MUSCLE).
  - Calcula distancias evolutivas y construye el árbol (método Neighbor-Joining).
  - Guarda el árbol como una imagen.

## analisis_riesgo_viral.py
- **Objetivo**: Analizar datos epidemiológicos para identificar y visualizar factores de riesgo asociados a enfermedades virales en cerdos.
  - **PEDv**: Analiza características de camiones (tipo de planta, etc.) para encontrar cuáles son más comunes en casos positivos.
  - **TTSuV2**: Relaciona la carga viral con el estado de salud y la edad de los animales.
  - Genera gráficos de barras y dispersión para presentar los hallazgos.

## linear_regression_prediction.py
- **Objetivo**: Crear un modelo predictivo simple (regresión lineal) para estimar la carga viral de TTSuV2 basándose en la edad del animal.
  - Entrena el modelo con datos de edad y carga viral.
  - Permite al usuario ingresar una edad y obtener una predicción de la carga viral.
  - Visualiza la relación lineal con un gráfico de dispersión y la línea de regresión.

## ml_risk_predictor.py
- **Objetivo**: Desarrollar un modelo de Machine Learning (Random Forest) para clasificar el riesgo de que un camión sea positivo a PEDv.
  - Entrena un modelo con datos históricos de camiones y sus factores de riesgo.
  - Guarda el modelo entrenado para su uso futuro.
  - Permite al usuario ingresar las características de un nuevo camión para predecir si el riesgo es 'Positivo' o 'Negativo'.