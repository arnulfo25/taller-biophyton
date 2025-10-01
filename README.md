**Manual del Curso Taller**
**Título:** Interpretación Biológica y Vigilancia Genómica para la Gestión de Brotes Mediante Biopython

---

### **Introducción General**

Bienvenidos al taller de "Interpretación Biológica y Vigilancia Genómica". En la era de la genómica, la capacidad de analizar grandes volúmenes de datos biológicos es fundamental para responder rápidamente a brotes de enfermedades. Este curso está diseñado para proporcionar a los participantes las habilidades computacionales necesarias para utilizar Biopython, una biblioteca esencial del lenguaje de programación Python, en el contexto de la vigilancia epidemiológica. A través de una serie de ejercicios prácticos, exploraremos cómo simular experimentos de laboratorio, analizar relaciones evolutivas, identificar factores de riesgo y construir modelos predictivos para la gestión eficaz de brotes.

---

### **Capítulo 1: Simulación de PCR in Silico (`analisis_pcr_insilico.py`)**

**Introducción:**
Antes de invertir tiempo y recursos en el laboratorio, es crucial validar que los cebadores (primers) diseñados para una reacción de PCR sean específicos y efectivos. La simulación *in silico* (realizada por computadora) nos permite predecir el resultado de una PCR, verificando si los cebadores se unirán al genoma objetivo y cuál será el tamaño del producto esperado (amplicón). Este capítulo introduce las herramientas de Biopython para automatizar esta validación.

**Objetivos de Aprendizaje:**
Al finalizar este capítulo, el alumno será capaz de:
*   Comprender la importancia de la validación de cebadores *in silico*.
*   Utilizar Biopython para realizar una búsqueda BLAST programática de secuencias de cebadores contra una base de datos de genomas.
*   Predecir el tamaño del amplicón de una reacción de PCR simulada.
*   Generar una visualización para comparar los resultados teóricos con los simulados.

Para el funcionamiento del script analisis_pcr_insilico.py, se deben utilizar las siguientes extensiones de archivo:

.fasta, .fa, .fna: Para los archivos de secuencias de genomas que se encuentran en el directorio genomes.
.png: Para el archivo de imagen de salida que contiene el gráfico de los resultados de la simulación de PCR.
Adicionalmente, el script interactúa con una base de datos BLAST (pig_viruses_db), la cual está compuesta por múltiples archivos con extensiones como .ndb, .nhr, .nin, .nsq, etc. Estos son generados por la herramienta makeblastdb.

---

### **Capítulo 2: Análisis Filogenético (`analisis_grafico.py`)**

**Introducción:**
El análisis filogenético es clave para la vigilancia genómica, ya que permite rastrear el origen y la diseminación de patógenos al revelar sus relaciones evolutivas. En este capítulo, aprenderemos a construir un árbol filogenético a partir de secuencias de genes virales. Este proceso implica alinear las secuencias para hacerlas comparables y luego usar modelos matemáticos para inferir el árbol que mejor represente su historia evolutiva.

**Objetivos de Aprendizaje:**
Al finalizar este capítulo, el alumno será capaz de:
*   Entender los fundamentos del análisis filogenético y su aplicación en la vigilancia de brotes.
*   Realizar un alineamiento múltiple de secuencias utilizando herramientas externas como MUSCLE a través de Biopython.
*   Calcular distancias evolutivas entre secuencias alineadas.
*   Construir y visualizar un árbol filogenético utilizando el método de Neighbor-Joining.

Para el funcionamiento del script analisis_grafico.py, se utilizan las siguientes extensiones de archivo:

.fasta: Para los archivos de entrada que contienen las secuencias genómicas.
.aln: Para los archivos de alineamiento generados por MUSCLE.
.png: Para la imagen de salida que contiene el árbol filogenético.

---

### **Capítulo 3: Identificación de Factores de Riesgo (`analisis_riesgo_viral.py`)**

**Introducción:**
La epidemiología molecular combina datos genómicos con datos epidemiológicos para obtener una comprensión más profunda de la dinámica de las enfermedades. Este capítulo se enfoca en analizar datos de campo para identificar factores de riesgo asociados a enfermedades virales en cerdos (PEDv y TTSuV2). Utilizaremos la potencia de las bibliotecas de análisis de datos de Python para encontrar patrones y correlaciones.

**Objetivos de Aprendizaje:**
Al finalizar este capítulo, el alumno será capaz de:
*   Manejar y analizar datos tabulares epidemiológicos con la biblioteca Pandas.
*   Identificar y cuantificar factores de riesgo asociados a casos positivos (ej. características de camiones para PEDv).
*   Analizar la correlación entre variables cuantitativas, como la carga viral y la edad.
*   Crear visualizaciones efectivas (gráficos de barras, diagramas de dispersión) para comunicar los hallazgos.


Para el funcionamiento del script analisis_riesgo_viral.py, la única extensión de archivo que se utiliza es:

.png: Para los archivos de imagen de salida que se generan (pedv_factores_riesgo.png y ttsuv2_carga_viral_analisis.png).
Este script no lee datos desde ningún archivo externo; la información se encuentra definida directamente en el código.
---

### **Capítulo 4: Modelado Predictivo con Regresión Lineal (`linear_regression_prediction.py`)**

**Introducción:**
El modelado predictivo puede ser una herramienta poderosa para la toma de decisiones en sanidad animal. En este capítulo, construiremos un modelo simple de regresión lineal para estimar la carga viral de TTSuV2 en un animal basándonos en su edad. Este ejercicio servirá como introducción a los conceptos básicos del Machine Learning y su aplicación en la biología.

**Objetivos de Aprendizaje:**
Al finalizar este capítulo, el alumno será capaz de:
*   Comprender los principios básicos de la regresión lineal.
*   Entrenar un modelo predictivo utilizando la biblioteca Scikit-learn.
*   Utilizar el modelo entrenado para hacer predicciones sobre nuevos datos.
*   Visualizar la relación entre las variables y la línea de regresión del modelo.



Para el funcionamiento del script linear_regression_prediction.py, la única extensión de archivo que se utiliza es:

.png: Para el archivo de imagen de salida que se genera (ttsuv2_linear_regression.png).
Este script no lee datos desde ningún archivo externo; la información se encuentra definida directamente en el código.
---

### **Capítulo 5: Clasificación de Riesgo con Machine Learning (`ml_risk_predictor.py`)**

**Introducción:**
Ampliando los conceptos del capítulo anterior, este módulo introduce un modelo de Machine Learning más avanzado: el clasificador Random Forest. Nuestro objetivo será desarrollar una herramienta capaz de predecir si un camión tiene un riesgo 'Positivo' o 'Negativo' de estar contaminado con PEDv, basándose en sus características. Este tipo de modelo es fundamental para crear sistemas de alerta temprana y optimizar los recursos de vigilancia.

**Objetivos de Aprendizaje:**
Al finalizar este capítulo, el alumno será capaz de:
*   Entender la diferencia entre problemas de regresión y clasificación en Machine Learning.
*   Entrenar un modelo clasificador Random Forest con datos de riesgo.
*   Guardar un modelo entrenado para su reutilización.
*   Implementar una función que permita a un usuario final obtener una predicción de riesgo para un nuevo caso.



Para el funcionamiento del script ml_risk_predictor.py, se utiliza la siguiente extensión de archivo:

.joblib: Para guardar y cargar el modelo de machine learning entrenado (pedv_risk_predictor_model.joblib).
Este script no lee datos de un archivo externo para el entrenamiento, ya que se encuentran definidos directamente en el código.
---

### **Bibliografía y Recursos Recomendados**

*   **Biopython Tutorial and Cookbook:** Cock, P. A., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., ... & de Hoon, M. J. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, 25(11), 1422-1423. (Disponible en línea: http://biopython.org/DIST/docs/tutorial/Tutorial.html)
*   **Documentación oficial de Pandas:** McKinney, W. (2010). Data structures for statistical computing in python. *Proceedings of the 9th Python in Science Conference*, 445, 51-56. (Disponible en línea: https://pandas.pydata.org/docs/)
*   **Documentación oficial de Scikit-learn:** Pedregosa, F., Varoquaux, G., Gramfort, A., Michel, V., Thirion, B., Grisel, O., ... & Duchesnay, E. (2011). Scikit-learn: Machine learning in Python. *Journal of machine learning research*, 12(Oct), 2825-2830. (Disponible en línea: https://scikit-learn.org/stable/documentation.html)
*   **Bioinformatics with Python Cookbook:** Tiago Antao (2018). *Bioinformatics with Python Cookbook - Second Edition*. Packt Publishing.
