import pandas as pd
from scipy.stats import fisher_exact

# --- Sección 0: Notas sobre Biopython y Epidemiología ---
print("---" + "Notas para el Estudiante de Epidemiología" + "---")
print("Biopython es una librería de Python fundamental en bioinformática, utilizada para trabajar con secuencias de ADN/ARN/Proteínas, estructuras 3D, bases de datos biológicas, etc.")
print("Para este ejercicio de epidemiología, las funcionalidades de Biopython no son directamente aplicables, ya que nos enfocamos en el análisis de datos tabulares (factores de riesgo).")
print("Sin embargo, Python es una herramienta muy potente y versátil en epidemiología para manejar, limpiar y analizar datos, así como para realizar cálculos estadísticos. Aquí usaremos librerías como Pandas para la manipulación de datos y SciPy para estadísticas.")
print("---------------------------------------------------" + "\n")

# --- Sección 1: Creación de la Base de Datos Simulada ---
# Los datos son SIMULADOS y se basan en las proporciones y asociaciones descritas en las fuentes.
# Este es un ejemplo simplificado para el ejercicio.
# En un estudio real, estos datos provendrían de encuestas y pruebas de laboratorio.
datos_simulados = [
    # Camiones con PEDv (mayoría, dada la prevalencia del 71.8% al ingreso)
    {'id': 1, 'tipo_planta': 'Nacional', 'zona_sacrificio': 'Mayor', 'uso_vehiculo': 'Cerdos_Otros_Productos', 'visita_concentrados': 'Si', 'desciende_vehiculo_operador': 'Si', 'limpieza_vehiculo': 'No', 'uso_desinfectante': 'No', 'pedv_presente': 'Si'},
    {'id': 2, 'tipo_planta': 'Nacional-Exportación', 'zona_sacrificio': 'Mayor', 'uso_vehiculo': 'Cerdos_Otros_Productos', 'visita_concentrados': 'Si', 'desciende_vehiculo_operador': 'Si', 'limpieza_vehiculo': 'No', 'uso_desinfectante': 'No', 'pedv_presente': 'Si'},
    {'id': 3, 'tipo_planta': 'Nacional', 'zona_sacrificio': 'Mayor', 'uso_vehiculo': 'Cerdos_Otros_Productos', 'visita_concentrados': 'Si', 'desciende_vehiculo_operador': 'Si', 'limpieza_vehiculo': 'No', 'uso_desinfectante': 'No', 'pedv_presente': 'Si'},
    {'id': 4, 'tipo_planta': 'Nacional', 'zona_sacrificio': 'Mayor', 'uso_vehiculo': 'Cerdos_Otros_Productos', 'visita_concentrados': 'Si', 'desciende_vehiculo_operador': 'Si', 'limpieza_vehiculo': 'No', 'uso_desinfectante': 'No', 'pedv_presente': 'Si'},
    {'id': 5, 'tipo_planta': 'Nacional-Exportación', 'zona_sacrificio': 'Mayor', 'uso_vehiculo': 'Cerdos_Otros_Productos', 'visita_concentrados': 'Si', 'desciende_vehiculo_operador': 'Si', 'limpieza_vehiculo': 'No', 'uso_desinfectante': 'No', 'pedv_presente': 'Si'},
    {'id': 6, 'tipo_planta': 'Nacional', 'zona_sacrificio': 'Mayor', 'uso_vehiculo': 'Cerdos_Otros_Productos', 'visita_concentrados': 'Si', 'desciende_vehiculo_operador': 'Si', 'limpieza_vehiculo': 'No', 'uso_desinfectante': 'No', 'pedv_presente': 'Si'},
    {'id': 7, 'tipo_planta': 'Nacional', 'zona_sacrificio': 'Mayor', 'uso_vehiculo': 'Exclusivo_Cerdos', 'visita_concentrados': 'No', 'desciende_vehiculo_operador': 'Si', 'limpieza_vehiculo': 'Si', 'uso_desinfectante': 'No', 'pedv_presente': 'Si'}, # PEDv Si a pesar de algunos buenos factores
    {'id': 8, 'tipo_planta': 'Nacional-Exportación', 'zona_sacrificio': 'Mayor', 'uso_vehiculo': 'Exclusivo_Cerdos', 'visita_concentrados': 'No', 'desciende_vehiculo_operador': 'Si', 'limpieza_vehiculo': 'Si', 'uso_desinfectante': 'No', 'pedv_presente': 'Si'}, # PEDv Si a pesar de algunos buenos factores
    {'id': 9, 'tipo_planta': 'Local', 'zona_sacrificio': 'Menor', 'uso_vehiculo': 'Exclusivo_Cerdos', 'visita_concentrados': 'No', 'desciende_vehiculo_operador': 'No', 'limpieza_vehiculo': 'Si', 'uso_desinfectante': 'Si', 'pedv_presente': 'Si'}, # PEDv Si a pesar de buenos factores
    {'id': 10, 'tipo_planta': 'Nacional', 'zona_sacrificio': 'Mayor', 'uso_vehiculo': 'Cerdos_Otros_Productos', 'visita_concentrados': 'Si', 'desciende_vehiculo_operador': 'Si', 'limpieza_vehiculo': 'No', 'uso_desinfectante': 'No', 'pedv_presente': 'Si'},
    {'id': 11, 'tipo_planta': 'Nacional', 'zona_sacrificio': 'Mayor', 'uso_vehiculo': 'Cerdos_Otros_Productos', 'visita_concentrados': 'Si', 'desciende_vehiculo_operador': 'Si', 'limpieza_vehiculo': 'No', 'uso_desinfectante': 'No', 'pedv_presente': 'Si'},
    {'id': 12, 'tipo_planta': 'Nacional', 'zona_sacrificio': 'Mayor', 'uso_vehiculo': 'Cerdos_Otros_Productos', 'visita_concentrados': 'Si', 'desciende_vehiculo_operador': 'Si', 'limpieza_vehiculo': 'No', 'uso_desinfectante': 'No', 'pedv_presente': 'Si'},
    {'id': 13, 'tipo_planta': 'Nacional-Exportación', 'zona_sacrificio': 'Mayor', 'uso_vehiculo': 'Cerdos_Otros_Productos', 'visita_concentrados': 'Si', 'desciende_vehiculo_operador': 'Si', 'limpieza_vehiculo': 'No', 'uso_desinfectante': 'No', 'pedv_presente': 'Si'},
    {'id': 14, 'tipo_planta': 'Nacional', 'zona_sacrificio': 'Mayor', 'uso_vehiculo': 'Cerdos_Otros_Productos', 'visita_concentrados': 'Si', 'desciende_vehiculo_operador': 'Si', 'limpieza_vehiculo': 'No', 'uso_desinfectante': 'No', 'pedv_presente': 'Si'},
    {'id': 15, 'tipo_planta': 'Nacional', 'zona_sacrificio': 'Mayor', 'uso_vehiculo': 'Cerdos_Otros_Productos', 'visita_concentrados': 'Si', 'desciende_vehiculo_operador': 'Si', 'limpieza_vehiculo': 'No', 'uso_desinfectante': 'No', 'pedv_presente': 'Si'},
    {'id': 16, 'tipo_planta': 'Nacional-Exportación', 'zona_sacrificio': 'Mayor', 'uso_vehiculo': 'Cerdos_Otros_Productos', 'visita_concentrados': 'Si', 'desciende_vehiculo_operador': 'Si', 'limpieza_vehiculo': 'No', 'uso_desinfectante': 'No', 'pedv_presente': 'Si'},
    {'id': 17, 'tipo_planta': 'Nacional', 'zona_sacrificio': 'Mayor', 'uso_vehiculo': 'Cerdos_Otros_Productos', 'visita_concentrados': 'Si', 'desciende_vehiculo_operador': 'Si', 'limpieza_vehiculo': 'No', 'uso_desinfectante': 'No', 'pedv_presente': 'Si'},
    {'id': 18, 'tipo_planta': 'Local', 'zona_sacrificio': 'Menor', 'uso_vehiculo': 'Exclusivo_Cerdos', 'visita_concentrados': 'No', 'desciende_vehiculo_operador': 'No', 'limpieza_vehiculo': 'Si', 'uso_desinfectante': 'Si', 'pedv_presente': 'Si'}, # PEDv Si a pesar de buenos factores
    {'id': 19, 'tipo_planta': 'Nacional', 'zona_sacrificio': 'Mayor', 'uso_vehiculo': 'Cerdos_Otros_Productos', 'visita_concentrados': 'Si', 'desciende_vehiculo_operador': 'Si', 'limpieza_vehiculo': 'No', 'uso_desinfectante': 'No', 'pedv_presente': 'Si'},
    {'id': 20, 'tipo_planta': 'Nacional-Exportación', 'zona_sacrificio': 'Mayor', 'uso_vehiculo': 'Cerdos_Otros_Productos', 'visita_concentrados': 'Si', 'desciende_vehiculo_operador': 'Si', 'limpieza_vehiculo': 'No', 'uso_desinfectante': 'No', 'pedv_presente': 'Si'},

    # Camiones sin PEDv (minoría)
    {'id': 21, 'tipo_planta': 'Local', 'zona_sacrificio': 'Menor', 'uso_vehiculo': 'Exclusivo_Cerdos', 'visita_concentrados': 'No', 'desciende_vehiculo_operador': 'No', 'limpieza_vehiculo': 'Si', 'uso_desinfectante': 'Si', 'pedv_presente': 'No'},
    {'id': 22, 'tipo_planta': 'Local', 'zona_sacrificio': 'Menor', 'uso_vehiculo': 'Exclusivo_Cerdos', 'visita_concentrados': 'No', 'desciende_vehiculo_operador': 'No', 'limpieza_vehiculo': 'Si', 'uso_desinfectante': 'Si', 'pedv_presente': 'No'},
    {'id': 23, 'tipo_planta': 'Local', 'zona_sacrificio': 'Menor', 'uso_vehiculo': 'Exclusivo_Cerdos', 'visita_concentrados': 'No', 'desciende_vehiculo_operador': 'No', 'limpieza_vehiculo': 'Si', 'uso_desinfectante': 'Si', 'pedv_presente': 'No'},
    {'id': 24, 'tipo_planta': 'Nacional', 'zona_sacrificio': 'Menor', 'uso_vehiculo': 'Exclusivo_Cerdos', 'visita_concentrados': 'No', 'desciende_vehiculo_operador': 'No', 'limpieza_vehiculo': 'Si', 'uso_desinfectante': 'Si', 'pedv_presente': 'No'}, # Nacional pero en zona menor sacrificio y buenas prácticas
    {'id': 25, 'tipo_planta': 'Local', 'zona_sacrificio': 'Menor', 'uso_vehiculo': 'Exclusivo_Cerdos', 'visita_concentrados': 'No', 'desciende_vehiculo_operador': 'No', 'limpieza_vehiculo': 'Si', 'uso_desinfectante': 'Si', 'pedv_presente': 'No'},
    {'id': 26, 'tipo_planta': 'Local', 'zona_sacrificio': 'Menor', 'uso_vehiculo': 'Exclusivo_Cerdos', 'visita_concentrados': 'No', 'desciende_vehiculo_operador': 'No', 'limpieza_vehiculo': 'Si', 'uso_desinfectante': 'Si', 'pedv_presente': 'No'},
    {'id': 27, 'tipo_planta': 'Nacional-Exportación', 'zona_sacrificio': 'Menor', 'uso_vehiculo': 'Exclusivo_Cerdos', 'visita_concentrados': 'No', 'desciende_vehiculo_operador': 'No', 'limpieza_vehiculo': 'Si', 'uso_desinfectante': 'Si', 'pedv_presente': 'No'}, # Nac-Exp pero en zona menor y buenas prácticas
    {'id': 28, 'tipo_planta': 'Local', 'zona_sacrificio': 'Menor', 'uso_vehiculo': 'Cerdos_Otros_Productos', 'visita_concentrados': 'No', 'desciende_vehiculo_operador': 'No', 'limpieza_vehiculo': 'Si', 'uso_desinfectante': 'Si', 'pedv_presente': 'No'}, # Cerdos y otros productos pero buenas prácticas y local
    {'id': 29, 'tipo_planta': 'Local', 'zona_sacrificio': 'Menor', 'uso_vehiculo': 'Exclusivo_Cerdos', 'visita_concentrados': 'No', 'desciende_vehiculo_operador': 'No', 'limpieza_vehiculo': 'Si', 'uso_desinfectante': 'Si', 'pedv_presente': 'No'},
    {'id': 30, 'tipo_planta': 'Local', 'zona_sacrificio': 'Menor', 'uso_vehiculo': 'Exclusivo_Cerdos', 'visita_concentrados': 'No', 'desciende_vehiculo_operador': 'No', 'limpieza_vehiculo': 'Si', 'uso_desinfectante': 'Si', 'pedv_presente': 'No'},
]

df = pd.DataFrame(datos_simulados)

# --- Sección 2: Cálculo de Prevalencia General ---
total_camiones = len(df)
camiones_pedv_positivos = df[df['pedv_presente'] == 'Si'].shape[0]
prevalencia_general = (camiones_pedv_positivos / total_camiones) * 100

print(f"---" + "Análisis de Prevalencia General" + "---")
print(f"Número total de camiones en la muestra: {total_camiones}")
print(f"Camiones con PEDv positivo: {camiones_pedv_positivos}")
print(f"**Prevalencia general de PEDv: {prevalencia_general:.2f}%**")
print(f"*(Según el estudio real, la prevalencia al ingreso fue del 71.8%)")
print("-------------------------------------------" + "\n")

# --- Sección 3: Análisis Bivariado (Prevalencia por Factor) ---
print("---" + "Prevalencia de PEDv por Categoría de Factor" + "---")
for factor in df.columns.drop(['id', 'pedv_presente']):
    print(f"\nFactor: {factor}")
    grupo_por_factor = df.groupby(factor)['pedv_presente'].value_counts(normalize=True).unstack(fill_value=0)
    if 'Si' in grupo_por_factor.columns:
        prevalencias_factor = grupo_por_factor['Si'] * 100
        print(prevalencias_factor.round(2))
    else:
        print("No se detectó PEDv en ninguna categoría de este factor en la muestra simulada.")
print("-------------------------------------------" + "\n")

# --- Sección 4: Cálculo de Odds Ratio (OR) ---
# Función para calcular el Odds Ratio (OR) y su intervalo de confianza (aproximado)
def calcular_odds_ratio(df, factor_riesgo, categoria_expuesta, categoria_no_expuesta):
    # Crear tabla de contingencia 2x2
    #                 PEDv_Si   PEDv_No
    # Expuesto (a)      a         b
    # No Expuesto (c)   c         d
    
    a = df[(df[factor_riesgo] == categoria_expuesta) & (df['pedv_presente'] == 'Si')].shape[0]
    b = df[(df[factor_riesgo] == categoria_expuesta) & (df['pedv_presente'] == 'No')].shape[0]
    c = df[(df[factor_riesgo] == categoria_no_expuesta) & (df['pedv_presente'] == 'Si')].shape[0]
    d = df[(df[factor_riesgo] == categoria_no_expuesta) & (df['pedv_presente'] == 'No')].shape[0]

    table = [[a, b], [c, d]]
    
    if a == 0 or b == 0 or c == 0 or d == 0:
        # Añadir 0.5 a todas las celdas si hay ceros para evitar divisiones por cero en el OR y log(0) en el IC
        # Esta es una corrección de Woolf, común en epidemiología para tablas con ceros.
        a += 0.5
        b += 0.5
        c += 0.5
        d += 0.5

    odds_expuestos = (a / b) if b != 0 else float('inf')
    odds_no_expuestos = (c / d) if d != 0 else float('inf')
    
    or_val = odds_expuestos / odds_no_expuestos

    # Cálculo del Intervalo de Confianza (95%) aproximado para el OR
    se_log_or = (1/a + 1/b + 1/c + 1/d)**0.5
    lower_ci = or_val * (2.71828 ** (-1.96 * se_log_or))
    upper_ci = or_val * (2.71828 ** (1.96 * se_log_or))
    
    # También se puede usar fisher_exact para obtener el p-value (aunque para OR el IC ya es buena indicación)
    _, p_value = fisher_exact(table)

    return or_val, lower_ci, upper_ci, p_value, table

print("---" + "Cálculo de Odds Ratios (OR)" + "---")

# Ejemplo 1: Tipo de planta (Nacional vs Local) [OR real: 15.95 para Nacional, 9.02 para Nacional-Exportación vs Local]
or_tipo_planta, lc_tipo_planta, uc_tipo_planta, pval_tipo_planta, tabla_tipo_planta = calcular_odds_ratio(df, 'tipo_planta', 'Nacional', 'Local')
print(f"\n**Factor: Tipo de Planta (Nacional vs Local)**")
print(f"  Tabla de contingencia 2x2 (PEDv Si/No vs Nacional/Local):")
print(f"    Nacional: {tabla_tipo_planta[0][0]} Si, {tabla_tipo_planta[0][1]} No")
print(f"    Local:    {tabla_tipo_planta[1][0]} Si, {tabla_tipo_planta[1][1]} No")
print(f"  Odds Ratio (OR): {or_tipo_planta:.2f} (IC 95%: {lc_tipo_planta:.2f} - {uc_tipo_planta:.2f})")
print(f"  P-valor (Fisher's Exact): {pval_tipo_planta:.3f}")
print(f"  *Interpretación: Los camiones que visitan plantas de tipo 'Nacional' tienen {or_tipo_planta:.2f} veces más probabilidades de ser positivos a PEDv que los que visitan plantas 'Local'. Esto concuerda con los hallazgos del estudio.")


# Ejemplo 2: Visita a plantas de concentrados (Si vs No) [OR real: 13.56]
or_visita_concentrados, lc_visita_concentrados, uc_visita_concentrados, pval_visita_concentrados, tabla_visita_concentrados = calcular_odds_ratio(df, 'visita_concentrados', 'Si', 'No')
print(f"\n**Factor: Visita a Plantas de Concentrados (Si vs No)**")
print(f"  Tabla de contingencia 2x2 (PEDv Si/No vs Visita concentrados Si/No):")
print(f"    Visita Si: {tabla_visita_concentrados[0][0]} Si, {tabla_visita_concentrados[0][1]} No")
print(f"    Visita No: {tabla_visita_concentrados[1][0]} Si, {tabla_visita_concentrados[1][1]} No")
print(f"  Odds Ratio (OR): {or_visita_concentrados:.2f} (IC 95%: {lc_visita_concentrados:.2f} - {uc_visita_concentrados:.2f})")
print(f"  P-valor (Fisher's Exact): {pval_visita_concentrados:.3f}")
print(f"  *Interpretación: Los camiones que visitan plantas de concentrados tienen {or_visita_concentrados:.2f} veces más probabilidades de ser positivos a PEDv. Este es un factor clave identificado en el estudio.")

# Ejemplo 3: Uso exclusivo del vehículo (Cerdos_Otros_Productos vs Exclusivo_Cerdos) [OR real: 3.75]
or_uso_vehiculo, lc_uso_vehiculo, uc_uso_vehiculo, pval_uso_vehiculo, tabla_uso_vehiculo = calcular_odds_ratio(df, 'uso_vehiculo', 'Cerdos_Otros_Productos', 'Exclusivo_Cerdos')
print(f"\n**Factor: Uso del Vehículo (Cerdos y Otros Productos vs Exclusivo para Cerdos)**")
print(f"  Tabla de contingencia 2x2 (PEDv Si/No vs Uso mixto/exclusivo):")
print(f"    Uso mixto:  {tabla_uso_vehiculo[0][0]} Si, {tabla_uso_vehiculo[0][1]} No")
print(f"    Exclusivo: {tabla_uso_vehiculo[1][0]} Si, {tabla_uso_vehiculo[1][1]} No")
print(f"  Odds Ratio (OR): {or_uso_vehiculo:.2f} (IC 95%: {lc_uso_vehiculo:.2f} - {uc_uso_vehiculo:.2f})")
print(f"  P-valor (Fisher's Exact): {pval_uso_vehiculo:.3f}")
print(f"  *Interpretación: Los camiones que transportan cerdos y otros productos tienen {or_uso_vehiculo:.2f} veces más probabilidades de ser positivos a PEDv que los de uso exclusivo. Esto se identificó como un factor de riesgo importante.")


# Ejemplo 4: Uso de desinfectante en la limpieza (No vs Si) [Factor protector, OR real < 1]
# Para un factor protector, exponemos "No uso" y comparamos con "Uso Si" para obtener un OR > 1.
or_uso_desinfectante, lc_uso_desinfectante, uc_uso_desinfectante, pval_uso_desinfectante, tabla_uso_desinfectante = calcular_odds_ratio(df, 'uso_desinfectante', 'No', 'Si')
print(f"\n**Factor: Uso de Desinfectante en la Limpieza (No uso vs Si uso)**")
print(f"  Tabla de contingencia 2x2 (PEDv Si/No vs No uso desinfectante/Si uso desinfectante):")
print(f"    No uso: {tabla_uso_desinfectante[0][0]} Si, {tabla_uso_desinfectante[0][1]} No")
print(f"    Si uso: {tabla_uso_desinfectante[1][0]} Si, {tabla_uso_desinfectante[1][1]} No")
print(f"  Odds Ratio (OR): {or_uso_desinfectante:.2f} (IC 95%: {lc_uso_desinfectante:.2f} - {uc_uso_desinfectante:.2f})")
print(f"  P-valor (Fisher's Exact): {pval_uso_desinfectante:.3f}")
print(f"  *Interpretación: Los camiones que NO usan desinfectante tienen {or_uso_desinfectante:.2f} veces más probabilidades de ser positivos a PEDv que los que SÍ lo usan. El uso de desinfectantes es un factor protector clave.*")

print("-------------------------------------------" + "\n")

# --- Sección 5: Conclusiones y Recomendaciones Basadas en los Hallazgos ---
print("---" + "Conclusiones y Medidas de Bioseguridad Sugeridas (basado en el estudio y ORs simulados)" + "---")
print("Este análisis simulado (y el estudio real) resalta que los camiones de transporte de cerdos son un vector crucial para la diseminación del PEDv en Colombia.")
print("\n**Factores de Riesgo Clave Confirmados (incrementan la probabilidad de PEDv):**")
print("1.  **Tipo de Planta:** Las plantas de tipo Nacional y Nacional-Exportación, posiblemente por su mayor volumen de sacrificio y flujo de vehículos, presentan mayor riesgo.")
print("2.  **Visita a Plantas de Concentrados:** La circulación de camiones entre granjas y plantas de concentrados actúa como un punto de contaminación, facilitando la propagación del virus.")
print("3.  **Uso No Exclusivo del Vehículo:** Los camiones que transportan cerdos y otros productos, o se usan para múltiples granjas, tienen un riesgo significativamente mayor de contaminación cruzada.")
print("4.  **Descenso del Personal del Vehículo:** El contacto del personal con el entorno de la planta/granja y luego con el camión es un riesgo significativo.")

print("\n**Factores Protectores Clave Confirmados (disminuyen la probabilidad de PEDv):**")
print("1.  **Limpieza y Desinfección Rigurosa del Vehículo:** El lavado diario o semanal y, crucialmente, el **uso adecuado de desinfectantes** son esenciales para inactivar el virus en las superficies del camión. La limpieza profunda de carrocería y cabina también es vital.")
print("2.  **Uso Exclusivo del Transporte Porcino:** Dedicar un vehículo únicamente al transporte de cerdos, idealmente para una sola granja, reduce drásticamente el riesgo de introducción y diseminación del virus.")
print("3.  **Medidas de Bioseguridad en Granjas:** El uso de arcos de desinfección, bombas de espalda con desinfectantes y agua a presión en las granjas antes del ingreso/egreso de camiones contribuye a la protección.")
print("4.  **Higiene del Personal:** Evitar que el personal descienda del vehículo en áreas de riesgo, así como el uso de dotación y elementos de protección proporcionados por la granja, son medidas importantes.")

print("\n**Recomendaciones para el Control del PEDv en el Transporte (basadas en los hallazgos):**")
print("Para controlar eficazmente la diseminación del PEDv, es imperativo:")
print("1.  **Establecer y hacer cumplir protocolos de lavado y desinfección estandarizados y rigurosos** en todas las plantas de beneficio y granjas, asegurando el uso de productos efectivos y tiempos de contacto adecuados.")
print("2.  **Promover el uso exclusivo de vehículos para el transporte de cerdos** y evitar la movilización de otros productos o la visita a plantas de concentrados con los mismos vehículos.")
print("3.  **Mejorar la bioseguridad en las plantas de beneficio y los puntos de carga/descarga**, ya que son cruciales para la contaminación de vehículos.")
print("4.  **Capacitar y concienciar a todo el personal** (transportistas, operarios de granja y planta) sobre la importancia de las medidas de bioseguridad, incluyendo la higiene personal y el uso de equipo de protección.")
print("5.  **Considerar la infraestructura de los camiones:** Las superficies lisas y lavables (metal) son preferibles a la madera para una desinfección efectiva.")

print("\nEste ejercicio subraya cómo el análisis de datos epidemiológicos, incluso simulados, puede guiar decisiones prácticas en salud animal y bioseguridad.")
