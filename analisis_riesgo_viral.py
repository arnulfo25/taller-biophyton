# --- Tarea 1: Análisis de Factores de Riesgo para PEDv en Camiones ---

# Datos simulados de camiones (PEDv)
# Cada diccionario representa un camión con sus características y el resultado de la prueba de PEDv.
# 'pcr_pedv_ingreso': 'Positivo' o 'Negativo' (resultado de RT-PCR al ingresar a la planta)
# 'tipo_planta': 'Nacional', 'Nacional-Exportacion', 'Local'
# 'zona_sacrificio': 'Mayor', 'Menor'
# 'uso_exclusivo_vehiculo': 'Si' (solo cerdos) o 'No' (cerdos y otros productos)
# 'visita_plantas_concentrado': 'Si' o 'No'
# 'conductor_desciende': 'Si' o 'No'

datos_camiones = [
    {'id_camion': 1, 'pcr_pedv_ingreso': 'Positivo', 'tipo_planta': 'Nacional', 'zona_sacrificio': 'Mayor', 'uso_exclusivo_vehiculo': 'No', 'visita_plantas_concentrado': 'Si', 'conductor_desciende': 'Si'},
    {'id_camion': 2, 'pcr_pedv_ingreso': 'Negativo', 'tipo_planta': 'Local', 'zona_sacrificio': 'Menor', 'uso_exclusivo_vehiculo': 'Si', 'visita_plantas_concentrado': 'No', 'conductor_desciende': 'No'},
    {'id_camion': 3, 'pcr_pedv_ingreso': 'Positivo', 'tipo_planta': 'Nacional-Exportacion', 'zona_sacrificio': 'Mayor', 'uso_exclusivo_vehiculo': 'No', 'visita_plantas_concentrado': 'Si', 'conductor_desciende': 'Si'},
    {'id_camion': 4, 'pcr_pedv_ingreso': 'Negativo', 'tipo_planta': 'Nacional', 'zona_sacrificio': 'Mayor', 'uso_exclusivo_vehiculo': 'Si', 'visita_plantas_concentrado': 'No', 'conductor_desciende': 'Si'}, # Conductor desciende
    {'id_camion': 5, 'pcr_pedv_ingreso': 'Positivo', 'tipo_planta': 'Nacional', 'zona_sacrificio': 'Mayor', 'uso_exclusivo_vehiculo': 'No', 'visita_plantas_concentrado': 'No', 'conductor_desciende': 'Si'},
    {'id_camion': 6, 'pcr_pedv_ingreso': 'Negativo', 'tipo_planta': 'Local', 'zona_sacrificio': 'Menor', 'uso_exclusivo_vehiculo': 'Si', 'visita_plantas_concentrado': 'No', 'conductor_desciende': 'No'},
    {'id_camion': 7, 'pcr_pedv_ingreso': 'Positivo', 'tipo_planta': 'Nacional-Exportacion', 'zona_sacrificio': 'Mayor', 'uso_exclusivo_vehiculo': 'No', 'visita_plantas_concentrado': 'Si', 'conductor_desciende': 'Si'},
    {'id_camion': 8, 'pcr_pedv_ingreso': 'Positivo', 'tipo_planta': 'Nacional', 'zona_sacrificio': 'Mayor', 'uso_exclusivo_vehiculo': 'No', 'visita_plantas_concentrado': 'Si', 'conductor_desciende': 'Si'},
    {'id_camion': 9, 'pcr_pedv_ingreso': 'Negativo', 'tipo_planta': 'Local', 'zona_sacrificio': 'Menor', 'uso_exclusivo_vehiculo': 'Si', 'visita_plantas_concentrado': 'No', 'conductor_desciende': 'No'},
    {'id_camion': 10, 'pcr_pedv_ingreso': 'Positivo', 'tipo_planta': 'Nacional', 'zona_sacrificio': 'Mayor', 'uso_exclusivo_vehiculo': 'No', 'visita_plantas_concentrado': 'Si', 'conductor_desciende': 'Si'},
    {'id_camion': 11, 'pcr_pedv_ingreso': 'Positivo', 'tipo_planta': 'Nacional', 'zona_sacrificio': 'Mayor', 'uso_exclusivo_vehiculo': 'Si', 'visita_plantas_concentrado': 'Si', 'conductor_desciende': 'Si'}, # Uso exclusivo SI, pero positivo, para mostrar variabilidad
    {'id_camion': 12, 'pcr_pedv_ingreso': 'Negativo', 'tipo_planta': 'Nacional-Exportacion', 'zona_sacrificio': 'Mayor', 'uso_exclusivo_vehiculo': 'No', 'conductor_desciende': 'No'},
]

# Datos simulados de animales (TTSuV2)
# 'id_animal': Identificador del animal
# 'estado_salud': 'PCV2-SD Afectado' o 'Sano'
# 'edad_semanas': Edad del animal en semanas (para ver infección temprana)
# 'carga_viral_ttsuv2_bm': Carga viral de TTSuV2 en médula ósea (log10 copias/mg), un valor más alto indica mayor carga
# Valores de carga viral altos > 7.5 log10 copias/mg [23, 26]
datos_animales = [
    {'id_animal': 101, 'estado_salud': 'PCV2-SD Afectado', 'edad_semanas': 3, 'carga_viral_ttsuv2_bm': 8.2}, # Infección temprana, alta carga
    {'id_animal': 102, 'estado_salud': 'Sano', 'edad_semanas': 3, 'carga_viral_ttsuv2_bm': 6.5},
    {'id_animal': 103, 'estado_salud': 'PCV2-SD Afectado', 'edad_semanas': 7, 'carga_viral_ttsuv2_bm': 8.5}, # Alta carga
    {'id_animal': 104, 'estado_salud': 'Sano', 'edad_semanas': 11, 'carga_viral_ttsuv2_bm': 7.0},
    {'id_animal': 105, 'estado_salud': 'PCV2-SD Afectado', 'edad_semanas': 1, 'carga_viral_ttsuv2_bm': 7.9}, # Infección temprana, alta carga
    {'id_animal': 106, 'estado_salud': 'Sano', 'edad_semanas': 8, 'carga_viral_ttsuv2_bm': 7.2},
    {'id_animal': 107, 'estado_salud': 'PCV2-SD Afectado', 'edad_semanas': 15, 'carga_viral_ttsuv2_bm': 9.0}, # Alta carga
    {'id_animal': 108, 'estado_salud': 'Sano', 'edad_semanas': 12, 'carga_viral_ttsuv2_bm': 6.8},
    {'id_animal': 109, 'estado_salud': 'PCV2-SD Afectado', 'edad_semanas': 5, 'carga_viral_ttsuv2_bm': 8.1}, # Alta carga
    {'id_animal': 110, 'estado_salud': 'Sano', 'edad_semanas': 10, 'carga_viral_ttsuv2_bm': 6.9},
]

print("--- Análisis de Factores de Riesgo para PEDv en Camiones ---")

total_camiones = len(datos_camiones)
print(f"Número total de camiones en el estudio: {total_camiones}")

# 1. Contar camiones positivos a PEDv
camiones_pedv_positivos = [c for c in datos_camiones if c['pcr_pedv_ingreso'] == 'Positivo']
num_pedv_positivos = len(camiones_pedv_positivos)
print(f"Número de camiones positivos a PEDv al ingreso: {num_pedv_positivos}")
if total_camiones > 0:
    porcentaje_pedv_positivos = (num_pedv_positivos / total_camiones) * 100
    print(f"Porcentaje de camiones positivos a PEDv: {porcentaje_pedv_positivos:.2f}%")
else:
    porcentaje_pedv_positivos = 0
    print("No hay camiones para analizar.")

if num_pedv_positivos > 0:
    print("\nFactores de riesgo en camiones positivos a PEDv:")
    # 2. Contar factores de riesgo específicos en camiones positivos
    factores_de_riesgo_pedv = {
        'Tipo de Planta (Nacional/Nacional-Exportacion)': 0,
        'Zona de Sacrificio (Mayor)': 0,
        'Uso Exclusivo del Vehículo (No)': 0,
        'Visita Plantas de Concentrado (Si)': 0,
        'Conductor desciende del vehículo (Si)': 0
    }

    for camion in camiones_pedv_positivos:
        if camion['tipo_planta'] in ['Nacional', 'Nacional-Exportacion']:
            factores_de_riesgo_pedv['Tipo de Planta (Nacional/Nacional-Exportacion)'] += 1
        if camion['zona_sacrificio'] == 'Mayor':
            factores_de_riesgo_pedv['Zona de Sacrificio (Mayor)'] += 1
        if camion['uso_exclusivo_vehiculo'] == 'No':
            factores_de_riesgo_pedv['Uso Exclusivo del Vehículo (No)'] += 1
        if camion['visita_plantas_concentrado'] == 'Si':
            factores_de_riesgo_pedv['Visita Plantas de Concentrado (Si)'] += 1
        if camion['conductor_desciende'] == 'Si':
            factores_de_riesgo_pedv['Conductor desciende del vehículo (Si)'] += 1

    # 3. Calcular porcentajes de estos factores
    for factor, cuenta in factores_de_riesgo_pedv.items():
        porcentaje_factor = (cuenta / num_pedv_positivos) * 100
        print(f"- {factor}: {cuenta} camiones ({porcentaje_factor:.2f}%)")
else:
    print("\nNo hay camiones positivos a PEDv para analizar factores de riesgo.")

# --- Tarea 2: Identificación de Condiciones de TTSuV2 de Alta Carga Viral ---

print("\n--- Identificación de Condiciones de TTSuV2 de Alta Carga Viral ---")

umbral_alta_carga_ttsuv2 = 7.5 # Basado en la discusión sobre cargas virales altas [23, 26]

animales_alta_carga_ttsuv2 = [
    a for a in datos_animales if a['carga_viral_ttsuv2_bm'] > umbral_alta_carga_ttsuv2
]

print(f"\nAnimales con carga viral de TTSuV2 en médula ósea > {umbral_alta_carga_ttsuv2} log10 copias/mg:")

if len(animales_alta_carga_ttsuv2) > 0:
    for animal in animales_alta_carga_ttsuv2:
        es_pcv2_sd_afectado = animal['estado_salud'] == 'PCV2-SD Afectado'
        infeccion_temprana = animal['edad_semanas'] <= 3 # Definimos 'temprana' como 3 semanas o menos

        print(f"- Animal ID: {animal['id_animal']}")
        print(f"  Estado de Salud: {animal['estado_salud']}")
        print(f"  Edad (semanas): {animal['edad_semanas']}")
        print(f"  Carga Viral TTSuV2 (BM): {animal['carga_viral_ttsuv2_bm']:.2f} log10 copias/mg")
        print(f"  ¿PCV2-SD Afectado? {'Sí' if es_pcv2_sd_afectado else 'No'}")
        print(f"  ¿Infección Temprana (<={infeccion_temprana} semanas)? {'Sí' if infeccion_temprana else 'No'}")
        print("-" * 30)

    # Contar y analizar específicamente aquellos con alta carga y PCV2-SD
    num_alta_carga_pcv2_sd = sum(1 for a in animales_alta_carga_ttsuv2 if a['estado_salud'] == 'PCV2-SD Afectado')
    num_alta_carga_infeccion_temprana = sum(1 for a in animales_alta_carga_ttsuv2 if a['edad_semanas'] <= 3)

    print(f"\nResumen:")
    print(f"- {num_alta_carga_pcv2_sd} animales con alta carga TTSuV2 están afectados por PCV2-SD (simulando la asociación).")
    print(f"- {num_alta_carga_infeccion_temprana} animales con alta carga TTSuV2 tuvieron infección temprana (simulando la relación con mayor propensión a PCV2-SD).")

else:
    print("No se encontraron animales con alta carga viral de TTSuV2 en médula ósea.")

# --- Paso 4: Visualización Gráfica de los Resultados ---
import matplotlib.pyplot as plt
import numpy as np

def graficar_factores_riesgo_pedv(factores, num_positivos):
    """Crea un gráfico de barras para los factores de riesgo de PEDv."""
    nombres_factores = list(factores.keys())
    conteos = list(factores.values())
    porcentajes = [(c / num_positivos) * 100 for c in conteos]

    fig, ax = plt.subplots(figsize=(12, 7))
    bars = ax.bar(nombres_factores, conteos, color='skyblue')

    ax.set_ylabel('Número de Camiones Positivos')
    ax.set_title('Factores de Riesgo en Camiones Positivos a PEDv')
    ax.set_xticks(np.arange(len(nombres_factores)))
    ax.set_xticklabels(nombres_factores, rotation=45, ha="right")

    for i, bar in enumerate(bars):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2.0, height, f'{porcentajes[i]:.1f}%', ha='center', va='bottom')

    plt.tight_layout()
    plt.savefig('pedv_factores_riesgo.png')
    print("\nGráfico de factores de riesgo para PEDv guardado como 'pedv_factores_riesgo.png'")

def graficar_analisis_ttsuv2(datos, umbral):
    """Crea un gráfico de dispersión para el análisis de TTSuV2."""
    edades = [d['edad_semanas'] for d in datos]
    cargas_virales = [d['carga_viral_ttsuv2_bm'] for d in datos]
    colores = ['red' if d['estado_salud'] == 'PCV2-SD Afectado' else 'green' for d in datos]
    tamaños = [200 if d['edad_semanas'] <= 3 else 50 for d in datos] # Más grande para infección temprana

    fig, ax = plt.subplots(figsize=(10, 6))
    scatter = ax.scatter(edades, cargas_virales, c=colores, s=tamaños, alpha=0.7)

    ax.axhline(y=umbral, color='gray', linestyle='--', label=f'Umbral Alta Carga ({umbral})')
    ax.set_xlabel('Edad (semanas)')
    ax.set_ylabel('Carga Viral TTSuV2 (log10 copias/mg)')
    ax.set_title('Análisis de Carga Viral de TTSuV2 vs. Edad y Estado de Salud')
    
    # Leyenda personalizada
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='PCV2-SD Afectado', markerfacecolor='red', markersize=10),
        Line2D([0], [0], marker='o', color='w', label='Sano', markerfacecolor='green', markersize=10),
        Line2D([0], [0], marker='o', color='w', label='Infección Temprana (<=3 sem)', markerfacecolor='gray', markersize=15, alpha=0.5),
        Line2D([0], [0], color='gray', linestyle='--', label=f'Umbral Alta Carga ({umbral})')
    ]
    ax.legend(handles=legend_elements)

    plt.tight_layout()
    plt.savefig('ttsuv2_carga_viral_analisis.png')
    print("Gráfico de análisis de TTSuV2 guardado como 'ttsuv2_carga_viral_analisis.png'")

print("--- Visualización Gráfica de los Resultados ---")
if num_pedv_positivos > 0:
    graficar_factores_riesgo_pedv(factores_de_riesgo_pedv, num_pedv_positivos)

if len(datos_animales) > 0:
    graficar_analisis_ttsuv2(datos_animales, umbral_alta_carga_ttsuv2)
