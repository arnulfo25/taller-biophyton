import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

# Datos simulados de animales (TTSuV2)
# Extraídos de analisis_riesgo_viral.py
datos_animales = [
    {'id_animal': 101, 'estado_salud': 'PCV2-SD Afectado', 'edad_semanas': 3, 'carga_viral_ttsuv2_bm': 8.2},
    {'id_animal': 102, 'estado_salud': 'Sano', 'edad_semanas': 3, 'carga_viral_ttsuv2_bm': 6.5},
    {'id_animal': 103, 'estado_salud': 'PCV2-SD Afectado', 'edad_semanas': 7, 'carga_viral_ttsuv2_bm': 8.5},
    {'id_animal': 104, 'estado_salud': 'Sano', 'edad_semanas': 11, 'carga_viral_ttsuv2_bm': 7.0},
    {'id_animal': 105, 'estado_salud': 'PCV2-SD Afectado', 'edad_semanas': 1, 'carga_viral_ttsuv2_bm': 7.9},
    {'id_animal': 106, 'estado_salud': 'Sano', 'edad_semanas': 8, 'carga_viral_ttsuv2_bm': 7.2},
    {'id_animal': 107, 'estado_salud': 'PCV2-SD Afectado', 'edad_semanas': 15, 'carga_viral_ttsuv2_bm': 9.0},
    {'id_animal': 108, 'estado_salud': 'Sano', 'edad_semanas': 12, 'carga_viral_ttsuv2_bm': 6.8},
    {'id_animal': 109, 'estado_salud': 'PCV2-SD Afectado', 'edad_semanas': 5, 'carga_viral_ttsuv2_bm': 8.1},
    {'id_animal': 110, 'estado_salud': 'Sano', 'edad_semanas': 10, 'carga_viral_ttsuv2_bm': 6.9},
]

# Preparar los datos para la regresión lineal
# X: Edad en semanas (variable independiente)
# y: Carga viral TTSuV2 (variable dependiente)
X = np.array([d['edad_semanas'] for d in datos_animales]).reshape(-1, 1)
y = np.array([d['carga_viral_ttsuv2_bm'] for d in datos_animales])

# Crear y entrenar el modelo de regresión lineal
model = LinearRegression()
model.fit(X, y)

# Realizar predicciones
y_pred = model.predict(X)

# Imprimir los coeficientes del modelo
print("--- Análisis de Regresión Lineal para Carga Viral TTSuV2 ---")
print(f"Coeficiente (pendiente): {model.coef_[0]:.2f}")
print(f"Intersección (ordenada al origen): {model.intercept_:.2f}")
print(f"R^2 (coeficiente de determinación): {model.score(X, y):.2f}")

# Función para predecir la carga viral para una nueva edad
def predecir_carga_viral(edad_semanas):
    return model.predict(np.array([[edad_semanas]]))[0]

# Ejemplo de predicción
edad_ejemplo = 9
carga_predicha = predecir_carga_viral(edad_ejemplo)
print(f"\nPredicción para una edad de {edad_ejemplo} semanas: {carga_predicha:.2f} log10 copias/mg")

# Visualización de los resultados
plt.figure(figsize=(10, 6))
plt.scatter(X, y, color='blue', label='Datos Reales')
plt.plot(X, y_pred, color='red', label='Línea de Regresión')
plt.title('Regresión Lineal: Edad vs. Carga Viral TTSuV2')
plt.xlabel('Edad (Semanas)')
plt.ylabel('Carga Viral TTSuV2 (log10 copias/mg)')
plt.legend()
plt.grid(True)
plt.savefig('ttsuv2_linear_regression.png')
print("\nGráfico de regresión lineal guardado como 'ttsuv2_linear_regression.png'")

# --- Sección para que el usuario pueda predecir ---
print("\n--- Predicción Interactiva ---")
while True:
    try:
        user_input = input("Ingrese una edad en semanas para predecir la carga viral (o 'salir' para terminar): ")
        if user_input.lower() == 'salir':
            break
        edad_usuario = float(user_input)
        if edad_usuario < 0:
            print("La edad no puede ser negativa. Intente de nuevo.")
            continue
        prediccion_usuario = predecir_carga_viral(edad_usuario)
        print(f"Para {edad_usuario} semanas, la carga viral predicha es: {prediccion_usuario:.2f} log10 copias/mg")
    except ValueError:
        print("Entrada inválida. Por favor, ingrese un número o 'salir'.")
    except Exception as e:
        print(f"Ocurrió un error: {e}")
