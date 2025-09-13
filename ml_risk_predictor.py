import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
import joblib
import os

# Datos simulados de camiones (PEDv)
# Extraídos de analisis_riesgo_viral.py
datos_camiones = [
    {'id_camion': 1, 'pcr_pedv_ingreso': 'Positivo', 'tipo_planta': 'Nacional', 'zona_sacrificio': 'Mayor', 'uso_exclusivo_vehiculo': 'No', 'visita_plantas_concentrado': 'Si', 'conductor_desciende': 'Si'},
    {'id_camion': 2, 'pcr_pedv_ingreso': 'Negativo', 'tipo_planta': 'Local', 'zona_sacrificio': 'Menor', 'uso_exclusivo_vehiculo': 'Si', 'visita_plantas_concentrado': 'No', 'conductor_desciende': 'No'},
    {'id_camion': 3, 'pcr_pedv_ingreso': 'Positivo', 'tipo_planta': 'Nacional-Exportacion', 'zona_sacrificio': 'Mayor', 'uso_exclusivo_vehiculo': 'No', 'visita_plantas_concentrado': 'Si', 'conductor_desciende': 'Si'},
    {'id_camion': 4, 'pcr_pedv_ingreso': 'Negativo', 'tipo_planta': 'Nacional', 'zona_sacrificio': 'Mayor', 'uso_exclusivo_vehiculo': 'Si', 'visita_plantas_concentrado': 'No', 'conductor_desciende': 'Si'},
    {'id_camion': 5, 'pcr_pedv_ingreso': 'Positivo', 'tipo_planta': 'Nacional', 'zona_sacrificio': 'Mayor', 'uso_exclusivo_vehiculo': 'No', 'visita_plantas_concentrado': 'No', 'conductor_desciende': 'Si'},
    {'id_camion': 6, 'pcr_pedv_ingreso': 'Negativo', 'tipo_planta': 'Local', 'zona_sacrificio': 'Menor', 'uso_exclusivo_vehiculo': 'Si', 'visita_plantas_concentrado': 'No', 'conductor_desciende': 'No'},
    {'id_camion': 7, 'pcr_pedv_ingreso': 'Positivo', 'tipo_planta': 'Nacional-Exportacion', 'zona_sacrificio': 'Mayor', 'uso_exclusivo_vehiculo': 'No', 'visita_plantas_concentrado': 'Si', 'conductor_desciende': 'Si'},
    {'id_camion': 8, 'pcr_pedv_ingreso': 'Positivo', 'tipo_planta': 'Nacional', 'zona_sacrificio': 'Mayor', 'uso_exclusivo_vehiculo': 'No', 'visita_plantas_concentrado': 'Si', 'conductor_desciende': 'Si'},
    {'id_camion': 9, 'pcr_pedv_ingreso': 'Negativo', 'tipo_planta': 'Local', 'zona_sacrificio': 'Menor', 'uso_exclusivo_vehiculo': 'Si', 'visita_plantas_concentrado': 'No', 'conductor_desciende': 'No'},
    {'id_camion': 10, 'pcr_pedv_ingreso': 'Positivo', 'tipo_planta': 'Nacional', 'zona_sacrificio': 'Mayor', 'uso_exclusivo_vehiculo': 'No', 'visita_plantas_concentrado': 'Si', 'conductor_desciende': 'Si'},
    {'id_camion': 11, 'pcr_pedv_ingreso': 'Positivo', 'tipo_planta': 'Nacional', 'zona_sacrificio': 'Mayor', 'uso_exclusivo_vehiculo': 'Si', 'visita_plantas_concentrado': 'Si', 'conductor_desciende': 'Si'},
    {'id_camion': 12, 'pcr_pedv_ingreso': 'Negativo', 'tipo_planta': 'Nacional-Exportacion', 'zona_sacrificio': 'Mayor', 'uso_exclusivo_vehiculo': 'No', 'conductor_desciende': 'No'},
]

# Convertir a DataFrame de pandas
df = pd.DataFrame(datos_camiones)

# Definir características (X) y variable objetivo (y)
X = df.drop(['id_camion', 'pcr_pedv_ingreso'], axis=1)
y = df['pcr_pedv_ingreso']

# Identificar columnas categóricas
categorical_features = X.select_dtypes(include=['object']).columns

# Crear un preprocesador para las columnas categóricas
preprocessor = ColumnTransformer(
    transformers=[
        ('cat', OneHotEncoder(handle_unknown='ignore'), categorical_features)
    ])

# Crear un pipeline con el preprocesador y el modelo
model_pipeline = Pipeline(steps=[
    ('preprocessor', preprocessor),
    ('classifier', RandomForestClassifier(random_state=42))
])

# Entrenar el modelo
print("--- Entrenamiento del Modelo de Predicción de Riesgo PEDv ---")
model_pipeline.fit(X, y)
print("Modelo entrenado exitosamente.")

# Guardar el modelo y el preprocesador
model_filename = 'pedv_risk_predictor_model.joblib'
joblib.dump(model_pipeline, model_filename)
print(f"Modelo guardado como '{model_filename}'")

# --- Función para predecir nuevos factores de riesgo ---
def predict_new_risk_factors(model, new_data):
    # Convertir la nueva data a DataFrame
    new_df = pd.DataFrame([new_data])
    prediction = model.predict(new_df)
    return prediction[0]

# --- Predicción para un nuevo camión ---
print("\n--- Predicción de Riesgo para un Nuevo Camión ---")

# Cargar el modelo
loaded_model = joblib.load(model_filename)

# Datos del nuevo camión
print("Por favor, ingrese los datos del nuevo camión:")
tipo_planta = input("Tipo de planta (e.g., Nacional, Local, Nacional-Exportacion): ")
zona_sacrificio = input("Zona de sacrificio (e.g., Mayor, Menor): ")
uso_exclusivo_vehiculo = input("¿Uso exclusivo del vehículo? (Si/No): ")
visita_plantas_concentrado = input("¿Visita plantas de concentrado? (Si/No): ")
conductor_desciende = input("¿El conductor desciende del vehículo? (Si/No): ")

new_truck_data = {
    'tipo_planta': tipo_planta,
    'zona_sacrificio': zona_sacrificio,
    'uso_exclusivo_vehiculo': uso_exclusivo_vehiculo,
    'visita_plantas_concentrado': visita_plantas_concentrado,
    'conductor_desciende': conductor_desciende
}

print(f"Datos del nuevo camión: {new_truck_data}")

# Predecir el riesgo
prediction = predict_new_risk_factors(loaded_model, new_truck_data)
print(f"\nLa predicción de riesgo de PEDv para este camión es: {prediction}")
print("--------------------------------------------------")
