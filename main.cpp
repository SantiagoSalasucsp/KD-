import time
import matplotlib.pyplot as plt
from wordcloud import WordCloud
import os

# Ruta al archivo de comunicación
PATH = "/Users/santiagosalas/Desktop/python/mensaje.txt"

# Estado de control
last_file_size = -1 
last_frequencies = {}

plt.ion() # modo interactivo para refrescar la nube
fig, ax = plt.subplots(figsize=(10, 5)) # Crear figura y eje una vez

def generar_nube(frecuencias):
    """Genera y muestra la nube de palabras a partir de un diccionario de frecuencias."""
    if not frecuencias:
        ax.text(0.5, 0.5, "Esperando datos del productor C++...", 
                ha='center', va='center', fontsize=20, color='gray')
        ax.axis("off")
        return

    # WordCloud pide un diccionario {palabra: frecuencia}
    wc = WordCloud(
        width=800,
        height=400,
        background_color="white",
        colormap="viridis"
    ).generate_from_frequencies(frecuencias)

    # Mostrar la nube
    ax.clear()
    ax.imshow(wc, interpolation="bilinear")
    ax.axis("off")
    ax.set_title("WordCloud de Temas (Acumulativo)", fontsize=16)
    plt.pause(0.01)

# Inicializar la nube con un mensaje de espera
generar_nube(last_frequencies)

print("--- INICIANDO CONSUMIDOR PYTHON ---")
print("Escuchando archivo:", PATH)

while True:
    try:
        # 1. Verificar si el archivo ha cambiado de tamaño
        current_file_size = os.stat(PATH).st_size
        
        if current_file_size > last_file_size:
            print(f"Cambiando de tamaño: {last_file_size} -> {current_file_size}. Leyendo...")
            last_file_size = current_file_size
            
            # 2. Leer TODAS las líneas del archivo (acumulación)
            # El archivo C++ está en modo append, por lo que leer todo acumula los datos.
            with open(PATH, "r") as f:
                lines = f.readlines()
            
            # 3. Procesar las líneas CSV
            current_frequencies = {}
            for line in lines:
                try:
                    # El formato es: palabra,frecuencia\n
                    parts = line.strip().split(',')
                    if len(parts) == 2:
                        palabra = parts[0]
                        frecuencia = int(parts[1])
                        # Acumulamos las frecuencias si la misma palabra aparece varias veces
                        current_frequencies[palabra] = current_frequencies.get(palabra, 0) + frecuencia
                except ValueError:
                    # Ignorar líneas mal formadas o vacías
                    continue

            # 4. Actualizar la nube de palabras
            last_frequencies = current_frequencies
            generar_nube(last_frequencies)
            
            print(f"Actualización exitosa. Total de palabras únicas: {len(last_frequencies)}")

        elif current_file_size == last_file_size and current_file_size > 0:
            # Archivo no ha cambiado, no actualizar la gráfica.
            # Este es el comportamiento deseado cuando C++ ha terminado de escribir.
            print("Archivo sin cambios. Esperando nuevo dato...")

        else:
             # El archivo está vacío o no existe aún
            print("Esperando a que C++ cree/escriba en el archivo...")


    except FileNotFoundError:
        # El archivo aún no existe (típico al inicio)
        print("Archivo no encontrado. Esperando a que el productor C++ inicie...")
    except Exception as e:
        print(f"Error inesperado: {e}")

    time.sleep(1) # Revisa cada segundo
