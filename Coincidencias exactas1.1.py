import pandas as pd
import re 

# Paso 1: Generar el diccionario de péptidos de alfa-sinucleína

# Definir la nueva secuencia de la proteína alfa-sinucleína
sequence = "MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA"

# Generar péptidos a partir de la secuencia
def generate_peptides(sequence, min_length=8, max_length=21):
    peptides = []
    seq_length = len(sequence)
    
    for length in range(min_length, max_length + 1):
        for start in range(seq_length - length + 1):
            peptide = sequence[start:start + length]
            peptides.append(peptide)
    
    return peptides

# Función para renombrar los péptidos
def rename_peptide(peptide_sequence, full_sequence):
    start_position = full_sequence.find(peptide_sequence)
    if start_position == -1:
        return "Not found"
    
    start_position += 1
    peptide_length = len(peptide_sequence)
    new_name = f"{start_position}P{peptide_length}pep"
    return new_name

# Generar péptidos y crear un DataFrame
peptides = generate_peptides(sequence)
data = {
    'New Name': [rename_peptide(peptide, sequence) for peptide in peptides],
    'Epitope - Name': peptides
}
df_peptides = pd.DataFrame(data)

# Guardar el DataFrame como archivo FASTA
output_fasta_path_1 = r"C:\Users\dennys\Desktop\Resultados oficiales\Peptidos AlphaSyn.fasta"

def save_as_fasta(df, output_fasta_path):
    with open(output_fasta_path, 'w') as fasta_file:
        for index, row in df.iterrows():
            new_name = row['New Name']
            sequence = row['Epitope - Name']
            fasta_file.write(f">{new_name}\n")
            fasta_file.write(f"{sequence}\n")

save_as_fasta(df_peptides, output_fasta_path_1)
print(f"Archivo FASTA guardado en: {output_fasta_path_1}")

# Convertir el archivo FASTA a diccionario y guardarlo como TXT
output_dict_path_1 = r"C:\Users\dennys\Desktop\Resultados oficiales\Peptidos AlphaSyn diccionario.txt"
diccionario_peptidos_1 = {}

with open(output_fasta_path_1, 'r') as archivo:
    nombre_peptido = None
    for linea in archivo:
        linea = linea.strip()
        if linea.startswith(">"):
            nombre_peptido = linea[1:]  # Guardar el nombre del péptido sin '>'
        else:
            diccionario_peptidos_1[nombre_peptido] = linea

# Guardar el diccionario en un archivo TXT
with open(output_dict_path_1, 'w') as archivo_txt:
    for nombre, secuencia in diccionario_peptidos_1.items():
        archivo_txt.write(f"{nombre}: {secuencia}\n")

print(f"El diccionario de péptidos se ha guardado en {output_dict_path_1}.")


# Paso 2: Procesar el archivo Excel de epítopos de Alfa-Sinucleína comprobados experimentalmente

# Ruta del archivo Excel
excel_file_path = r"C:\Users\dennys\Desktop\Resultados oficiales\AlphaSyn epítopos MHC II mouse + PTM.xlsx"

# Cargar el archivo Excel
df = pd.read_excel(excel_file_path)

# Función para eliminar modificaciones PTM de la secuencia
def clean_peptide_sequence(peptide_sequence):
    # Limpiar el epítopo, eliminando todo después del "+"
    cleaned_sequence = re.sub(r'\+.*', '', peptide_sequence).strip()
    return cleaned_sequence

# Función para extraer las modificaciones PTM (todo lo que está después del "+")
def extract_ptm(peptide_sequence):
    ptm_match = re.search(r'\+\s*(.*)', peptide_sequence)
    if ptm_match:
        return ptm_match.group(1).strip()  # Captura la modificación después del "+"
    return ""

# Función para renombrar los péptidos con PTM y posición inicial
def rename_peptide_with_ptm_and_position(peptide_sequence, starting_position):
    # Limpiar el epítopo para quitar las modificaciones PTM
    cleaned_peptide = clean_peptide_sequence(peptide_sequence)
    
    # Extraer las modificaciones PTM
    ptm = extract_ptm(peptide_sequence)
    
    # Obtener la longitud del epítopo limpio
    peptide_length = len(cleaned_peptide)
    
    # Verificar si hay una posición de inicio proporcionada
    if pd.notna(starting_position):
        # Generar el nombre en formato "PosiciónPLongitudpep"
        new_name = f"{int(starting_position)}P{peptide_length}pep"
    else:
        # Si no hay posición de inicio, solo usa la longitud
        new_name = f"{peptide_length}pep"
    
    # Añadir la modificación PTM si existe (sin espacio)
    if ptm:
        new_name += ptm
    
    return new_name, cleaned_peptide

# Aplicar la función a cada fila del DataFrame
df[['New Name', 'Cleaned Sequence']] = df.apply(
    lambda row: pd.Series(rename_peptide_with_ptm_and_position(row['Epitope - Name'], row['Epitope - Starting Position'])),
    axis=1
)

# Guardar el archivo modificado como Excel
output_excel_path = r"C:\Users\dennys\Desktop\Resultados oficiales\Epítopos comprobados MHC I humano renombrados.xlsx"
df.to_excel(output_excel_path, index=False)
print(f"Archivo Excel guardado en: {output_excel_path}")

# Guardar el DataFrame modificado como archivo FASTA
output_fasta_path_2 = r"C:\Users\dennys\Desktop\Resultados oficiales\Alpha Syn MHC II mouse PTM.fasta"

def save_as_fasta(df, output_fasta_path):
    with open(output_fasta_path, 'w') as fasta_file:
        for index, row in df.iterrows():
            new_name = row['New Name']
            sequence = row['Cleaned Sequence']
            fasta_file.write(f">{new_name}\n")
            fasta_file.write(f"{sequence}\n")

save_as_fasta(df, output_fasta_path_2)
print(f"Archivo FASTA guardado en: {output_fasta_path_2}")

# Convertir el archivo FASTA de epítopos a diccionario
output_dict_path_2 = r"C:\Users\dennys\Desktop\Resultados oficiales\Epítopos comprobados MHC I humano diccionario.txt"
diccionario_peptidos_2 = {}

with open(output_fasta_path_2, 'r') as archivo:
    nombre_peptido = None
    for linea in archivo:
        linea = linea.strip()
        if linea.startswith(">"):
            nombre_peptido = linea[1:]  # Guardar el nombre del péptido sin '>'
        else:
            diccionario_peptidos_2[nombre_peptido] = linea

# Guardar el diccionario en un archivo TXT
with open(output_dict_path_2, 'w') as archivo_txt:
    for nombre, secuencia in diccionario_peptidos_2.items():
        archivo_txt.write(f"{nombre}: {secuencia}\n")

print(f"El diccionario de péptidos se ha guardado en {output_dict_path_2}.")

# Paso 3: Buscar ocurrencias exactas entre los diccionarios de peptidos de alfa sinucleina generados con epítopos comprobados experimentalmente

# Rutas de los diccionarios y resultados
archivo_diccionario1 = output_dict_path_2  # Diccionario del paso 2
archivo_diccionario2 = output_dict_path_1  # Diccionario del paso 1
archivo_salida_txt = r"C:\Users\dennys\Desktop\Resultados oficiales\Resultados comparacion AlphaSyn MHC II mouse PTM.txt"
archivo_salida_excel = r"C:\Users\dennys\Desktop\Resultados oficiales\Resultados comparacion AlphaSyn MHC II mouse PTM.xlsx"
archivo_salida_diccionario_coincidentes = r"C:\Users\dennys\Desktop\Resultados oficiales\Diccionario Coincidencias epitopos AlphaSyn.txt"

# Función para cargar un diccionario desde un archivo
def cargar_diccionario(archivo):
    diccionario = {}
    try:
        with open(archivo, 'r') as f:
            for linea in f:
                linea = linea.strip()
                if linea and ": " in linea:
                    nombre, secuencia = linea.split(": ", 1)
                    diccionario[nombre] = secuencia
    except FileNotFoundError as e:
        print(f"Error: {e}")
    return diccionario

# Función para remover modificaciones PTM
def remover_modificaciones(secuencia):
    return re.sub(r"\s*\+.*", "", secuencia)

# Cargar los diccionarios
diccionario1 = cargar_diccionario(archivo_diccionario1)
diccionario2 = cargar_diccionario(archivo_diccionario2)

# Crear un diccionario sin modificaciones para diccionario2
diccionario2_sin_modificaciones = {nombre: remover_modificaciones(secuencia) for nombre, secuencia in diccionario2.items()}

# Comparar las secuencias
secuencias_coincidentes = []
secuencias_no_en_diccionario2 = []
diccionario_coincidentes = {}  # Nuevo diccionario para las secuencias coincidentes

for nombre, secuencia in diccionario1.items():
    secuencia_sin_modificaciones = remover_modificaciones(secuencia)
    if secuencia_sin_modificaciones in diccionario2_sin_modificaciones.values():
        secuencias_coincidentes.append((nombre, secuencia))
        diccionario_coincidentes[nombre] = secuencia  # Agregar a diccionario de coincidentes
    else:
        secuencias_no_en_diccionario2.append((nombre, secuencia))

# Guardar resultados en un archivo TXT
with open(archivo_salida_txt, 'w') as f:
    f.write("Péptidos con secuencias coincidentes:\n")
    for nombre, secuencia in secuencias_coincidentes:
        f.write(f"{nombre}: {secuencia}\n")
    
    f.write("\nPéptidos que NO coinciden en el segundo diccionario:\n")
    for nombre, secuencia in secuencias_no_en_diccionario2:
        f.write(f"{nombre}: {secuencia}\n")

# Guardar resultados en un archivo Excel
resultados = {
    "Coincidentes": pd.DataFrame(secuencias_coincidentes, columns=["Nombre", "Secuencia"]),
    "No Coincidentes": pd.DataFrame(secuencias_no_en_diccionario2, columns=["Nombre", "Secuencia"]),
}

with pd.ExcelWriter(archivo_salida_excel, engine='openpyxl') as writer:
    resultados["Coincidentes"].to_excel(writer, sheet_name="Coincidentes", index=False)
    resultados["No Coincidentes"].to_excel(writer, sheet_name="No Coincidentes", index=False)

# Guardar el diccionario de coincidentes en un archivo TXT
with open(archivo_salida_diccionario_coincidentes, 'w') as archivo_coincidentes:
    for nombre, secuencia in diccionario_coincidentes.items():
        archivo_coincidentes.write(f"{nombre}: {secuencia}\n")

# Mensaje de confirmación
print(f"Los resultados se han guardado en {archivo_salida_txt}, {archivo_salida_excel}, y {archivo_salida_diccionario_coincidentes}.")


# Paso 4: Comparación resultados coincidencias de epítopos de alfa sinucleina con proteoma de bacteria
# Función para convertir un archivo FASTA a diccionario
def fasta_to_dict(file_path):
    fasta_dict = {}
    with open(file_path, 'r') as file:
        sequence_id = None
        sequence = []

        for line in file:
            line = line.strip()
            if line.startswith('>'):  # Línea de encabezado
                if sequence_id is not None:  # Guardar la secuencia anterior
                    fasta_dict[sequence_id] = ''.join(sequence)
                sequence_id = line[1:].split()[0]  # Obtener el ID de la secuencia
                sequence = []  # Reiniciar la lista de secuencia
            else:
                sequence.append(line)  # Agregar líneas de secuencia

        # Guardar la última secuencia en el diccionario
        if sequence_id is not None:
            fasta_dict[sequence_id] = ''.join(sequence)

    return fasta_dict

# Función para cargar un diccionario desde un archivo TXT
def cargar_diccionario(archivo):
    diccionario = {}
    try:
        with open(archivo, 'r') as f:
            for linea in f:
                linea = linea.strip()
                if linea and ": " in linea:
                    nombre, secuencia = linea.split(": ", 1)
                    diccionario[nombre] = secuencia
    except FileNotFoundError as e:
        print(f"Error: {e}")
    return diccionario

# Función para remover modificaciones PTM
def remover_modificaciones(secuencia):
    return re.sub(r"\s*\+.*", "", secuencia)

# Convertir el archivo FASTA a diccionario
input_file_path = r"C:\Users\dennys\Desktop\Resultados oficiales\UP000001031 A. muciniphila.fasta"
proteome_dict = fasta_to_dict(input_file_path)

# Guardar el diccionario de secuencias proteicas 
output_file_path_bacteria = r"C:\Users\dennys\Desktop\Resultados oficiales\Proteoma_diccionario A. muciniphila.txt"
with open(output_file_path_bacteria, 'w') as file:
    for key, value in proteome_dict.items():
        file.write(f"{key}: {value}\n")  # Formato: "nombre: secuencia"

# Cargar los diccionarios y comparar
archivo_diccionario1 = archivo_salida_diccionario_coincidentes  # Diccionario del paso 2
archivo_diccionario2 = output_file_path_bacteria  # Usar el diccionario que acabamos de generar

# Rutas para guardar los resultados
archivo_salida_txt = r"C:\Users\dennys\Desktop\Resultados oficiales\Resultados comparacion A. muciniphila.txt"
archivo_salida_excel = r"C:\Users\dennys\Desktop\Resultados oficiales\Resultados comparacion A. muciniphila.xlsx"


# Cargar diccionarios
diccionario1 = cargar_diccionario(archivo_diccionario1)
diccionario2 = cargar_diccionario(archivo_diccionario2)

# Crear diccionarios sin modificaciones PTM
diccionario1_sin_modificaciones = {nombre: remover_modificaciones(secuencia) for nombre, secuencia in diccionario1.items()}
diccionario2_sin_modificaciones = {nombre: remover_modificaciones(secuencia) for nombre, secuencia in diccionario2.items()}

# Comparar las secuencias
secuencias_coincidentes_exactas = []
secuencias_coincidentes_parciales = []
secuencias_no_en_diccionario2 = []
diccionario_coincidentes = {}

for nombre1, secuencia1 in diccionario1_sin_modificaciones.items():
    coincidencia_exacta = False
    coincidencia_parcial = False
    
    for nombre2, secuencia2 in diccionario2_sin_modificaciones.items():
        if secuencia1 == secuencia2:
            secuencias_coincidentes_exactas.append((nombre1, diccionario1[nombre1]))  # Usamos la secuencia original
            diccionario_coincidentes[nombre1] = diccionario1[nombre1]
            coincidencia_exacta = True
            break
        elif secuencia1 in secuencia2 or secuencia2 in secuencia1:
            secuencias_coincidentes_parciales.append((nombre1, diccionario1[nombre1], nombre2, diccionario2[nombre2]))  # Guardar las originales
            diccionario_coincidentes[nombre1] = diccionario1[nombre1]
            coincidencia_parcial = True
    
    if not coincidencia_exacta and not coincidencia_parcial:
        secuencias_no_en_diccionario2.append((nombre1, diccionario1[nombre1]))  # Guardar la original

# Guardar resultados en un archivo TXT
with open(archivo_salida_txt, 'w') as f:
    f.write("Péptidos con secuencias coincidentes exactas:\n")
    for nombre, secuencia in secuencias_coincidentes_exactas:
        f.write(f"{nombre}: {secuencia}\n")
    
    f.write("\nPéptidos con coincidencias parciales:\n")
    for nombre1, secuencia1, nombre2, secuencia2 in secuencias_coincidentes_parciales:
        f.write(f"{nombre1}: {secuencia1} (Coincide parcialmente con {nombre2}: {secuencia2})\n")
    
    f.write("\nPéptidos que NO coinciden en el segundo diccionario:\n")
    for nombre, secuencia in secuencias_no_en_diccionario2:
        f.write(f"{nombre}: {secuencia}\n")

# Guardar resultados en un archivo Excel
resultados = {
    "Coincidentes Exactos": pd.DataFrame(secuencias_coincidentes_exactas, columns=["Nombre", "Secuencia"]),
    "Coincidencias Parciales": pd.DataFrame(secuencias_coincidentes_parciales, columns=["Nombre Dic1", "Secuencia Dic1", "Nombre Dic2", "Secuencia Dic2"]),
    "No Coincidentes": pd.DataFrame(secuencias_no_en_diccionario2, columns=["Nombre", "Secuencia"]),
}

with pd.ExcelWriter(archivo_salida_excel, engine='openpyxl') as writer:
    resultados["Coincidentes Exactos"].to_excel(writer, sheet_name="Coincidentes Exactos", index=False)
    resultados["Coincidencias Parciales"].to_excel(writer, sheet_name="Coincidencias Parciales", index=False)
    resultados["No Coincidentes"].to_excel(writer, sheet_name="No Coincidentes", index=False)

# Mensaje de confirmación
print(f"Los resultados se han guardado en {archivo_salida_txt}, {archivo_salida_excel}, y {archivo_salida_diccionario_coincidentes}.")
