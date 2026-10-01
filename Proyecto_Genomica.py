#Proyecto: Analizador de Variantes Genéticas (mini-pipeline de genómica)

#Contexto real: en genómica, uno de los análisis más comunes es comparar una secuencia de un paciente/muestra contra una 
#secuencia de referencia, para detectar mutaciones puntuales (SNPs) — esto es la base de cómo se detectan variantes asociadas 
#a enfermedades, resistencia a antibióticos, o mutaciones virales (como las variantes de SARS-CoV-2).

#Objetivo del proyecto: construir un script que reciba una secuencia de referencia y una secuencia de muestra 
#(ambas de la misma longitud), las compare, detecte las mutaciones, las clasifique, y genere un reporte.

import matplotlib.pyplot as plt
import pandas as pd 
import seaborn as sns 

#Funcionalidades a implementar:
#1. Validación

#Función que verifique que ambas secuencias tengan la misma longitud y solo contengan bases válidas.

def validar_secuencias_dna(referencia,muestra):
    ref = referencia.upper()
    mue = muestra.upper()
    
    #Verificar que tengan la misma longitud
    if len(ref) != len(mue):
        print(f"Error. Las longitudes no coinciden ({len(ref)} vs {len(mue)})")
        return False
    
    #Verificar que solo contengan bases validas A,T,G,C
    bases_validas = {"A", "T", "G", "C"}
    
    if not set(ref).issubset(bases_validas):
        print(f"Error: La secuencia de referencia contiene caracteres no validos.")
        return False
    
    if not set(mue).issubset(bases_validas):
            print(f"Error: La secuencia de muestre contiene caracteres no validos.")
            return False
    
    return True


#2. Detección de mutaciones (SNPs)
#Crea una función que compare ambas secuencias posición por posición
#Cuando encuentre una diferencia, guarda esa información (posición, base de referencia, base mutada) en una lista de diccionarios
#Debe retornar esa lista completa de mutaciones encontradas
def detectar_mutaciones_snps(referencia,muestra):
    mutaciones = []
    for i, (base_ref, base_mue) in enumerate(zip(referencia,muestra)):
        if base_ref != base_mue:
            snp = {
                "posicion": i,
                "referencia": base_ref,
                "mutada": base_mue
            }
            mutaciones.append(snp)
    return mutaciones
            
#3. Clasificación de mutaciones
#Para cada mutación detectada, identifica a qué codón pertenece esa posición
#Traduce el codón original y el codón mutado a sus aminoácidos correspondientes (usando tu código genético)
#Clasifica la mutación como:
#Silenciosa (mismo aminoácido)
#Missense (aminoácido distinto)
#Nonsense (el nuevo codón es STOP)

def clasificar_mutaciones(mutaciones,referencia,muestra):
    codigo_genetico = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGT':'R', 'CGG':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'STOP', 'TAG':'STOP',
        'TGC':'C', 'TGT':'C', 'TGA':'STOP', 'TGG':'W',
    }
    
    mutaciones_clasificadas = []
    
    for mutacion in mutaciones:
        pos = mutacion["posicion"]
        inicio_codon = (pos // 3)*3
        
        codon_ref = referencia[inicio_codon : inicio_codon + 3]
        codon_mue = muestra[inicio_codon : inicio_codon + 3]
        
        aa_ref = codigo_genetico.get(codon_ref, "?")
        aa_mue = codigo_genetico.get(codon_mue, "?")
        
        if aa_mue == "STOP" and aa_ref != "STOP":
            tipo = "Nonsense"
        elif aa_ref == aa_mue:
            tipo = "Silenciosa"
        else:
            tipo = "Missense"
            
        mutaciones_clasificadas.append({
            "posicion": pos,
            "ref_base": mutacion["referencia"],
            "mue_base": mutacion["mutada"],
            "codon_ref": codon_ref,
            "codon_mue": codon_mue,
            "aa_ref": aa_ref,
            "aa_mue": aa_mue,
            "tipo": tipo
        })
        
    return mutaciones_clasificadas


# ----- PRUEBA DEL CODIGO -----

referencia = "ATGGCCGACTTGCAACCGTGA"
muestra    = "ATGGCAGACTTGCAACCGTAA"

# PASO 0: VALIDAR (te falta este paso)
if not validar_secuencias_dna(referencia, muestra):
    print("No se puede continuar: secuencias inválidas.")
    exit()  # o usa un return si esto estuviera dentro de una función main()

# PASO 1: DETECTAR
snps = detectar_mutaciones_snps(referencia,muestra)

# PASO 2: CLASIFICAR
resultado = clasificar_mutaciones(snps,referencia,muestra)
tasa_mutacion = (len(resultado) / len(referencia)) * 100
print(f"\nTasa de mutación: {tasa_mutacion:.2f}%")

#MOSTRAR RESULTADOS
for m in resultado:
    print(m)

import pandas as pd

#Convertir la lista de diccionarios a un DataFrame de Pandas
df_mutaciones = pd.DataFrame(resultado)

#Mostrar tabla foramteada en consola
print("\n--- REPORTE DE MUTACIONES (SNPs) ---")
print(df_mutaciones)

#Exportar el archivo a un CSV
df_mutaciones.to_csv("reporte_mutaciones.csv", index=False)

#Convertir los resultados a un DataFrame
df = pd.DataFrame(resultado)

#Configurar el estilo visual
sns.set_theme(style="whitegrid")
plt.figure(figsize=(5,8))

#Crear Grafico de barras contando el Tipo de mutacion
# ...existing code...
ax = sns.countplot(
    data=df,
    x="tipo",
    hue="tipo",
    palette="Set2",
    order=["Silenciosa", "Missense", "Nonsense"]
)

plt.title("Frecuencia por tipo de mutación (SNPs)", fontsize=14, fontweight="bold")
plt.xlabel("Tipo de mutación", fontsize=12)
plt.ylabel("Cantidad de mutaciones", fontsize=12)

#Configurar eje Y para mostrar numeros enteros
plt.yticks(range(0, len(df) + 2))

#Guardar la iamgen en la carpeta de trabajo
plt.tight_layout()
plt.savefig("grafico_mutaciones.png", dpi=300)
plt.show()


#GENERAR REPORTE A PARTIR DE UNA FUNCION
def generar_reporte_completo(referencia,muestra):
    referencia = referencia.upper()
    muestra = muestra.upper()
    
    #PASO 1: VALIDAR
    if not validar_secuencias_dna(referencia,muestra):
        print("No se puede contoinuar: secuencias invalidas.")
        return None
    
    #DETECTAR
    snps = detectar_mutaciones_snps(referencia, muestra)
    
    #CLASIFICAR
    resultado = clasificar_mutaciones(snps, referencia, muestra)
    
    #TASA DE MUTACION
    tasa_mutacion = (len(resultado) / len(referencia)) * 100
    
    #REPORTE EN CONSOLA
    print("\n" + "=" * 50)
    print("REPORTE DE ANALISIS DE VARIANTES GENÉTICAS")
    print("=" * 50)
    print(f"Longitud de secuencia: {len(referencia)} bases")
    print(f"Total de mutaciones detectada: {len(resultado)}")
    print(f"Tasa de mutación: {tasa_mutacion:.2f}%")
    print("-" * 50)
    for m in resultado:
        print(m)
    print("-" * 50)
    
    #EXPORTAR A CSV
    df_mutaciones = pd.DataFrame(resultado)
    df_mutaciones.to_csv("reporte_mutaciones.csv", index=False)
    print("\nReporte exportado a 'reporte_mutaciones.csv'")
    
    #GRAFICAR: SOLO SI HAY MUTACIONES QUE MOSTRAR
    if leb(resultado) > 0:
        sns.set_theme(style="Whitegrid")
        plt.figure(figsize=(5,8))
        
        ax = sns.countplot(
            data=df_mutaciones,
            x="tipo",
            hue="tipo",
            palette="Set2",
            order=["Silenciosa", "Missense", "Nonsense"],
            legend=False
        )
        
        plt.title("Frecuencia por tipo de mutación (SNPs)", fontsize=14, fontweight="bold")
        plt.xlabel("Tipo de mutación", fontsize=12)
        plt.ylabel("Cantidad de mutaciones", fontsize=12)
        plt.yticks(range(0, len(df_mutaciones) + 2))
    
        plt.tight_layout()
        plt.savefig("grafico_mutaciones.png", dpi=300)
        plt.show()
    
    else:
        print("NO hay mutaciones para graficar.")
    
#RETORNAR EL RESUTLADO POR SI QUIERO SEGUIR USANDOLO
    return{
        "mutaciones": resultado,
        "tasa_mutacion": tasa_mutacion,
        "dataframe": df_mutaciones
    }
    
# ----- PRUEBA DEL CÓDIGO -----
referencia = "ATGGCCGACTTGCAACCGTGA"
muestra    = "ATGGCAGACTTGCAACCGTAA"

reporte = generar_reporte_completo(referencia, muestra)