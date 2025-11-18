# Proyecto_DESeq2

# Resultados de análisis — comparacion_R_vs_ATLAS_resultados.csv
Antes de correr los comandos asegurese de estar en la carpeta donde tiene guardado el archivo que se va a analizar.   Este archivo es el archivo resultante luego de correr el código de comparación de datos:
```bash
comparacion.ipynb
```

## Conteos en R

Cuenta cuantos Up hay en R  
```bash
awk -F',' '$7 ~ /Up/ {c++} END{print c}' comparacion_R_vs_ATLAS_resultados.csv  
295
```
Cuenta cuantos Down hay en R
```bash
awk -F';' '$7=="Down" {c++} END{print c}' comparacion_R_vs_ATLAS_resultados.csv  
678
```
Cuenta cuantos NS hay en R  
```bash
awk -F’, '$7=="NS" {c++} END{print c}' comparacion_R_vs_ATLAS_resultados.csv  
325
```

## Coincidencias (Match)

Cuenta cuantos match FALSE hay  
```bash
grep -i -c 'false' comparacion_R_vs_ATLAS_resultados.csv  
326
```

Cuenta cuantos match TRUE hay  
```bash
grep -i -c 'true' comparacion_R_vs_ATLAS_resultados.csv  
972
```

## Conteos en ATLAS

Cuenta cuantos genes significativos fueron hallados en atlas  
```bash
wc -l comparacion_R_vs_ATLAS_resultados.csv  
1298
```

Cuenta cuantos Up hay en Atlas  
```bash
awk -F’,’ '$8=="Up" {c++} END{print c}' comparacion_R_vs_ATLAS_resultados.csv  
392
```

Cuenta cuantos Down hay en Atlas  
```bash
awk -F’,’ '$8=="Down" {c++} END{print c}' comparacion_R_vs_ATLAS_resultados.csv  
906
```
