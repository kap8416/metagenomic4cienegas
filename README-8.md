
# Viroma Functional Annotation – Cuatro Ciénegas, Coahuila

Este repositorio contiene scripts y visualizaciones para el análisis funcional de virus identificados mediante metagenómica en muestras de tapetes microbianos y estromatolitos de Cuatro Ciénegas, Coahuila.

## Descripción del Proyecto

- **Objetivo**: Caracterizar el perfil funcional de virus detectados en diferentes muestras ambientales.
- **Ubicación**: Cuatro Ciénegas, Coahuila, México
- **Tipo de muestras**:
  - 6 muestras de **tapete microbiano**
  - 4 muestras de **estromatolito**

Cada muestra ha sido analizada mediante un flujo de trabajo de metagenómica viral y anotada funcionalmente con *eggNOG-mapper*. Se presentan las visualizaciones para cuatro tipos de anotaciones:

- Categorías funcionales **COG**
- Términos **KEGG Orthologs (KO)**
- **Vías metabólicas KEGG**
- Dominios **Pfam**

## Contenido del Repositorio

```
.
├── viroma_annotation_plots.py      # Script para generar paneles funcionales
├── figures/                        # Carpeta de salida con figuras PNG generadas
├── data/                           # Carpeta sugerida para almacenar archivos tabulares de eggNOG
└── README.md                       # Descripción del proyecto
```

## Requisitos

- Python 3.7+
- pandas
- matplotlib

Instalación de dependencias:
```bash
pip install pandas matplotlib
```

## Uso del Script

Coloca los archivos `.tabular` de salida de eggNOG-mapper en una carpeta. Luego, ejecuta el script pasando el nombre de archivo y el identificador de la muestra:

```python
from viroma_annotation_plots import generate_panel_plot

# Ejemplo de uso
generate_panel_plot("data/Galaxy42-[eggNOG Mapper on data 40_ annotations].tabular", sample_name="40")
```

Esto generará una figura PNG en la carpeta `figures/` con las anotaciones funcionales representadas como paneles de barras.

## Autores
- Dra. Katia Aviña Padilla
- Proyecto Viroma – Biocomplejidad, Cuatro Ciénegas

## Licencia
Este repositorio está disponible bajo una licencia abierta. Puedes reutilizarlo citando a los autores del proyecto.
