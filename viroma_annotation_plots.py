
# viroma_annotation_plots.py

import pandas as pd
import matplotlib.pyplot as plt
import os

# Diccionario de descripciones COG
COG_DESCRIPTIONS = {
    'A': 'RNA processing and modification',
    'B': 'Chromatin structure and dynamics',
    'C': 'Energy production and conversion',
    'D': 'Cell cycle control, cell division, chromosome partitioning',
    'E': 'Amino acid transport and metabolism',
    'F': 'Nucleotide transport and metabolism',
    'G': 'Carbohydrate transport and metabolism',
    'H': 'Coenzyme transport and metabolism',
    'I': 'Lipid transport and metabolism',
    'J': 'Translation, ribosomal structure and biogenesis',
    'K': 'Transcription',
    'L': 'Replication, recombination and repair',
    'M': 'Cell wall/membrane/envelope biogenesis',
    'N': 'Cell motility',
    'O': 'Posttranslational modification, protein turnover, chaperones',
    'P': 'Inorganic ion transport and metabolism',
    'Q': 'Secondary metabolites biosynthesis, transport and catabolism',
    'R': 'General function prediction only',
    'S': 'Function unknown',
    'T': 'Signal transduction mechanisms',
    'U': 'Intracellular trafficking, secretion, and vesicular transport',
    'V': 'Defense mechanisms',
    'W': 'Extracellular structures',
    'Y': 'Nuclear structure',
    'Z': 'Cytoskeleton'
}

def count_terms(series):
    return pd.Series(
        term.strip() for row in series if row != '-' for term in str(row).split(',')
    ).value_counts().head(15)

def generate_panel_plot(file_path, sample_name, output_dir="figures"):
    os.makedirs(output_dir, exist_ok=True)

    df = pd.read_csv(file_path, sep="\t", comment="#", low_memory=False)
    cog_col, ko_col, kegg_col, pfam_col = df.iloc[:, 6], df.iloc[:, 11], df.iloc[:, 12], df.iloc[:, 20]

    top_cog = cog_col.value_counts().head(15)
    top_kos = count_terms(ko_col.dropna())
    top_kegg = count_terms(kegg_col.dropna())
    top_pfam = count_terms(pfam_col.dropna())

    fig, axs = plt.subplots(2, 2, figsize=(18, 14), dpi=300)

    # Panel 1 - COG
    axs[0, 0].barh(top_cog.index[::-1], top_cog.values[::-1], color='steelblue')
    axs[0, 0].set_title(f'Top 15 COG Categories in Sample {sample_name}', fontsize=12, weight='bold')
    axs[0, 0].set_xlabel('Count')
    axs[0, 0].tick_params(labelsize=8)
    cog_desc_labels = [COG_DESCRIPTIONS.get(c, 'Multiple') for c in top_cog.index[::-1]]
    for i, (label, count) in enumerate(zip(cog_desc_labels, top_cog.values[::-1])):
        axs[0, 0].text(count + 3, i, label, va='center', fontsize=8, color='black')

    # Panel 2 - KO
    axs[0, 1].barh(top_kos.index[::-1], top_kos.values[::-1], color='lightgreen')
    axs[0, 1].set_title('Top 15 KEGG Orthologs (KO)', fontsize=12, weight='bold')
    axs[0, 1].set_xlabel('Count')
    axs[0, 1].tick_params(labelsize=8)

    # Panel 3 - KEGG Pathways
    axs[1, 0].barh(top_kegg.index[::-1], top_kegg.values[::-1], color='salmon')
    axs[1, 0].set_title('Top 15 KEGG Pathways', fontsize=12, weight='bold')
    axs[1, 0].set_xlabel('Count')
    axs[1, 0].tick_params(labelsize=8)

    # Panel 4 - Pfam Domains
    axs[1, 1].barh(top_pfam.index[::-1], top_pfam.values[::-1], color='plum')
    axs[1, 1].set_title('Top 15 Pfam Domains', fontsize=12, weight='bold')
    axs[1, 1].set_xlabel('Count')
    axs[1, 1].tick_params(labelsize=8)

    plt.tight_layout()
    output_path = os.path.join(output_dir, f"Functional_Profile_Sample{sample_name}_MultiPanel_COG_Labeled.png")
    plt.savefig(output_path, format='png', dpi=300)
    plt.close()
    print(f"Saved: {output_path}")
