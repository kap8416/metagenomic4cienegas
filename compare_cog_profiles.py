
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import os

COG_DESCRIPTIONS = {
    'A': 'RNA processing and modification', 'B': 'Chromatin structure and dynamics',
    'C': 'Energy production and conversion', 'D': 'Cell cycle control and chromosome partitioning',
    'E': 'Amino acid transport and metabolism', 'F': 'Nucleotide transport and metabolism',
    'G': 'Carbohydrate transport and metabolism', 'H': 'Coenzyme transport and metabolism',
    'I': 'Lipid transport and metabolism', 'J': 'Translation, ribosomal structure and biogenesis',
    'K': 'Transcription', 'L': 'Replication, recombination and repair',
    'M': 'Cell wall/membrane/envelope biogenesis', 'N': 'Cell motility',
    'O': 'Posttranslational modification, protein turnover, chaperones',
    'P': 'Inorganic ion transport and metabolism', 'Q': 'Secondary metabolites biosynthesis',
    'R': 'General function prediction only', 'S': 'Function unknown',
    'T': 'Signal transduction mechanisms', 'U': 'Intracellular trafficking, secretion, and vesicular transport',
    'V': 'Defense mechanisms', 'W': 'Extracellular structures', 'Y': 'Nuclear structure', 'Z': 'Cytoskeleton'
}

def load_cog_counts(file_list, sample_names):
    cog_matrix = {}
    for file_path, sample_name in zip(file_list, sample_names):
        df = pd.read_csv(file_path, sep="\t", comment="#", low_memory=False)
        cog_counts = df.iloc[:, 6].value_counts()
        for cog, count in cog_counts.items():
            if cog not in cog_matrix:
                cog_matrix[cog] = {}
            cog_matrix[cog][sample_name] = count
    return pd.DataFrame(cog_matrix).T.fillna(0).astype(int)

def plot_cog_stacked(df, title, output_prefix):
    prop_df = df.div(df.sum(axis=0), axis=1)
    top_cogs = df.sum(axis=1).sort_values(ascending=False).head(10).index
    plot_df = prop_df.loc[top_cogs].T
    fig, ax = plt.subplots(figsize=(20, 12), dpi=300)
    plot_df.plot(kind='bar', stacked=True, ax=ax, colormap='tab20')
    ax.set_title(title, fontsize=16, weight='bold')
    ax.set_ylabel('Proportion', fontsize=12)
    ax.set_xlabel('Sample', fontsize=12)
    ax.set_xticklabels(plot_df.index, rotation=0, ha='center', fontsize=11)
    ax.legend(title='COG Category', bbox_to_anchor=(1.28, 1), loc='upper left', fontsize=10, title_fontsize=11)

    legend_text = "\n".join([f"{k}: {v}" for k, v in COG_DESCRIPTIONS.items() if k in plot_df.columns])
    props = dict(boxstyle='round', facecolor='white', alpha=0.95, edgecolor='gray')
    plt.gcf().text(1.32, 0.5, legend_text, fontsize=8, verticalalignment='center', bbox=props)
    plt.subplots_adjust(left=0.08, right=0.75)

    fig.savefig(f"{output_prefix}.png", bbox_inches='tight')
    fig.savefig(f"{output_prefix}.pdf", bbox_inches='tight')
    fig.savefig(f"{output_prefix}.svg", bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    os.makedirs("figures", exist_ok=True)

    wells = [
        ("data/Galaxy10-[eggNOG Mapper on data 8_ annotations].tabular", "Well 1"),
        ("data/Galaxy14-[eggNOG Mapper on data 12_ annotations].tabular", "Well 2"),
        ("data/Galaxy18-[eggNOG Mapper on data 16_ annotations].tabular", "Well 3"),
        ("data/Galaxy22-[eggNOG Mapper on data 20_ annotations].tabular", "Well 4"),
        ("data/Galaxy26-[eggNOG Mapper on data 24_ annotations].tabular", "Well 5"),
        ("data/Galaxy30-[eggNOG Mapper on data 28_ annotations].tabular", "Well 6"),
    ]
    stroms = [
        ("data/Galaxy34-[eggNOG Mapper on data 32_ annotations].tabular", "Stromatolite 1"),
        ("data/Galaxy38-[eggNOG Mapper on data 36_ annotations].tabular", "Stromatolite 2"),
        ("data/Galaxy42-[eggNOG Mapper on data 40_ annotations].tabular", "Stromatolite 3"),
        ("data/Galaxy46-[eggNOG Mapper on data 44_ annotations].tabular", "Stromatolite 4"),
    ]

    well_df = load_cog_counts([f[0] for f in wells], [f[1] for f in wells])
    strom_df = load_cog_counts([f[0] for f in stroms], [f[1] for f in stroms])

    plot_cog_stacked(well_df, "Relative Abundance of Top COG Categories in Well Samples", "figures/Well_COG_Profile")
    plot_cog_stacked(strom_df, "Relative Abundance of Top COG Categories in Stromatolite Samples", "figures/Stromatolite_COG_Profile")
