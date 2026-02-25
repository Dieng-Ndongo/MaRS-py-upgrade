"""
#!/usr/bin/env python3
import matplotlib
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import argparse

matplotlib.use('Agg')
parser = argparse.ArgumentParser(description='filename')
parser.add_argument('-n', dest='filename', type=str, help="name of snp summary file")

args = parser.parse_args()
filename=args.filename

try:
    Novel = pd.read_csv(filename)
except pd.errors.EmptyDataError:
    print('Empty Novel snps csv file!')
    sys.exit

df = Novel.groupby(['CHROM','VOI','Type','Annotation']).size().reset_index(name='counts')

df_pv= df.pivot_table(values='counts', index=['CHROM','VOI','Annotation'], columns='Type', aggfunc='first')
df_pv = df_pv.fillna(0).reset_index()

column_names = ['Mixed','Mutant']
df_pv['Total']= df_pv[column_names].sum(axis=1)
df_pv["Snps"] = df_pv["CHROM"] + ":" + df_pv["VOI"] + ":N=" + df_pv["Total"].astype(str)

df_pv = df_pv.rename(columns={'Mixed': 'Minor', 'Mutant':'Major'})
SNPvals=df_pv[["Snps",'Minor','Major','Total','Annotation']]

# Separate synonymous and missense:
SNPs_NS  = SNPvals[SNPvals['Annotation']  == "missense_variant"]


######################## Novel missense SNPS graph

#Setup for loading
Totes = SNPs_NS.groupby('Snps')['Total'].sum().reset_index()
Minor = SNPs_NS.groupby('Snps')['Minor'].sum().reset_index()
Major = SNPs_NS.groupby('Snps')['Major'].sum().reset_index()

#Math and definition of SNPratio
Minor['SNPratio'] = [i / j for i,j in zip(Minor['Minor'], Totes['Total'])]
Major['SNPratio'] = [i / j for i,j in zip(Major['Major'], Totes['Total'])]

AllTogether = pd.concat([Minor.Snps, Minor.SNPratio, Major.SNPratio], axis=1)

AllTogether.columns = ['Snps','Minor: AF < 95%', 'Major: AF >= 95%']

df_table_SNP=AllTogether.sort_values(by=['Snps'])

df_table_SNP["index"]=df_table_SNP.Snps.str.split(":").str[1].str[1:-1]
df_table_SNP["index"] = df_table_SNP["index"].str.extract('(\d+)').astype(int)
df_table_SNP["index2"]=df_table_SNP.Snps.str.split(":").str[0]

plot = df_table_SNP.sort_values(by = ['index2', 'index'],ascending=False)[['Snps',  'Minor: AF < 95%','Major: AF >= 95%']].plot(x='Snps', kind='barh', stacked=True, title='Novel missense Mutations', figsize=(20,20), color={"Minor: AF < 95%": "#F3ABA8", "Major: AF >= 95%": "#98DAA7"})

plot.legend(ncol = 2, loc = 'lower right')
plot.set(ylabel="SNPs")
plot.set(xlabel="SNP ratio")
plot.legend(loc=(1,0))
plt.savefig('SNPs-Novel-missense.pdf')




######################## Novel Synonymous SNPS graph

SNPs_S  = SNPvals[SNPvals['Annotation']  == "synonymous_variant"]

#Setup for loading
Totes = SNPs_S.groupby('Snps')['Total'].sum().reset_index()
Minor = SNPs_S.groupby('Snps')['Minor'].sum().reset_index()
Major = SNPs_S.groupby('Snps')['Major'].sum().reset_index()

#Math and definition of SNPratio
Minor['SNPratio'] = [i / j for i,j in zip(Minor['Minor'], Totes['Total'])]
Major['SNPratio'] = [i / j for i,j in zip(Major['Major'], Totes['Total'])]

AllTogether = pd.concat([Minor.Snps, Minor.SNPratio, Major.SNPratio], axis=1)

AllTogether.columns = ['Snps','Minor: AF < 95%', 'Major: AF >= 95%']

df_table_SNP=AllTogether.sort_values(by=['Snps'])

df_table_SNP["index"]=df_table_SNP.Snps.str.split(":").str[1].str[1:-1]
df_table_SNP["index"] = df_table_SNP["index"].str.extract('(\d+)').astype(int)
df_table_SNP["index2"]=df_table_SNP.Snps.str.split(":").str[0]

plot = df_table_SNP.sort_values(by = ['index2', 'index'],ascending=False)[['Snps', 'Minor: AF < 95%','Major: AF >= 95%']].plot(x='Snps', kind='barh', stacked=True, title='Novel synonymous Mutations', figsize=(20,20), color={"Minor: AF < 95%": "#F3ABA8", "Major: AF >= 95%": "#98DAA7"})
#plot.legend(bbox_to_anchor=(0.97, 0.1))

plot.legend(ncol = 2, loc = 'lower right')
#sns.despine(left = True, bottom = True)
plot.set(ylabel="SNPs")
plot.set(xlabel="SNP ratio")
plot.legend(loc=(1,0))
plt.savefig('SNPs-Novel-synonymous.pdf')

"""
#!/usr/bin/env python3
import matplotlib
matplotlib.use('Agg')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys


# =========================
# Arguments
# =========================
parser = argparse.ArgumentParser(description='Novel SNPs dataviz by site')
parser.add_argument('-n', dest='filename', type=str, required=True,
                    help="Novel SNP summary CSV file")
args = parser.parse_args()


# =========================
# Lecture du fichier
# =========================
try:
    Novel = pd.read_csv(args.filename)
except pd.errors.EmptyDataError:
    print("Empty Novel SNPs CSV file!")
    sys.exit(1)


# =========================
# Extraction du site
# =========================
# Site = caractères [4:6] du nom de l’échantillon
Novel["Site"] = Novel["Sample_name"].astype(str).str[4:6]


# =========================
# Résumé SNPs
# =========================
df = (
    Novel
    .groupby(['Site', 'CHROM', 'VOI', 'Type', 'Annotation'])
    .size()
    .reset_index(name='counts')
)

df_pv = (
    df.pivot_table(
        values='counts',
        index=['Site', 'CHROM', 'VOI', 'Annotation'],
        columns='Type',
        aggfunc='first'
    )
    .fillna(0)
    .reset_index()
)

# Renommage
df_pv = df_pv.rename(columns={'Mixed': 'Minor', 'Mutant': 'Major'})

# Total
df_pv['Total'] = df_pv[['Minor', 'Major']].sum(axis=1)

# Label SNP
df_pv['Snps'] = (
    df_pv['CHROM'] + ":" +
    df_pv['VOI'] + ":N=" +
    df_pv['Total'].astype(int).astype(str)
)

SNPvals = df_pv[['Site', 'Snps', 'Minor', 'Major', 'Total', 'Annotation']]


# =========================
# Fonction utilitaire
# =========================
def prepare_plot_dataframe(df):
    # Prépare le dataframe pour affichage par site
    Totes = df.groupby(['Site', 'Snps'])['Total'].sum().reset_index()
    Minor = df.groupby(['Site', 'Snps'])['Minor'].sum().reset_index()
    Major = df.groupby(['Site', 'Snps'])['Major'].sum().reset_index()

    Minor['SNPratio'] = Minor['Minor'] / Totes['Total']
    Major['SNPratio'] = Major['Major'] / Totes['Total']

    plot_df = pd.concat(
        [Minor[['Site', 'Snps']], Minor['SNPratio'], Major['SNPratio']],
        axis=1
    )

    plot_df.columns = ['Site', 'Snps', 'Minor: AF < 95%', 'Major: AF >= 95%']

    # Extraction position SNP pour tri
    plot_df['index'] = (
        plot_df['Snps']
        .str.split(":").str[1]
        .str.extract(r'(\d+)')
        .astype(int)
    )
    plot_df['index2'] = plot_df['Snps'].str.split(":").str[0]

    plot_df = plot_df.sort_values(
        by=['Site', 'index2', 'index'],
        ascending=[True, False, False]
    )

    # Création labels Y avec affichage du site
    ylabels = []
    last_site = None
    for _, row in plot_df.iterrows():
        if row['Site'] != last_site:
            ylabels.append(f"SITE {row['Site']} | {row['Snps']}")
            last_site = row['Site']
        else:
            ylabels.append(f"      {row['Snps']}")

    plot_df['Ylabel'] = ylabels

    return plot_df


def add_site_separators(ax, df):
    #Ajoute une ligne horizontale entre chaque site
    sites = df['Site'].values
    last_site = sites[0]
    for i, site in enumerate(sites):
        if site != last_site:
            ax.axhline(i - 0.5, color='black', linewidth=1)
            last_site = site


# =========================
# MISSENSE
# =========================
SNPs_NS = SNPvals[SNPvals['Annotation'] == "missense_variant"]

if not SNPs_NS.empty:
    df_plot = prepare_plot_dataframe(SNPs_NS)

    ax = df_plot[
        ['Ylabel', 'Minor: AF < 95%', 'Major: AF >= 95%']
    ].plot(
        x='Ylabel',
        kind='barh',
        stacked=True,
        figsize=(20, 20),
        title='Novel missense Mutations (by site)',
        color={
            "Minor: AF < 95%": "#F3ABA8",
            "Major: AF >= 95%": "#98DAA7"
        }
    )

    add_site_separators(ax, df_plot)
    
    # Mettre les noms des sites en gras sur l’axe Y
    for label in ax.get_yticklabels():
        if label.get_text().startswith("SITE"):
            label.set_fontweight("bold")
    
    ax.legend(ncol = 2, loc = 'lower right')
    ax.set_ylabel("SNPs")
    ax.set_xlabel("SNP ratio")
    ax.legend(loc=(1, 0))


    plt.tight_layout()
    plt.savefig('SNPs-Novel-missense.pdf')
    plt.close()


# =========================
# SYNONYMOUS
# =========================
SNPs_S = SNPvals[SNPvals['Annotation'] == "synonymous_variant"]

if not SNPs_S.empty:
    df_plot = prepare_plot_dataframe(SNPs_S)

    ax = df_plot[
        ['Ylabel', 'Minor: AF < 95%', 'Major: AF >= 95%']
    ].plot(
        x='Ylabel',
        kind='barh',
        stacked=True,
        figsize=(20, 20),
        title='Novel synonymous Mutations (by site)',
        color={
            "Minor: AF < 95%": "#F3ABA8",
            "Major: AF >= 95%": "#98DAA7"
        }
    )

    add_site_separators(ax, df_plot)
    
     # Mettre les noms des sites en gras sur l’axe Y
    for label in ax.get_yticklabels():
        if label.get_text().startswith("SITE"):
            label.set_fontweight("bold")

    ax.legend(ncol = 2, loc = 'lower right')
    ax.set_ylabel("SNPs")
    ax.set_xlabel("SNP ratio")
    ax.legend(loc=(1, 0))
    
    plt.tight_layout()
    plt.savefig('SNPs-Novel-synonymous.pdf')
    plt.close()

