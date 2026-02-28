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

# ─── AJOUT : extraction du site depuis Sample_name[4:6] ──────────────────────
Novel['Site'] = Novel['Sample_name'].str[4:6]
# ─────────────────────────────────────────────────────────────────────────────

df = Novel.groupby(['CHROM','VOI','Type','Annotation']).size().reset_index(name='counts')

df_pv= df.pivot_table(values='counts', index=['CHROM','VOI','Annotation'], columns='Type', aggfunc='first')
df_pv = df_pv.fillna(0).reset_index()

column_names = ['Mixed','Mutant']
df_pv['Total']= df_pv[column_names].sum(axis=1)
df_pv["Snps"] = df_pv["CHROM"] + ":" + df_pv["VOI"] + ":N=" + df_pv["Total"].astype(str)

df_pv = df_pv.rename(columns={'Mixed': 'Minor', 'Mutant':'Major'})
SNPvals=df_pv[["Snps",'Minor','Major','Total','Annotation']]

# ─── AJOUT : table SNP → Site (via Novel) ────────────────────────────────────
Snps_site = Novel.groupby(['CHROM','VOI'])['Site'].first().reset_index()
Snps_site["Snps_key"] = Snps_site["CHROM"] + ":" + Snps_site["VOI"]
# ─────────────────────────────────────────────────────────────────────────────

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

# ─── AJOUT : jointure site + tri par site puis position ──────────────────────
df_table_SNP["Snps_key"] = df_table_SNP["index2"] + ":" + df_table_SNP.Snps.str.split(":").str[1]
df_table_SNP = df_table_SNP.merge(Snps_site[["Snps_key","Site"]], on="Snps_key", how="left")
df_table_SNP_sorted = df_table_SNP.sort_values(by=['Site', 'index2', 'index'], ascending=[True, True, False]).reset_index(drop=True)
# ─────────────────────────────────────────────────────────────────────────────

plot = df_table_SNP_sorted[['Snps',  'Minor: AF < 95%','Major: AF >= 95%']].plot(x='Snps', kind='barh', stacked=True, title='Novel missense Mutations', figsize=(20,20), color={"Minor: AF < 95%": "#F3ABA8", "Major: AF >= 95%": "#98DAA7"})

plot.legend(ncol = 2, loc = 'lower right')
plot.set(ylabel="SNPs")
plot.set(xlabel="SNP ratio")
plot.legend(loc=(1,0))

# ─── AJOUT : lignes de séparation et nom unique par site ─────────────────────
site_per_bar = df_table_SNP_sorted["Site"].values

current_site = None
site_start   = 0
for i, site in enumerate(site_per_bar):
    if site != current_site:
        if current_site is not None:
            # Ligne de séparation
            plot.axhline(y=i - 0.5, color='black', linewidth=1.5, linestyle='--', alpha=0.7)
            # Nom du site centré sur son groupe
            plot.text(1.01, (site_start + i - 1) / 2, f"Site: {current_site}",
                      transform=plot.get_yaxis_transform(),
                      va='center', ha='left', fontsize=9, fontweight='bold')
        current_site = site
        site_start   = i

# Dernier site
plot.text(1.01, (site_start + len(site_per_bar) - 1) / 2, f"Site: {current_site}",
          transform=plot.get_yaxis_transform(),
          va='center', ha='left', fontsize=9, fontweight='bold')
# ─────────────────────────────────────────────────────────────────────────────

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

# ─── AJOUT : jointure site + tri par site puis position ──────────────────────
df_table_SNP["Snps_key"] = df_table_SNP["index2"] + ":" + df_table_SNP.Snps.str.split(":").str[1]
df_table_SNP = df_table_SNP.merge(Snps_site[["Snps_key","Site"]], on="Snps_key", how="left")
df_table_SNP_sorted = df_table_SNP.sort_values(by=['Site', 'index2', 'index'], ascending=[True, True, False]).reset_index(drop=True)
# ─────────────────────────────────────────────────────────────────────────────

plot = df_table_SNP_sorted[['Snps', 'Minor: AF < 95%','Major: AF >= 95%']].plot(x='Snps', kind='barh', stacked=True, title='Novel synonymous Mutations', figsize=(20,20), color={"Minor: AF < 95%": "#F3ABA8", "Major: AF >= 95%": "#98DAA7"})

plot.legend(ncol = 2, loc = 'lower right')
plot.set(ylabel="SNPs")
plot.set(xlabel="SNP ratio")
plot.legend(loc=(1,0))

# ─── AJOUT : lignes de séparation et nom unique par site ─────────────────────
site_per_bar = df_table_SNP_sorted["Site"].values

current_site = None
site_start   = 0
for i, site in enumerate(site_per_bar):
    if site != current_site:
        if current_site is not None:
            plot.axhline(y=i - 0.5, color='black', linewidth=1.5, linestyle='--', alpha=0.7)
            plot.text(1.01, (site_start + i - 1) / 2, f"Site: {current_site}",
                      transform=plot.get_yaxis_transform(),
                      va='center', ha='left', fontsize=9, fontweight='bold')
        current_site = site
        site_start   = i

# Dernier site
plot.text(1.01, (site_start + len(site_per_bar) - 1) / 2, f"Site: {current_site}",
          transform=plot.get_yaxis_transform(),
          va='center', ha='left', fontsize=9, fontweight='bold')
# ─────────────────────────────────────────────────────────────────────────────

plt.savefig('SNPs-Novel-synonymous.pdf')
