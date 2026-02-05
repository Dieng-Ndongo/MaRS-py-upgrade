# 🧬 MaRS-Py-upgrade 
**Réimplémentation et extension en Python du pipeline MaRS pour l’analyse des marqueurs moléculaires de résistance de *Plasmodium falciparum***

---

## 📌 Présentation générale
**MaRS-Py-upgrade** est un pipeline bioinformatique modulaire développé en **Python**, destiné à l’analyse des données de séquençage NGS de *Plasmodium falciparum* afin d’identifier et de caractériser les marqueurs moléculaires associés à la résistance aux antipaludiques.

Ce pipeline s’inscrit dans un cadre académique et de recherche, notamment pour l’analyse des gènes **pfcrt**, **pfmdr1**,**pfk13**, **pfdhfr** et **pfdhps**, utilisés comme marqueurs de résistance aux traitements antipaludiques.

Il s’agit d’une réimplémentation et d’une extension du pipeline **MaRS**, avec une architecture plus lisible, reproductible et automatisée.

---

## 🎯 Objectifs
- Automatiser l’analyse bioinformatique des données NGS
- Identifier les variants génétiques associés à la résistance aux antipaludiques
- Comparer les résultats issus de plusieurs outils d’appel de variants
- Analyser les haplotypes par gène et par site
- Générer des rapports de synthèse et des visualisations exploitables
- Garantir la traçabilité des analyses via des fichiers de logs

---

## 🔬 Données analysées
- Données de séquençage NGS (FASTQ compressés)
- Échantillons individuels et/ou poolés
- Génome de référence : *Plasmodium falciparum* 3D7

---

## 🔄 Workflow général
Le pipeline est structuré sous forme de modules fonctionnels indépendants, exécutés de manière séquentielle :

1. Préparation et contrôle des données FASTQ  
2. Alignement des lectures sur le génome de référence (*Pf3D7*)  
3. Traitement des fichiers BAM  
4. Appel de variants avec plusieurs outils :
   - Samtools
   - FreeBayes
   - GATK HaplotypeCaller
   - VarDict
5. Fusion et harmonisation des fichiers VCF  
6. Filtrage et annotation des variants  
7. Analyse des haplotypes par gène  
8. Génération de rapports et de graphiques

---

## 📁 Organisation du projet

```text
HOME/
└── pipeline/
    ├── data/                     # Données brutes (FASTQ)
    │   └── *.fastq.gz
    │
    ├── bin/                      # Scripts Python du pipeline
    │   ├── fastq_processing.py
    │   ├── alignment.py
    │   ├── variant_calling.py
    │   ├── csv_merge.py
    │   └── reporting.py
    │
    ├── output/                   # Résultats générés
    │   ├── fastq/
    │   ├── bam/
    │   ├── variants/
    │   ├── haplotypes/
    │   └── reports/
    │
    ├── logs/                     # Logs d’exécution
    │   └── *.log
    │
    ├── pf_3D7/                   # Génome de référence
    │
    ├── requirements.txt
    └── README.md
---

## ⚙️ Prérequis

- Linux (recommandé)

- macOS (non testé)
