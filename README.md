# MaRS-py-upgrade 
**Réimplémentation et extension en Python du pipeline MaRS pour l’analyse des marqueurs moléculaires de résistance de *Plasmodium falciparum***

---

## Présentation générale
**MaRS-py-upgrade** est un pipeline bioinformatique modulaire développé en **Python**, destiné à l’analyse des données de séquençage NGS de *Plasmodium falciparum* afin d’identifier et de caractériser les marqueurs moléculaires associés à la résistance aux antipaludiques.

Ce pipeline s’inscrit dans un cadre académique et de recherche, notamment pour l’analyse des gènes **pfcrt**, **pfmdr1**,**pfk13**, **pfdhfr** et **pfdhps**, utilisés comme marqueurs de résistance aux traitements antipaludiques.

Il s’agit d’une réimplémentation et d’une extension du pipeline **MaRS**, avec une architecture lisible, reproductible et automatisée.

---

## Objectifs
- Automatiser l’analyse bioinformatique des données NGS de *Plasmodium falciparum* 
- Identifier les variants génétiques associés à la résistance aux antipaludiques en utilisant plusieurs outils d’appel de variants (samtools, GATK, freebayes et vardict)
- Calcul des VAF (fréquence allélique du variant)  par gène et par site
- Analyser les haplotypes par gène et par site
- Générer de rapport de synthèse et de visualisation exploitable
- Garantir la traçabilité des analyses via des fichiers de logs

---

## Données analysées
- Données de séquençage NGS (FASTQ compressés) utiisant la nommanclature **AMD_ID**
- Échantillons individuels et/ou poolés
- Génome de référence : *Plasmodium falciparum* 3D7

---

## Workflow général
Le pipeline est structuré sous forme de modules fonctionnels indépendants, exécutés de manière séquentielle :

1. Contrôle des données FASTQ  
2. Trimming 
3. Alignement des lectures sur le génome de référence (*Pf3D7*)
4. Appel de variants avec plusieurs outils :
   - Samtools
   - FreeBayes
   - GATK HaplotypeCaller
   - VarDict
5. Fusion et harmonisation des fichiers VCF  
6. Filtrage et annotation des variants
7. Calcul de VAF (fréquence allélique du variant) 
8. Analyse des haplotypes par gène  
9. Génération de rapports et de graphiques

---

## 📁 Organisation du projet

```text
HOME/
└── pipeline/
    ├── data/                     # Données brutes (FASTQ)
    │   └── *.fastq.gz
    │
    ├── bin/                      # Scripts secondaires appelés dans le script principale
    │   
    │
    ├── output/                   # Résultats générés
    │   ├── QC/
    │   ├── bam/
    │   ├── variants/
    │   └── haplotypes/, etc.
    │
    ├── logs/                     # Logs d’exécution
    │   └── *.log
    │
    ├── pf_3D7/                   # Génome de référence
    │
    ├── pf_3D7_snpEff_db          # Création du base d'annotation
    |
    |
    ├── pipeline_python.py        # Script du pipeline
    ├── Dockerfile                # Fichier docker
    └── README.md                 # Fichier README

```
---

## Prérequis
- Système d'exploitation de type Unix (Linux, macOS, etc.) 
- Docker

---

## Installation
### 1. Installer Docker
Pour installer Docker, executer succesivement les codes suivants dans le terminale :
```bash

sudo apt update

sudo apt install -y ca-certificates curl gnupg

sudo install -m 0755 -d /etc/apt/keyrings

curl -fsSL https://download.docker.com/linux/ubuntu/gpg | \ 
sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg

sudo chmod a+r /etc/apt/keyrings/docker.gpg

echo \
"deb [arch=$(dpkg --print-architecture) \
signed-by=/etc/apt/keyrings/docker.gpg] \
https://download.docker.com/linux/ubuntu \
$(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

sudo apt update

sudo apt install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

sudo usermod -aG docker $USER

newgrp docker
```
### 2. Télécharger le dépôt git

```bash

mkdir -p ~/pipeline

git clone 
