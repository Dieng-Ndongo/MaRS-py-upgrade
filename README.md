# MaRS-py-upgrade 
**Réimplémentation et extension du pipeline MaRS en Python pour l’analyse des marqueurs moléculaires de résistance de *Plasmodium falciparum***

---

## Présentation générale
**MaRS-py-upgrade** est un pipeline bioinformatique modulaire développé en **Python**, destiné à l’analyse des données de séquençage NGS de *Plasmodium falciparum* afin d’identifier et de caractériser les marqueurs moléculaires associés à la résistance aux antipaludiques.

Ce pipeline s’inscrit dans un cadre académique et de recherche, notamment pour l’analyse des gènes ***pfcrt***, ***pfmdr1***,***pfk13***, ***pfdhfr*** et ***pfdhps***, utilisés comme marqueurs de résistance aux traitements antipaludiques.

Il s’agit d’une réimplémentation et d’une extension du pipeline **MaRS**, avec une architecture lisible, reproductible et automatisée. 
Le pipeline peux etre utiliser via une interface graphique.

---

## Objectifs
- Automatiser l’analyse bioinformatique des données NGS de *Plasmodium falciparum* 
- Identifier les variants génétiques associés à la résistance aux antipaludiques en utilisant plusieurs outils d’appel de variants (samtools, GATK, freebayes et vardict)
- Calcul des VAF (fréquence allélique du variant)  par gène et par site
- Analyser les haplotypes par gène et par site
- Générer de rapport de synthèse et de visualisation exploitable
- Garantir la traçabilité des analyses via des fichiers de logs
- Rendre le pipeline MaRS accessible et facile d’utilisation

---

## Données
- Données de séquençage NGS (fastq compressés)

### Nomenclature des echantillons

![Nomenclature AMD_ID pour les échantillons individuels](images/individual_AMD_ID.png)
**Exemple :** 23SNKG28A0004PfB4721_KG48D28_S8_L001_R1_001.fastq.gz 

![Nomenclature AMD_ID pour les échantillons poolés](images/pooled_AMD_ID.png)
**Exemple :** 23SNKG00I033P10B4721_POOL33_S51_L001_R1_001.fastq.gz

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
└── MaRS-py-upgrade/              # Dossier principal
    |          
    ├── data/                     # Données brutes (FASTQ)
    │   └── *.fastq.gz
    │
    ├── bin/                      # Scripts secondaires appelés dans le script principale
    │   
    │
    ├── output/                   # Résultats générés
    │
    │
    ├── logs/                     # Logs d’exécution
    │   └── *.log
    │
    ├── pf_3D7/                   # Génome de référence
    │
    ├── pf_3D7_snpEff_db          # Création du base d'annotation
    |
    |  
    ├── images                    # dossier contenant des images
    |    
    |
    ├── pipeline_python.py        # Script du pipeline  
    ├── environment.yml           # environment       
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
Pour installer Docker, executer succesivement les commandes suivantes dans le terminale :
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

git clone https://github.com/Dieng-Ndongo/MaRS-py-upgrade.git

cd MaRS-py-upgrade

```
### 3. Construction de l'image docker

```bash

docker build -t bioinfo_pipeline .

```
### 4. Exécution du pipeline

```bash

docker run --rm -it -v $(pwd):/app -w /app bioinfo_pipeline

```
