# MaRS-py-upgrade 
**Réimplémentation et extension du pipeline bioinformatique MaRS en language Python pour l’analyse des marqueurs moléculaires de résistance de *Plasmodium falciparum***

---

## Présentation générale
**MaRS-py-upgrade** est un pipeline bioinformatique modulaire développé en **Python** avec interface web, destiné à l’analyse des données de séquençage NGS de *Plasmodium falciparum* afin d’identifier et de caractériser les marqueurs moléculaires associés à la résistance aux antipaludiques.

Ce pipeline s’inscrit dans un cadre académique et de recherche, notamment pour l’analyse des gènes ***pfcrt***, ***pfmdr1***,***pfk13***, ***pfdhfr*** et ***pfdhps***, utilisés comme marqueurs de résistance aux traitements antipaludiques.

Il s’agit d’une réimplémentation et d’une extension du pipeline **MaRS**, avec une architecture lisible, reproductible et automatisée. 
Le pipeline est utiliser via une interface web.

---

## Objectifs
- Automatiser l’analyse bioinformatique des données NGS de *Plasmodium falciparum* 
- Identifier les variants génétiques associés à la résistance aux antipaludiques en utilisant plusieurs outils d’appel de variants (samtools, GATK, freebayes et vardict)
- Calcul des VAF (fréquence allélique du variant)  par gène et par site
- Analyser les haplotypes par gène et par site
- Générer de rapport de synthèse et de visualisation exploitable
- Garantir la traçabilité des analyses via des fichiers de logs
- Rendre le pipeline MaRS accessible et facile d’utilisation via une interface web

---

## Données
- Données de séquençage NGS (fastq compressés(.fastq.gz))

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

## 📁 Structure du projet

```text
MaRS-py-upgrade/
│
├── 📂 bin/        
│      └── ...                  # Scripts secondaires
|
├── 📂 pf_3D7/
│   └── ...                     # Génome de référence de Plasmodium falciparum
│
├── 📂 pf_3D7_snpEff_db/
│   └── ...                     # Base de données d'annotation SnpEff
│
├── 📂 images/                  # Images utilisées dans l'interface
│         └── ...           
│
├── 📂 runs/
│   └── run_id/
│       ├── 📂 data/
│       │   └── *.fastq.gz      # Données de séquençage d'entrée
│       │
│       ├── 📂 output/
│       │   └── ...             # Résultats de l'analyse
│       │
│       └── 📂 log/
│           └── ...             # Journaux d'exécution
│
│
├── 🐍 pipeline_python.py       # Script principal du pipeline
├── 🖥️ app.py                   # Interface utilisateur Streamlit
├── 🐳 Dockerfile               # Configuration de l'image Docker
├── 📦 environment.yml          # Environnement et dépendances Conda du pipeline
├── 📄 requirements.txt         # Dépendances Python de l'application web
├── 🚀 start.sh                 # Script de lancement du pipeline
├── 📄 .dockerignore            # Fichiers exclus de l'image Docker
├── 📄 .gitignore               # Fichiers exclus du dépôt Git
└── 📖 README.md                # Documentation du projet

```markdown
## 📌 Description des principaux composants

| Élément | Description |
|---|---|
| `pipeline_python.py` | Implémentation principale du pipeline bioinformatique en Python. |
| `app.py` | Interface graphique permettant de lancer et suivre les analyses avec Streamlit. |
| `bin/` | Regroupe les scripts secondaires utilisés par le pipeline et l'interface. |
| `pf_3D7/` | Contient le génome de référence de *Plasmodium falciparum* 3D7. |
| `pf_3D7_snpEff_db/` | Contient les fichiers nécessaires à l'annotation des variants avec SnpEff. |
| `runs/` | Regroupe les résultats générés pour chaque analyse. |
| `environment.yml` | Définit l'environnement Conda et les dépendances du pipeline. |
| `requirements.txt` | Liste les dépendances Python nécessaires à l'interface. |
| `Dockerfile` | Permet de construire l'environnement d'exécution du pipeline avec Docker. |
| `start.sh` | Automatise le lancement du pipeline et/ou de l'interface. |
| `.gitignore` | Définit les fichiers et répertoires qui ne doivent pas être suivis par Git. |

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

git clone https://github.com/Dieng-Ndongo/MaRS-py-upgrade

cd MaRS-py-upgrade

```
### 3. Construction de l'image docker

```bash

docker build -t bioinfo_pipeline .

```
### 4. Exécution du pipeline

```bash

./start.sh

```
