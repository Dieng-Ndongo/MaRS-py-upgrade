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
├── 📂 pf_3D7_Ref/
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

```

## 📌 Description des principaux composants

| Élément | Description |
|---|---|
| `pipeline_python.py` | Implémentation principale du pipeline bioinformatique en Python. |
| `app.py` | Interface graphique permettant de lancer et suivre les analyses avec Streamlit. |
| `bin/` | Regroupe les scripts secondaires utilisés par le pipeline et l'interface. |
| `pf_3D7_Ref/` | Contient le génome de référence de *Plasmodium falciparum* 3D7. |
| `pf_3D7_snpEff_db/` | Contient les fichiers nécessaires à l'annotation des variants avec SnpEff. |
| `runs/` | Regroupe les résultats générés pour chaque analyse. |
| `environment.yml` | Définit l'environnement Conda et les dépendances du pipeline. |
| `requirements.txt` | Liste les dépendances Python nécessaires à l'interface. |
| `Dockerfile` | Permet de construire l'environnement d'exécution du pipeline avec Docker. |
| `start.sh` | Automatise le lancement du pipeline et/ou de l'interface. |
| `.gitignore` | Définit les fichiers et répertoires qui ne doivent pas être suivis par Git. |

---

## Prérequis
- Système d'exploitation de type Unix (Linux, macOS, etc.) 
- Docker
- Python 3.12 (pour l'interface Streamlit ; l'environnement conda du pipeline reste en Python 3.10 dans le conteneur Docker)

---

## Installation
### 1. Installer Python 3.12
Ubuntu 22.04 n'installe pas Python 3.12 par défaut, il faut passer par le PPA `deadsnakes` :
```bash

sudo apt update && sudo apt install -y software-properties-common

sudo add-apt-repository -y ppa:deadsnakes/ppa

sudo apt update && sudo apt install -y python3.12 python3.12-venv

```
### 2. Installer Docker
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
### 3. Télécharger le dépôt git

```bash

git clone https://github.com/Dieng-Ndongo/MaRS-py-upgrade

cd ~MaRS-py-upgrade

```
### 4. Construction de l'image docker

```bash

docker build -t bioinfo_pipeline .

```
### 5. Installer requirements.txt

```bash

cd ~/MaRS-py-upgrade
pip install -r requirements.txt

```

### 6. Exécution du pipeline

```bash

./start.sh

```

---

## Déploiement en production (systemd + Cloudflare Tunnel)

Pour un déploiement persistant sur le serveur (survit aux déconnexions SSH et aux redémarrages), l'app tourne comme service systemd sous un utilisateur dédié, exposée via un tunnel Cloudflare existant plutôt qu'un reverse proxy public.

### 1. Créer l'utilisateur dédié et déployer le dépôt
```bash

sudo useradd -r -M -d /opt/mars-py-upgrade -s /usr/sbin/nologin mars

sudo git clone https://github.com/Dieng-Ndongo/MaRS-py-upgrade /opt/mars-py-upgrade

sudo usermod -aG docker mars

sudo chown -R mars:mars /opt/mars-py-upgrade

```
Copiez `.streamlit/secrets.toml.example` vers `.streamlit/secrets.toml` dans `/opt/mars-py-upgrade` et remplissez `APP_PASSWORD` et, si besoin, la section `[EMAIL]`.

### 2. Installer le service systemd
```bash

sudo cp /opt/mars-py-upgrade/deploy/mars-streamlit.service /etc/systemd/system/
sudo cp /opt/mars-py-upgrade/deploy/mars-streamlit.timer /etc/systemd/system/

sudo systemctl daemon-reload

sudo systemctl enable --now mars-streamlit.timer

cd /opt/mars-py-upgrade
sudo docker build -t bioinfo_pipeline .

sudo systemctl start mars-streamlit
sudo systemctl status mars-streamlit

journalctl -u mars-streamlit -f

df -h /

```
Le démarrage automatique après un redémarrage du serveur est volontairement retardé de 15 minutes (`mars-streamlit.timer`, `OnBootSec=15min`) pour laisser le reste du serveur (réseau, Docker, autres services) se stabiliser avant de lancer l'app — ce délai ne s'applique qu'au déclenchement après un boot. Pour un démarrage immédiat (première installation, ou reprise manuelle après un arrêt), utilisez `sudo systemctl start mars-streamlit` directement ; la reprise sur échec (`Restart=on-failure`) reste également instantanée (5s), ce délai ne concerne que le boot.

Le premier démarrage crée le venv Python 3.12, installe `requirements.txt` et construit l'image Docker `bioinfo_pipeline` (voir `start.sh`) — cette dernière étape peut prendre plusieurs minutes (résolution Conda d'une vingtaine d'outils bio-informatiques), ce qui est normal. Les démarrages suivants sont quasi instantanés, sauf si `Dockerfile`/`environment.yml` ont changé (reconstruction automatique de l'image).

### 3. Exposer l'app via Cloudflare Tunnel
L'app écoute uniquement en local (`127.0.0.1:8501`, voir `.streamlit/config.toml`) — aucun port public n'est ouvert. Ajoutez une règle d'ingress dans la config `cloudflared` existante du serveur (`/etc/cloudflared/config.yml`) :
```yaml

ingress:
  - hostname: mars.votredomaine.tld
    service: http://localhost:8501
  # ... règles des autres apps déjà en place ...
  - service: http_status:404

```
Puis routez le DNS et relancez le tunnel :
```bash

cloudflared tunnel route dns <nom-du-tunnel> mars.votredomaine.tld

sudo systemctl restart cloudflared

```
