import os
import glob
import subprocess
import shutil
import re
from pathlib import Path
import sys
import warnings
from tqdm import tqdm  # <--- Barre de progression
import time 
import logging 
from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib import colors
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import PageTemplate, Frame
from reportlab.lib.units import mm
from functools import reduce


#======================================================================
# Debut du pipeline
#======================================================================

def QC_pre_trimming(input_dir="data", output_dir="output/QC_pre_trimming", logs_dir="logs"):
    """
    Exécute FastQC sur tous les fichiers .fastq.gz dans input_dir
    et enregistre les rapports dans output_dir.
    Toutes les sorties (stdout/stderr) sont stockées dans un seul fichier log global.
    """
    print("********************** QC_pre_trimming : Début *****************************")

    # Création des répertoires nécessaires
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)

    # Fichier log global
    log_file = os.path.join(logs_dir, "QC_pre_trimming.log")

    # Récupération des fichiers FASTQ
    fastq_files = glob.glob(os.path.join(input_dir, "*.fastq.gz"))

    # Calcul du nombre d'échantillons
    n_samples = len(fastq_files)

    with open(log_file, "w") as lf:
        lf.write("===== LOG GLOBAL - PRE_FASTQC =====\n")
        lf.write(f"Début : {datetime.now()}\n\n")

        if n_samples == 0:
            msg = f"[WARNING] Aucun fichier trouvé dans {input_dir}\n"
            print(msg.strip())
            lf.write(msg)
            return

        print(f"[PIPELINE] Démarrage de l’étape QC_pre_trimming ({n_samples} échantillons détectés)\n")
        lf.write(f"[PIPELINE] Démarrage de l’étape QC_pre_trimming ({n_samples} échantillons détectés)\n")

        start_time = datetime.now()

        # Boucle sur les fichiers avec barre de progression
        for fq in tqdm(fastq_files, desc="Exécution de FastQC", unit="fichier"):
            fq_name = os.path.basename(fq)
            command = ["fastqc", fq, "--outdir", output_dir]

            lf.write(f"\n===== Traitement de {fq_name} =====\n")
            lf.write(f"Commande : {' '.join(command)}\n\n")
            lf.flush()

            try:
                subprocess.run(command, stdout=lf, stderr=lf, text=True, check=True)
                lf.write(f"[INFO] QC_pre_trimming terminé avec succès pour {fq_name}\n")
            except subprocess.CalledProcessError as e:
                error_msg = f"[ERROR] FastQC a échoué pour {fq_name} : {e}\n"
                lf.write(error_msg)
                print(error_msg.strip())

        end_time = datetime.now()
        duration = end_time - start_time

        lf.write(f"\n===== Fin du processus =====\n")
        lf.write(f"Durée totale : {duration.seconds} secondes\n")
        lf.write(f"Fin : {datetime.now()}\n")

    print(f"\n[PIPELINE] Étape QC_pre_trimming terminée en {duration.seconds} secondes.")
    print(f"[INFO] Log global disponible ici : {log_file}")
    print("********************** QC_pre_trimming : Fin *****************************")



def MultiQC_pre_trimming(output_dir="output/QC_pre_trimming", report_name="MultiQC_pre_trimming_report.html", logs_dir="logs"):
    """
    Lance MultiQC sur le dossier de sortie (output_dir)
    et génère un rapport global.
    Toutes les sorties (stdout/stderr) sont enregistrées dans un fichier log global.
    """
    print("********************** MultiQC_pre_trimming : Début *****************************")

    # Création du dossier logs s’il n’existe pas
    os.makedirs(logs_dir, exist_ok=True)

    # Fichier log global
    log_file = os.path.join(logs_dir, "MultiQC_pre_trimming.log")

    # Commande MultiQC
    command = ["multiqc", output_dir, "-o", output_dir, "-n", report_name, "--force"]

    start_time = datetime.now()
    with open(log_file, "w") as lf:
        lf.write("===== LOG GLOBAL - PRE_MULTIQC =====\n")
        lf.write(f"Début : {start_time}\n\n")
        lf.write(f"[INFO] Commande exécutée : {' '.join(command)}\n\n")

        print(f"[PIPELINE] Démarrage de l’étape MultiQC_pre_trimming sur le dossier : {output_dir}\n")

        # Barre de progression simulée (MultiQC est une seule commande)
        for _ in tqdm(range(1), desc="Exécution de MultiQC", unit="rapport"):
            try:
                subprocess.run(command, stdout=lf, stderr=lf, text=True, check=True)
                lf.write(f"\n[INFO] MultiQC_pre_trimming terminé avec succès.\n")
            except subprocess.CalledProcessError as e:
                lf.write(f"\n[ERROR] Échec de MultiQC : {e}\n")
                print("[ERROR] MultiQC_pre_trimming a échoué (voir log).")
                break

        end_time = datetime.now()
        duration = end_time - start_time

        lf.write(f"\n===== Fin du processus =====\n")
        lf.write(f"Durée totale : {duration.seconds} secondes\n")
        lf.write(f"Fin : {end_time}\n")

    print(f"\n[PIPELINE] Étape MultiQC_pre_trimming terminée en {duration.seconds} secondes.")
    print(f"[INFO] Rapport généré : {os.path.join(output_dir, report_name)}")
    print(f"[INFO] Log global disponible ici : {log_file}")
    print("********************** MultiQC_pre_trimming : Fin *****************************")




def trimming(input_dir="data", output_dir="output/trimmed_reads", adapter_file="pf_3D7_Ref/adapters.fa", logs_dir="logs"):
    """
    Exécute BBduk pour le trimming des reads paired-end dans input_dir.
    Enregistre les sorties dans output_dir et un log global dans logs_dir.
    Affiche la progression et le temps total d'exécution.
    """
    print("********************** TRIMMING : Début *****************************")

    # Vérification BBduk
    if shutil.which("bbduk.sh") is None:
        print("[ERROR] BBduk n'est pas installé ou pas dans le PATH.")
        return []

    # Création des dossiers
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)

    # Fichier log 
    log_file = os.path.join(logs_dir, f"trimming.log")

    reads_R1 = sorted(glob.glob(os.path.join(input_dir, "*_R1*.fastq*")))
    reads_R2 = sorted(glob.glob(os.path.join(input_dir, "*_R2*.fastq*")))

    if not reads_R1 or not reads_R2:
        print(f"[WARNING] Aucun fichier R1 ou R2 trouvé dans {input_dir}")
        return []

    if len(reads_R1) != len(reads_R2):
        print("[WARNING] Le nombre de fichiers R1 et R2 ne correspond pas !")
        return []

    total_samples = len(reads_R1)
    print(f"[PIPELINE] Démarrage de l’étape Trimming ({total_samples} échantillons détectés)\n")

    trimmed_files = []
    start_time = time.time()

    with open(log_file, "w") as lf:
        lf.write("===== LOG GLOBAL - TRIMMING =====\n")
        lf.write(f"Début : {datetime.now()}\n")
        lf.write(f"Nombre d'échantillons : {total_samples}\n\n")

        for r1, r2 in tqdm(zip(reads_R1, reads_R2), total=total_samples, desc="Trimming en cours", unit="échantillon"):
            sample_id = os.path.basename(r1)
            sample_id = re.sub(r"_R1.*\.fastq.*", "", sample_id)

            out_r1 = os.path.join(output_dir, f"{sample_id}_trimmed_R1.fastq.gz")
            out_r2 = os.path.join(output_dir, f"{sample_id}_trimmed_R2.fastq.gz")
            stats_file = os.path.join(output_dir, f"{sample_id}.stats.txt")

            cmd = [
                "bbduk.sh", "-Xmx1g", "ktrim=r", "k=27", "hdist=1", "edist=0",
                f"ref={adapter_file}", "qtrim=rl", "trimq=30", "minlength=50",
                "trimbyoverlap=t", "minoverlap=24", "qin=33",
                f"in={r1}", f"in2={r2}",
                f"out={out_r1}", f"out2={out_r2}", f"stats={stats_file}"
            ]

            lf.write(f"\n[INFO] Trimming de l'échantillon : {sample_id}\n")
            lf.write(f"Commande : {' '.join(cmd)}\n")

            try:
                subprocess.run(cmd, stdout=lf, stderr=lf, text=True, check=True)
                lf.write(f"[OK] Trimming terminé pour {sample_id}\n")
                trimmed_files.extend([out_r1, out_r2])
            except subprocess.CalledProcessError as e:
                lf.write(f"[ERROR] Échec du trimming pour {sample_id} : {e}\n")

        total_time = time.time() - start_time
        lf.write(f"\n===== Fin du processus =====\n")
        lf.write(f"Durée totale : {total_time/60:.2f} minutes\n")
        lf.write(f"Fin : {datetime.now()}\n")

    print(f"\n[PIPELINE] Étape Trimming terminée en {total_time/60:.2f} minutes.")
    print(f"[INFO] Log global disponible ici : {log_file}")
    print("********************** TRIMMING : Fin *****************************")

    return trimmed_files


def QC_post_trimming(input_dir="output/trimmed_reads", output_dir="output/QC_post_trimming", logs_dir="logs"):
    """
    Exécute FastQC sur tous les fichiers .fastq.gz dans input_dir
    et enregistre les rapports dans output_dir + un log global dans logs_dir.
    Affiche la progression et le temps total d'exécution.
    """
    print("********************** QC_post_trimming : Début *****************************")

    # Création des dossiers nécessaires
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)

    # Récupération des fichiers FASTQ
    fastq_files = glob.glob(os.path.join(input_dir, "*.fastq.gz"))

    if not fastq_files:
        print(f"[WARNING] Aucun fichier trouvé dans {input_dir}")
        return

    total_samples = len(fastq_files)
    print(f"[PIPELINE] Démarrage de l’étape QC_post_trimming ({total_samples} fichiers détectés)\n")

    # Préparation du fichier log
    log_file = os.path.join(logs_dir, f"QC_post_trimming.log")

    start_time = time.time()

    with open(log_file, "w") as lf:
        lf.write("===== LOG GLOBAL - QC_post_trimming =====\n")
        lf.write(f"Début : {datetime.now()}\n")
        lf.write(f"Nombre de fichiers : {total_samples}\n\n")

        # Boucle avec barre de progression
        for fq in tqdm(fastq_files, desc="Exécution FastQC", unit="fichier"):
            file_name = os.path.basename(fq)
            lf.write(f"\n[INFO] Analyse de : {file_name}\n")

            command = ["fastqc", fq, "--outdir", output_dir]
            lf.write(f"Commande : {' '.join(command)}\n")

            try:
                subprocess.run(command, stdout=lf, stderr=lf, text=True, check=True)
                lf.write(f"[OK] QC_post_trimming terminé pour {file_name}\n")
            except subprocess.CalledProcessError as e:
                lf.write(f"[ERROR] QC_post_trimming a échoué pour {file_name} : {e}\n")

        total_time = time.time() - start_time
        lf.write(f"\n===== Fin du processus =====\n")
        lf.write(f"Durée totale : {total_time/60:.2f} minutes\n")
        lf.write(f"Fin : {datetime.now()}\n")

    print(f"\n[PIPELINE] Étape QC_post_trimming terminée en {total_time/60:.2f} minutes.")
    print(f"[INFO] Log global disponible ici : {log_file}")
    print("********************** QC_post_trimming : Fin *****************************")



def MultiQC_post_trimming(output_dir="output/QC_post_trimming", report_name="MultiQC_post_trimming_report.html", logs_dir="logs"):
    """
    Lance MultiQC sur le dossier de sortie (output_dir),
    génère un rapport global et crée un log complet dans logs_dir.
    """
    print("********************** MultiQC_post_trimming : Début *****************************")

    # Création des répertoires nécessaires
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)

    log_file = os.path.join(logs_dir, "MultiQC_post_trimming.log")

    start_time = time.time()

    with open(log_file, "w") as lf:
        lf.write("===== LOG GLOBAL - MultiQC_post_trimming =====\n")
        lf.write(f"Début : {datetime.now()}\n")
        lf.write(f"Dossier d’analyse : {output_dir}\n\n")

        command = ["multiqc", output_dir, "-o", output_dir, "-n", report_name, "--force"]
        lf.write(f"[INFO] Commande exécutée : {' '.join(command)}\n\n")

        print(f"[PIPELINE] Démarrage de l’étape MultiQC_post_trimming sur le dossier {output_dir}\n")

        try:
            # Barre de progression (même si une seule tâche)
            for _ in tqdm(range(1), desc="Exécution de MultiQC", unit="rapport"):
                subprocess.run(command, stdout=lf, stderr=lf, text=True, check=True)

            lf.write(f"[OK] Rapport MultiQC_post_trimming généré : {os.path.join(output_dir, report_name)}\n")
            print(f"[INFO] Rapport MultiQC_post_trimming généré : {os.path.join(output_dir, report_name)}")

        except subprocess.CalledProcessError as e:
            lf.write(f"[ERROR] MultiQC a échoué : {e}\n")
            print(f"[ERROR] MultiQC a échoué : {e}")

        total_time = time.time() - start_time
        lf.write(f"\n===== Fin du processus =====\n")
        lf.write(f"Durée totale : {total_time/60:.2f} minutes\n")
        lf.write(f"Fin : {datetime.now()}\n")

    print(f"\n[PIPELINE] Étape MultiQC_post_trimming terminée en {total_time/60:.2f} minutes.")
    print(f"[INFO] Log global disponible ici : {log_file}")
    print("********************** MultiQC_post_trimming : Fin *****************************")



def run_cmd(cmd, log_file):
    """Exécute une commande en redirigeant stdout/stderr vers log_file"""
    with open(log_file, "a") as lf:
        lf.write(f"Commande : {' '.join(cmd)}\n")
        try:
            subprocess.run(cmd, stdout=lf, stderr=lf, text=True, check=True)
            lf.write("[OK] Commande terminée\n")
        except subprocess.CalledProcessError as e:
            lf.write(f"[ERROR] Échec de la commande : {e}\n")
            raise



def bwa_index(ref, output_dir="output/index_ref", logs_dir="logs"):
    """
    Indexation de la séquence de référence :
      - BWA index
      - Samtools faidx
      - Picard CreateSequenceDictionary (via Conda)
    """
    print("********************** BWA_INDEX : Début *****************************")
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)

    log_file = os.path.join(logs_dir, "bwa_index.log")
    start_time = time.time()

    with open(log_file, "w") as lf:
        lf.write("===== LOG GLOBAL - BWA_INDEX =====\n")
        lf.write(f"Début : {datetime.now()}\n")
        lf.write(f"Fichier de référence : {ref}\n\n")

    # Copier la référence
    ref_copy = os.path.join(output_dir, os.path.basename(ref))
    if not os.path.exists(ref_copy):
        shutil.copy(ref, ref_copy)
    print(f"[INFO] Référence copiée dans {output_dir}")

    steps = [
        ("BWA index", ["bwa", "index", ref_copy]),
        ("Samtools faidx", ["samtools", "faidx", ref_copy]),
        ("Picard CreateSequenceDictionary", ["picard", "CreateSequenceDictionary"])
    ]

    for step_name, cmd in tqdm(steps, desc="Indexation en cours", unit="étape"):
        print(f"[INFO] {step_name} en cours...")

        if step_name != "Picard CreateSequenceDictionary":
            run_cmd(cmd, log_file)
        else:
            dict_output = ref_copy.replace(".fasta", ".dict").replace(".fa", ".dict")
            if os.path.exists(dict_output):
                os.remove(dict_output)

            picard_cmd = shutil.which("picard")
            if not picard_cmd:
                with open(log_file, "a") as lf:
                    lf.write("[ERROR] Picard non trouvé dans l'environnement Conda\n")
                raise FileNotFoundError("Picard non trouvé dans l'environnement Conda")

            run_cmd([
                picard_cmd,
                "CreateSequenceDictionary",
                f"R={ref_copy}",
                f"O={dict_output}"
            ], log_file)

    total_time = time.time() - start_time
    print(f"[PIPELINE] Étape bwa_index terminée en {total_time/60:.2f} minutes")
    print(f"[INFO] Log global disponible ici : {log_file}")
    print("********************** BWA_INDEX : Fin *****************************")




def bwa_index_0(ref, output_dir="output/index_ref", local_picard_jar="tools/picard.jar", logs_dir="logs"):
    """
    Indexation de la séquence de référence :
      - BWA index
      - Samtools faidx
      - Picard CreateSequenceDictionary
    Affiche uniquement la progression et le temps final dans le terminal.
    Tout le reste (stdout/stderr, erreurs) va dans logs/bwa_index.log
    """
    print("********************** BWA_INDEX : Début *****************************")
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)

    log_file = os.path.join(logs_dir, "bwa_index.log")
    start_time = time.time()

    with open(log_file, "w") as lf:
        lf.write("===== LOG GLOBAL - BWA_INDEX =====\n")
        lf.write(f"Début : {datetime.now()}\n")
        lf.write(f"Fichier de référence : {ref}\n\n")

    # Copier la référence
    ref_copy = os.path.join(output_dir, os.path.basename(ref))
    if not os.path.exists(ref_copy):
        shutil.copy(ref, ref_copy)
    print(f"[INFO] Référence copiée dans {output_dir}")

    # Étapes
    steps = [
        ("BWA index", ["bwa", "index", ref_copy]),
        ("Samtools faidx", ["samtools", "faidx", ref_copy]),
        ("Picard CreateSequenceDictionary", None)
    ]

    for step_name, cmd in tqdm(steps, desc="Indexation en cours", unit="étape"):
        print(f"[INFO] {step_name} en cours...")
        if cmd:
            run_cmd(cmd, log_file)
        else:
            # Picard
            dict_output = ref_copy.replace(".fasta", ".dict").replace(".fa", ".dict")
            if os.path.exists(dict_output):
                os.remove(dict_output)

            picard_cmd = shutil.which("picard")
            if picard_cmd:
                run_cmd([picard_cmd, "CreateSequenceDictionary", f"R={ref_copy}", f"O={dict_output}"], log_file)
            elif os.path.exists(local_picard_jar):
                run_cmd(["java", "-jar", local_picard_jar, "CreateSequenceDictionary", f"R={ref_copy}", f"O={dict_output}"], log_file)
            else:
                with open(log_file, "a") as lf:
                    lf.write("[ERROR] Picard non trouvé !\n")
                raise FileNotFoundError("Picard non trouvé")

    total_time = time.time() - start_time
    print(f"[PIPELINE] Étape bwa_index terminée en {total_time/60:.2f} minutes")
    print(f"[INFO] Log global disponible ici : {log_file}")
    print("********************** BWA_INDEX : Fin *****************************")


def bwa_align(input_dir="output/trimmed_reads", sam_dir="output/sam_files", bam_dir="output/bam_files", 
                reference="pf_3D7_Ref/mars_pf_ref.fasta", logs_dir="logs"):
    """
    Aligne des reads paired-end sur un génome de référence avec BWA-MEM,
    produit SAM et BAM trié et indexé.
    Log global dans logs_dir, affichage minimal dans le terminal.
    """
    print("********************** BWA_ALIGN : Début *****************************")
    Path(sam_dir).mkdir(parents=True, exist_ok=True)
    Path(bam_dir).mkdir(parents=True, exist_ok=True)
    Path(logs_dir).mkdir(parents=True, exist_ok=True)

    log_file = os.path.join(logs_dir, "bwa_align.log")
    fastq_files = sorted(Path(input_dir).glob("*_trimmed_R1.fastq.gz"))
    total_samples = len(fastq_files)

    print(f"[PIPELINE] Démarrage de l’étape bwa_align ({total_samples} échantillons détectés)\n")
    start_time = time.time()

    with open(log_file, "w") as lf:
        lf.write("===== LOG GLOBAL - BWA_ALIGN =====\n")
        lf.write(f"Début : {datetime.now()}\n")
        lf.write(f"Nombre d'échantillons : {total_samples}\n")
        lf.write(f"Référence : {reference}\n\n")

    for r1_file in tqdm(fastq_files, desc="Alignement BWA", unit="échantillon"):
        sample_name = r1_file.name.replace("_trimmed_R1.fastq.gz", "")
        r2_file = r1_file.with_name(r1_file.name.replace("_R1", "_R2"))

        if not r2_file.exists():
            print(f"[WARNING] Paire pour {r1_file.name} non trouvée. Ignoré.")
            with open(log_file, "a") as lf:
                lf.write(f"[WARNING] Paire pour {r1_file.name} non trouvée. Ignoré.\n")
            continue

        sam_file = Path(sam_dir) / f"{sample_name}.sam"
        bam_file = Path(bam_dir) / f"{sample_name}.sorted.bam"
        bai_file = Path(bam_dir) / f"{sample_name}.sorted.bam.bai"

        # Étape 1 : SAM
        with open(log_file, "a") as lf:
            cmd_sam = ["bwa", "mem", "-t", "4", reference, str(r1_file), str(r2_file)]
            lf.write(f"\n[INFO] Génération SAM pour {sample_name}\n")
            lf.write(f"Commande : {' '.join(cmd_sam)}\n")
            try:
                subprocess.run(cmd_sam, stdout=open(sam_file, "w"), stderr=lf, text=True, check=True)
                lf.write(f"[OK] SAM généré : {sam_file}\n")
            except subprocess.CalledProcessError as e:
                lf.write(f"[ERROR] SAM échoué pour {sample_name} : {e}\n")
        print(f"[INFO] SAM généré : {sam_file}")

        # Étape 2 : BAM trié
        cmd_bam = f"samtools view -@ 4 -bS {sam_file} | samtools sort -@ 4 -o {bam_file}"
        with open(log_file, "a") as lf:
            lf.write(f"\n[INFO] Tri BAM pour {sample_name}\n")
            lf.write(f"Commande : {cmd_bam}\n")
            try:
                subprocess.run(cmd_bam, shell=True, stdout=lf, stderr=lf, text=True, check=True, executable="/bin/bash")
                lf.write(f"[OK] BAM trié généré : {bam_file}\n")
            except subprocess.CalledProcessError as e:
                lf.write(f"[ERROR] BAM trié échoué pour {sample_name} : {e}\n")
        print(f"[INFO] BAM trié généré : {bam_file}")

        # Étape 3 : BAM index
        with open(log_file, "a") as lf:
            cmd_index = ["samtools", "index", str(bam_file)]
            lf.write(f"\n[INFO] Index BAM pour {sample_name}\n")
            lf.write(f"Commande : {' '.join(cmd_index)}\n")
            try:
                subprocess.run(cmd_index, stdout=lf, stderr=lf, text=True, check=True)
                lf.write(f"[OK] BAM index généré : {bai_file}\n")
            except subprocess.CalledProcessError as e:
                lf.write(f"[ERROR] BAM index échoué pour {sample_name} : {e}\n")
        print(f"[INFO] BAM index généré : {bai_file}")

    total_time = time.time() - start_time
    print(f"\n[PIPELINE] Étape bwa_align terminée en {total_time/60:.2f} minutes")
    print(f"[INFO] Log global disponible ici : {log_file}")
    print("********************** BWA_ALIGN : Fin *****************************")


def picard_add_readgroups(input_dir="output/sam_files", output_dir="output/bam_picard_readgroups", 
                             local_picard_jar="tools/picard.jar", logs_dir="logs"):
    """
    Ajoute les read groups à tous les fichiers SAM d'un dossier avec Picard AddOrReplaceReadGroups.
    Log global dans logs_dir, affichage minimal dans le terminal.
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    Path(logs_dir).mkdir(parents=True, exist_ok=True)
    log_file = os.path.join(logs_dir, "picard_add_readgroups.log")

    sam_files = sorted(Path(input_dir).glob("*.sam"))
    total_samples = len(sam_files)

    if total_samples == 0:
        print(f"[WARNING] Aucun fichier SAM trouvé dans {input_dir}")
        return

    print(f"[PIPELINE] Démarrage de l’étape picard_add_readgroups ({total_samples} fichiers SAM détectés)\n")
    start_time = time.time()

    with open(log_file, "w") as lf:
        lf.write("===== LOG GLOBAL - PICARD_ADD_READGROUPS =====\n")
        lf.write(f"Début : {datetime.now()}\n")
        lf.write(f"Nombre de fichiers SAM : {total_samples}\n\n")

    for sam_file in tqdm(sam_files, desc="Ajout Read Groups", unit="fichier"):
        sample_id = sam_file.stem
        bam_out = Path(output_dir) / f"{sample_id}_picard_readgroup.bam"

        # Détection de Picard (conda ou jar local)
        picard_cmd = shutil.which("picard")
        if picard_cmd:
            cmd = [picard_cmd, "AddOrReplaceReadGroups",
                   f"I={sam_file}", f"O={bam_out}",
                   "SORT_ORDER=coordinate",
                   "RGLB=ExomeSeq", "RGPL=Illumina", "RGPU=NextSeq",
                   f"RGSM={sample_id}", "CREATE_INDEX=true"]
        elif os.path.exists(local_picard_jar):
            cmd = ["java", "-jar", local_picard_jar, "AddOrReplaceReadGroups",
                   f"I={sam_file}", f"O={bam_out}",
                   "SORT_ORDER=coordinate",
                   "RGLB=ExomeSeq", "RGPL=Illumina", "RGPU=NextSeq",
                   f"RGSM={sample_id}", "CREATE_INDEX=true"]
        else:
            raise FileNotFoundError(
                "Picard non trouvé ! Installez-le via conda ou placez picard.jar dans tools/picard.jar"
            )

        # Exécution avec log global
        with open(log_file, "a") as lf:
            lf.write(f"\n[INFO] Traitement {sam_file.name}\n")
            lf.write(f"Commande : {' '.join(cmd)}\n")
            try:
                subprocess.run(cmd, stdout=lf, stderr=lf, text=True, check=True)
                lf.write(f"[OK] BAM avec read groups généré : {bam_out}\n")
            except subprocess.CalledProcessError as e:
                lf.write(f"[ERROR] Picard a échoué pour {sam_file.name} : {e}\n")

        print(f"[INFO] BAM avec read groups généré : {bam_out}")

    total_time = time.time() - start_time
    print(f"\n[PIPELINE] Étape picard_add_readgroups terminée en {total_time/60:.2f} minutes")
    print(f"[INFO] Log global disponible ici : {log_file}")
    print("********************** PICARD_ADD_READGROUPS : Fin *****************************")



def get_bed(input_dir="pf_3D7_Ref", output_dir="output/BED", bed_name="mars_genes.bed", logs_dir="logs"):
    """
    Convertit tous les fichiers GFF d'un dossier en fichiers BED contenant uniquement les gènes.
    Les sorties (stdout, stderr) sont enregistrées dans un fichier .log.
    Le terminal affiche uniquement les informations essentielles, la progression et le temps d’exécution.
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    start_time = time.time()
    Path(logs_dir).mkdir(parents=True, exist_ok=True)
    log_file = os.path.join(logs_dir, "get_bed.log")

    # Lister tous les fichiers GFF dans le dossier
    gff_files = sorted(Path(input_dir).glob("*.gff"))
    nb_samples = len(gff_files)

    if nb_samples == 0:
        print(f"[WARNING] Aucun fichier GFF trouvé dans {input_dir}")
        return

    print(f"\n[PIPELINE] Démarrage de l’étape get_bed : {nb_samples} échantillons détectés\n")

    # Exécution pour chaque fichier GFF
    for gff_file in tqdm(gff_files, desc="[PROGRESSION] Conversion GFF → BED", ncols=100):
        bed_file = Path(output_dir) / f"{gff_file.stem}_{bed_name}"

        cmd = (
            f"awk 'BEGIN {{ OFS=\"\\t\" }} "
            f"{{if ($3==\"gene\") {{print $1,$4-1,$5,$10,$16,$7}}}}' {gff_file} > {bed_file}"
        )

    # Exécution silencieuse avec redirection des logs
    with open(log_file, "a") as log:
        log.write(f"\n[COMMAND] {cmd}\n")
        subprocess.run(cmd, shell=True, stdout=log, stderr=log, check=True)

        
    end_time = time.time()
    duration = end_time - start_time

    print(f"\n[INFO] Conversion terminée pour {nb_samples} échantillons.")
    print(f"[INFO] Fichiers BED enregistrés dans : {output_dir}")
    print(f"[TEMPS] Durée totale d’exécution : {duration:.2f} secondes\n")

    return output_dir



def vcf_call_all_samples(bam_dir="output/bam_picard_readgroups", ref="output/index_ref/mars_pf_ref.fasta", bed_file="output/BED/mars_pf_mars_genes.bed",
                            output_dir="output/vcf_files", logs_dir="logs"):
    """
    Appelle les variants pour tous les échantillons dans bam_dir.
    Génère pour chaque échantillon les fichiers VCF : samtools, Freebayes, GATK, VarDict.
    Toutes les sorties sont enregistrées dans un fichier log.
    """
    print("********************** VCF_CALL_ALL_SAMPLES : Début *****************************")

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    Path(logs_dir).mkdir(parents=True, exist_ok=True)

    bam_files = sorted(Path(bam_dir).glob("*_picard_readgroup.bam"))

    if not bam_files:
        print(f"[WARNING] Aucun BAM trouvé dans {bam_dir}")
        return

    total_samples = len(bam_files)
    print(f"[PIPELINE] Démarrage de l’étape vcf_call_all_samples : {total_samples} échantillons détectés\n")

    log_file = Path(logs_dir) / "vcf_call_all_samples.log"
    start_time = time.time()

    with open(log_file, "w") as lf:
        lf.write("===== LOG GLOBAL - VCF_CALL_ALL_SAMPLES =====\n")
        lf.write(f"Début : {datetime.now()}\n")
        lf.write(f"Nombre d’échantillons : {total_samples}\n\n")

        for bam in tqdm(bam_files, desc="Appel des variants", unit="échantillon"):
            sample_id = bam.stem.replace("_picard_readgroup", "")
            lf.write(f"\n===== Échantillon : {sample_id} =====\n")

            samtools_vcf = Path(output_dir) / f"{sample_id}_samtools.vcf"
            freebayes_vcf = Path(output_dir) / f"{sample_id}_Freebayes.vcf"
            gatk_vcf = Path(output_dir) / f"{sample_id}_gatk.vcf"
            vardict_vcf = Path(output_dir) / f"{sample_id}_vardict.vcf"
            mpileup_file = Path(output_dir) / f"{sample_id}.mpileup"

            try:
                # Indexation BAM
                subprocess.run(["samtools", "index", str(bam)],
                               stdout=lf, stderr=lf, text=True, check=True)

                # Samtools variant calling
                subprocess.run(["bcftools", "mpileup", "-O", "b", "-d", "10000", "-o", str(mpileup_file),
                                "-f", str(ref), str(bam)],
                                stdout=lf, stderr=lf, text=True, check=True)
                subprocess.run(["bcftools", "call", "-mv", "-P", "0.05",
                                "-o", str(samtools_vcf), str(mpileup_file)],
                               stdout=lf, stderr=lf, text=True, check=True)

                # Freebayes
                subprocess.run(["freebayes", "-f", str(ref), "-F", "0.05", "-E", "3",
                                "--report-all-haplotype-alleles", str(bam)],
                               stdout=open(freebayes_vcf, "w"), stderr=lf, text=True, check=True)

                # GATK
                subprocess.run(["gatk", "HaplotypeCaller", "--native-pair-hmm-threads", "8",
                                "-R", str(ref), "-I", str(bam),
                                "--min-base-quality-score", "0", "-O", str(gatk_vcf)],
                                stdout=lf, stderr=lf, text=True, check=True)

                # VarDict
                vardict_cmd = f"vardict -G {ref} -f 0.05 -N {sample_id} -b {bam} -c 1 -S 2 -E 3 -g 4 {bed_file} | var2vcf_valid.pl -N {sample_id} -E -f 0.05 > {vardict_vcf}"
                subprocess.run(vardict_cmd, shell=True, stdout=lf, stderr=lf, text=True, check=True)

            
            
                lf.write(f"[OK] Variants générés pour {sample_id}\n")
            except subprocess.CalledProcessError as e:
                lf.write(f"[ERROR] Appel de variants échoué pour {sample_id} : {e}\n")
                continue

            lf.write(f"[INFO] Fichiers VCF générés pour {sample_id} :\n"
                     f"  SAMTOOLS : {samtools_vcf}\n"
                     f"  FREEBAYES : {freebayes_vcf}\n"
                     f"  GATK : {gatk_vcf}\n"
                     f"  VARDICT : {vardict_vcf}\n")

        total_time = time.time() - start_time
        lf.write(f"\n===== Fin du processus =====\n")
        lf.write(f"Durée totale : {total_time/60:.2f} minutes\n")
        lf.write(f"Fin : {datetime.now()}\n")

    print(f"\n[PIPELINE] Étape vcf_call_all_samples terminée en {total_time/60:.2f} minutes.")
    print(f"[INFO] Log global disponible ici : {log_file}")
    print("********************** VCF_CALL_ALL_SAMPLES : Fin *****************************")



def build_snpeff_db(db_name, snpeff_config, output_dir=None):
    """
    Construit une base de données snpEff en utilisant snpEff installé via conda.
    
    Arguments :
        db_name : nom de la base de données
        snpeff_config : répertoire contenant snpEff.config et genes.gbk
        output_dir : dossier de sortie (par défaut : répertoire courant)
    """
    start_time = time.time()
    step_name = "build_snpeff_db"

    # Création du dossier de logs
    log_dir = Path("logs")
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{step_name}.log"

    with open(log_file, "w") as lf:
        print(f"[PIPELINE] Démarrage de l’étape {step_name}")
        lf.write(f"===== Lancement de l’étape {step_name} =====\n")

        if output_dir is None:
            output_dir = os.getcwd()
        db_dir = Path(output_dir) / db_name 
        inner_db_dir = db_dir / db_name

        # Création des répertoires
        inner_db_dir.mkdir(parents=True, exist_ok=True)

        # Copier les fichiers nécessaires
        config_src = Path(snpeff_config) / "snpEff.config"
        gbk_src = Path(snpeff_config) / "genes.gbk"

        config_dst = db_dir / "snpEff.config"
        gbk_dst = inner_db_dir / "genes.gbk"

        shutil.copy(config_src, config_dst)
        shutil.copy(gbk_src, gbk_dst)

        # Commande snpEff
        cmd = [
            "snpEff",
            "build",
            "-c", "snpEff.config",
            "-genbank", "-v", db_name
        ]

        lf.write(f"[INFO] Commande exécutée : {' '.join(cmd)}\n")

        # Barre de progression
        with tqdm(total=1, desc="Construction de la base snpEff", ncols=100) as pbar:
            try:
                subprocess.run(cmd, cwd=str(db_dir), stdout=lf, stderr=lf, text=True, check=True)
                pbar.update(1)
                print(f"[SUCCESS] Base de données snpEff '{db_name}' construite dans {db_dir}")
                lf.write(f"[SUCCESS] Base de données snpEff '{db_name}' construite dans {db_dir}\n")
            except subprocess.CalledProcessError as e:
                pbar.close()
                print(f"[ERROR] Erreur pendant la construction de la DB snpEff : {e}")
                lf.write(f"[ERROR] Erreur pendant la construction : {e}\n")

        elapsed = time.time() - start_time 
        print(f"[INFO] Temps d’exécution : {elapsed:.2f} secondes")
        lf.write(f"Temps d’exécution : {elapsed:.2f} secondes\n")

    return db_dir



def annotate_all_vcfs(vcf_dir = "output/vcf_files", db_name = "pf_3D7_snpEff_db", db_dir = Path("output") / "pf_3D7_snpEff_db", 
                        output_dir= "output/vcf_annotated", logs_dir="logs"):
    
    """
    Annoter tous les fichiers VCF dans vcf_dir en utilisant une base snpEff déjà créée.
    Toutes les sorties (stdout + stderr) sont redirigées vers un fichier log.
    Seules les informations essentielles, la progression et le temps d'exécution
    s'affichent dans le terminal.
    """

    start_time = time.time()

    vcf_dir = Path(vcf_dir).resolve()
    db_dir = Path(db_dir).resolve()

    # Dossier de sortie
    if output_dir is None:
        output_dir = vcf_dir
    else:
        output_dir = Path(output_dir).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)

    # Dossier de logs
    logs_dir = Path(logs_dir)
    logs_dir.mkdir(parents=True, exist_ok=True)
    log_file = logs_dir / "annotate_all_vcfs.log"

    annotated_files = {}
    config_file = db_dir / "snpEff.config"

    vcf_files = sorted(vcf_dir.glob("*.vcf"))
    n_samples = len(vcf_files)

    print(f"[PIPELINE] Démarrage de l’étape annotate_all_vcfs ({n_samples} échantillons détectés)")

    with open(log_file, "a") as lf:
        for vcf_file in tqdm(vcf_files, desc="Annotation des variants", unit="échantillon"):
            sample_id = vcf_file.stem
            annotated_vcf = output_dir / f"{sample_id}_annot.vcf"

            cmd = [
                "snpEff",
                "-c", str(config_file),
                "-hgvs1LetterAa",
                "-noShiftHgvs",
                db_name,
                str(vcf_file.resolve())
            ]

            lf.write(f"\n[INFO] Annotation de {vcf_file.name}...\n")
            lf.write(f"Commande : {' '.join(cmd)}\n")

            try:
                with open(annotated_vcf, "w") as out:
                    subprocess.run(cmd, stdout=out, stderr=lf, text=True, check=True)

                annotated_files[sample_id] = str(annotated_vcf)
                print(f"[INFO] Fichier annoté généré : {annotated_vcf}")
                lf.write(f"[SUCCESS] Annotation terminée pour {vcf_file.name}\n")

            except subprocess.CalledProcessError as e:
                print(f"[ERROR] Échec de l’annotation pour {vcf_file.name}")
                lf.write(f"[ERROR] Échec de l’annotation pour {vcf_file.name} : {e}\n")

    elapsed_time = time.time() - start_time
    print(f"[PIPELINE] Étape annotate_all_vcfs terminée en {elapsed_time:.2f} secondes.")
    print(f"[LOG] Détails complets enregistrés dans : {log_file}")

    return annotated_files




import json
import subprocess
from pathlib import Path

from pathlib import Path

def find_snpsift_jar(env_name="pipeline_env"):
    """
    Retourne le chemin absolu du SnpSift.jar dans l'environnement Micromamba
    """
    base_path = Path("/opt/conda/envs") / env_name / "share"
    # Chercher un dossier contenant 'snpsift'
    snpsift_dirs = [d for d in base_path.glob("snpsift*") if d.is_dir()]
    if not snpsift_dirs:
        raise FileNotFoundError(f"SnpSift.jar introuvable dans {base_path}")
    snpsift_jar = snpsift_dirs[0] / "SnpSift.jar"
    if not snpsift_jar.exists():
        raise FileNotFoundError(f"SnpSift.jar introuvable à : {snpsift_jar}")
    return str(snpsift_jar)





def run_vartype_all(input_dir, output_dir, env_path=None, logs_dir=None):
    """
    Applique SnpSift varType sur tous les fichiers VCF annotés présents dans input_dir.
    Affiche la progression et enregistre toutes les sorties dans un log global.
    """
    start_time = time.time()
    os.makedirs(output_dir, exist_ok=True)

    # Dossier de logs
    if logs_dir is None:
        logs_dir = Path("logs")
    logs_dir = Path(logs_dir)
    logs_dir.mkdir(parents=True, exist_ok=True)
    log_file = os.path.join(logs_dir, "vartype.log")

    # Ouvrir le log global
    with open(log_file, "w") as global_log:
        global_log.write("[PIPELINE] --- DÉBUT DE L'ÉTAPE VARTYPE ---\n")
        global_log.write(f"Début : {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        # Localisation du jar
        """
        if env_path:
            snpsift_jar = Path(env_path) / "share" / "snpsift-5.3.0a-0" / "SnpSift.jar"
        else:
            raise FileNotFoundError("Veuillez spécifier env_path vers votre environnement conda contenant SnpSift.jar")
        if not snpsift_jar.exists():
            raise FileNotFoundError(f"SnpSift.jar introuvable à : {snpsift_jar}")
        """
        # Localisation du jar
        # Localisation du jar SnpSift
        snpsift_jar = find_snpsift_jar(env_name="pipeline_env")
        print(f"[INFO] SnpSift.jar trouvé à : {snpsift_jar}")



        

        # Liste des fichiers annotés
        annotated_vcfs = glob.glob(os.path.join(input_dir, "*_annot.vcf"))

        # Extraction des IDs d’échantillons
        sample_ids = set()
        for vcf in annotated_vcfs:
            fname = os.path.basename(vcf)
            for suffix in ["_samtools_annot.vcf", "_Freebayes_annot.vcf", "_gatk_annot.vcf", "_vardict_annot.vcf"]:
                if fname.endswith(suffix):
                    sample_ids.add(fname.replace(suffix, ""))

        print(f"[PIPELINE] Démarrage de l’étape VarType ({len(sample_ids)} échantillons détectés)\n")
        global_log.write(f"{len(sample_ids)} échantillons détectés.\n\n")

        results = {}

        # Barre de progression sur les échantillons
        for sample_id in tqdm(sample_ids, desc="Traitement des échantillons", unit="échantillon"):
            generated_files = []

            callers = {
                "samtools": "_samtools_annot.vcf",
                "freebayes": "_Freebayes_annot.vcf",
                "HaplotypeCaller": "_gatk_annot.vcf",
                "vardict": "_vardict_annot.vcf"
            }

            for caller_name, suffix in callers.items():
                input_vcf = Path(input_dir) / f"{sample_id}{suffix}"
                if not input_vcf.exists():
                    global_log.write(f"[WARN] Fichier manquant : {input_vcf}\n")
                    continue

                output_vcf = Path(output_dir) / f"{sample_id}_{caller_name}_vartype.vcf"
                cmd = ["java", "-jar", str(snpsift_jar), "varType", str(input_vcf)]

                global_log.write(f"[CMD] {' '.join(cmd)}\n")

                try:
                    # Capture sortie complète de SnpSift
                    result = subprocess.run(
                        cmd,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        text=True,
                        check=True
                    )
                    with open(output_vcf, "w") as f_out:
                        f_out.write(result.stdout)
                    generated_files.append(str(output_vcf))
                    global_log.write(f"[OK] {output_vcf} généré avec succès.\n\n")
                except subprocess.CalledProcessError as e:
                    global_log.write(f"[ERROR] Échec sur {sample_id} - {caller_name}\n{e.output}\n\n")

            results[sample_id] = generated_files

        total_time = time.time() - start_time
        log_summary = f"\n[PIPELINE] Étape VarType terminée en {total_time/60:.2f} minutes.\n"
        print(log_summary)
        global_log.write(log_summary)
        global_log.write(f"Fin : {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        global_log.write("[PIPELINE] --- FIN DE L'ÉTAPE VARTYPE ---\n")
        print(f"[INFO] Log global disponible ici : {log_file}")

    return results


def get_coverage(bam_dir="output/bam_picard_readgroups", output_dir="output/samtools_coverage", logs_dir="logs"):
    """
    Calcule la couverture et la profondeur pour chaque échantillon à partir des fichiers BAM.

    Étapes :
      - Indexation du fichier BAM
      - Calcul de la couverture avec samtools coverage
      - Calcul de la profondeur avec samtools depth
      - Génération d’un fichier combiné finalcoverageall.txt

    Paramètres :
        bam_dir (str): dossier contenant les fichiers *_picard_readgroup.bam
        output_dir (str): dossier de sortie pour les fichiers de couverture
        logs_dir (str): dossier pour enregistrer les fichiers logs
    """
    start_time = time.time()

    bam_dir = Path(bam_dir)
    output_dir = Path(output_dir)
    logs_dir = Path(logs_dir)

    output_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    log_file = logs_dir / "get_coverage.log"

    bam_files = sorted(bam_dir.glob("*_picard_readgroup.bam"))
    n_samples = len(bam_files)

    if not bam_files:
        print("[WARNING] Aucun fichier BAM trouvé dans", bam_dir)
        return

    print(f"[PIPELINE] Démarrage de l’étape get_coverage : {n_samples} fichiers détectés")

    with open(log_file, "a") as lf:
        for bam in tqdm(bam_files, desc="Calcul de couverture", unit="échantillon"):
            sample_id = bam.stem.replace("_picard_readgroup", "")
            sample_dir = output_dir / sample_id
            sample_dir.mkdir(parents=True, exist_ok=True)

            coverage_txt = sample_dir / f"{sample_id}_coverage.txt"
            depth_txt = sample_dir / f"{sample_id}_depth.txt"
            final_cov = sample_dir / f"{sample_id}_finalcoverageall.txt"

            try:
                # 1. Indexation du BAM
                subprocess.run(["samtools", "index", str(bam)], stdout=lf, stderr=lf, text=True, check=True)

                # 2. Calcul de la couverture
                subprocess.run(
                    ["samtools", "coverage", str(bam), "-o", str(coverage_txt)],
                    stdout=lf, stderr=lf, text=True, check=True)

                # 3. Calcul de la profondeur
                subprocess.run(
                    ["samtools", "depth", str(bam), "-aa", "-o", str(depth_txt)],
                    stdout=lf, stderr=lf, text=True, check=True)

                # 4. Fusion des fichiers coverage
                with open(final_cov, "w") as f_out, open(coverage_txt) as f_in:
                    f_out.write(f_in.read())

                print(f"[INFO] Couverture générée pour {sample_id}")
                lf.write(f"[SUCCESS] Étape terminée pour {sample_id}\n")

            except subprocess.CalledProcessError as e:
                print(f"[ERROR] Échec de get_coverage pour {sample_id}")
                lf.write(f"[ERROR] Échec de get_coverage pour {sample_id} : {e}\n")
                continue

    elapsed_time = time.time() - start_time
    print(f"[PIPELINE] Étape get_coverage terminée en {elapsed_time:.2f} secondes.")
    print(f"[LOG] Détails complets enregistrés dans : {log_file}")


def run_wt_cov(ref, gff, voi, depth_dir="output/samtools_coverage", output_dir="output/WT_cov", logs_dir="logs"):
    """
    Exécute wt_cov.py automatiquement pour tous les fichiers *_depth.txt trouvés dans depth_dir.
    Les fichiers résultats (_coverage.csv) sont déplacés dans output/WT_cov/.
    Un seul log global est créé dans logs/.
    """

    depth_dir = Path(depth_dir).resolve()
    output_dir = Path(output_dir).resolve()
    logs_dir = Path(logs_dir).resolve()

    # Création des dossiers si nécessaires
    output_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    log_file = logs_dir / "WT_cov.log"
    script_path = Path("bin/wt_cov.py").resolve()

    # Trouver tous les fichiers *_depth.txt
    depth_files = sorted(depth_dir.rglob("*_depth.txt"))
    if not depth_files:
        print(f"[ERROR] Aucun fichier *_depth.txt trouvé dans {depth_dir}")
        return {}

    print(f"[PIPELINE] Démarrage de l’étape WT_cov : {len(depth_files)} échantillons détectés\n") 

    results = {}
    start_time = time.time()

    with open(log_file, "w") as log:
        log.write(f"[PIPELINE] WT_cov démarré le {time.ctime()}\n")
        log.write(f"Nombre de fichiers détectés : {len(depth_files)}\n\n")

        for depth_file in tqdm(depth_files, desc="WT_cov en cours", unit="échantillon"):
            sample_id = depth_file.stem.replace("_depth", "")
            expected_name = f"{sample_id}_coverage.csv"
            output_file = output_dir / expected_name

            cmd = [
                "python3", str(script_path),
                "-R", ref,
                "-G", gff,
                "-N", str(depth_file),
                "-V", voi
            ]

            try:
                # Exécuter le script
                subprocess.run(cmd, stdout=log, stderr=log, text=True, check=True)

                # Chercher le fichier généré (par défaut dans le même dossier que pipeline)
                local_csv = Path(expected_name)
                script_csv = script_path.parent / expected_name

                if local_csv.exists():
                    shutil.move(str(local_csv), str(output_file))
                elif script_csv.exists():
                    shutil.move(str(script_csv), str(output_file))

                results[sample_id] = str(output_file)
                log.write(f"[SUCCESS] {sample_id} terminé -> {output_file}\n")
                print(f"[SUCCESS] {sample_id} terminé")

            except subprocess.CalledProcessError as e:
                log.write(f"[ERROR] WT_cov échoué pour {sample_id} : {e}\n")
                print(f"[ERROR] WT_cov échoué pour {sample_id} : {e}")

        log.write("\n[PIPELINE] WT_cov terminé le {}\n".format(time.ctime()))

    elapsed = time.time() - start_time
    print(f"[PIPELINE] Étape WT_cov terminée en {elapsed:.2f} secondes pour {len(results)} échantillons.\n")
    print(f"[LOG] Fichier global : {log_file}\n")

    return results



def run_trim_stats(stats_dir="output/trimmed_reads", coverage_dir="output/samtools_coverage",   output_dir="output/Readscoverage", logs_dir="logs"):
    """
    Exécute reads_stats.py pour tous les échantillons détectés.
    - Les fichiers .stats.txt sont dans output/trimmed_reads
    - Les coverage.txt correspondants dans output/samtools_coverage/<echX>/
    - Produit *_readcoverage.csv dans output/Readscoverage
    - Tous les logs détaillés sont redirigés vers logs/Trim_Stats.log
    """

    stats_dir = Path(stats_dir).resolve()
    coverage_dir = Path(coverage_dir).resolve()
    output_dir = Path(output_dir).resolve()
    logs_dir = Path(logs_dir).resolve()

    output_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    log_file = logs_dir / "Trim_Stats.log"
    script_path = Path("bin/reads_stats.py").resolve()

    stats_files = sorted(stats_dir.glob("*.stats.txt"))
    if not stats_files:
        print(f"[ERREUR] Aucun fichier .stats.txt trouvé dans {stats_dir}")
        return {}

    results = {}
    start_time = time.time()

    print(f"[PIPELINE] Démarrage de l’étape Trim_Stats : {len(stats_files)} échantillons détectés\n")

    with open(log_file, "w") as log:
        log.write(f"[PIPELINE] Trim_Stats démarré le {time.ctime()}\n")
        log.write(f"Nombre d'échantillons : {len(stats_files)}\n\n")

        for stats_file in tqdm(stats_files, desc="Analyse Trim_Stats", unit="échantillon", ncols=80):
            sample_id = stats_file.stem.replace(".stats", "")
            cov_file = coverage_dir / sample_id / f"{sample_id}_coverage.txt"
            output_file = output_dir / f"{sample_id}_readcoverage.csv"

            if not cov_file.exists():
                msg = f"[WARNING] Fichier coverage introuvable pour {sample_id} : {cov_file}\n"
                log.write(msg)
                continue

            cmd = [
                "python3", str(script_path),
                "-S", str(stats_file),
                "-N", str(cov_file),
                "-o", str(output_dir)
            ]

            # Redirection complète stdout + stderr vers le log
            process = subprocess.run(cmd, stdout=log, stderr=log, text=True)

            if process.returncode == 0:
                results[sample_id] = str(output_file)
                log.write(f"[SUCCESS] {sample_id} terminé avec succès\n")
            else:
                log.write(f"[ERROR] Échec pour {sample_id} (code {process.returncode})\n")

        elapsed = time.time() - start_time
        log.write(f"\n[PIPELINE] Trim_Stats terminé le {time.ctime()}\n")
        log.write(f"Durée totale : {elapsed:.2f} secondes\n")

    # Terminal propre : résumé minimal
    print(f"\n[PIPELINE] Étape Trim_Stats terminée ({len(results)} échantillons traités)")
    print(f"[LOG] Tous les détails enregistrés dans : {log_file}\n")

    return results


def run_reads_merge(input_dir="output/Readscoverage", output_dir="output/Summary", logs_dir="logs"):
    """
    Fusionne tous les fichiers *_readcoverage.csv en un seul fichier Reads_Metrics_Samples.csv.
    - Supprime les en-têtes dupliqués (comme awk)
    - Sauvegarde dans output/Summary
    - Génère un log complet et affiche une progression dans le terminal
    """

    input_dir = Path(input_dir).resolve()
    output_dir = Path(output_dir).resolve()
    logs_dir = Path(logs_dir).resolve()

    output_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    log_file = logs_dir / "Reads_merge.log"
    output_file = output_dir / "Reads_Metrics_Samples.csv"

    csv_files = sorted(input_dir.glob("*.csv"))

    if not csv_files:
        print(f"[WARNING] Aucun fichier CSV trouvé dans {input_dir}")
        return None

    start_time = time.time()

    with open(log_file, "w") as log:
        log.write(f"[PIPELINE] Étape Reads_merge démarrée le {time.ctime()}\n")
        log.write(f"Nombre de fichiers détectés : {len(csv_files)}\n\n")

        df_list = []

        for file in tqdm(csv_files, desc="Fusion des fichiers CSV", unit="fichier"):
            try:
                df = pd.read_csv(file)
                df_list.append(df)
                log.write(f"[OK] Fichier fusionné : {file.name}\n")
            except Exception as e:
                log.write(f"[ERROR] Erreur lors de la lecture de {file.name} : {e}\n")

        # Fusion finale
        merged_df = pd.concat(df_list, ignore_index=True)
        merged_df.to_csv(output_file, index=False)

        elapsed = time.time() - start_time
        log.write(f"\n[PIPELINE] Reads_merge terminé le {time.ctime()} en {elapsed:.2f} secondes\n")
        log.write(f"[OUTPUT] Fichier final : {output_file}\n")

    print(f"\n[PIPELINE] Étape Reads_merge terminée en {elapsed:.2f} secondes.")
    print(f"[LOG] Fichier de log : {log_file}")
    print(f"[OUTPUT] Résultat : {output_file}")

    return output_file




def run_vcf_to_df(input_dir="output/vartype", output_dir="output/vcf_to_DF", script_path="bin/vcf_merge.py", logs_dir="logs"):
    """
    Exécute vcf_merge.py sur tous les échantillons détectés dans output/vartype/.
    - Chaque échantillon doit avoir 4 fichiers *_vartype.vcf
    - Crée un sous-dossier par sample dans output/vcf_to_DF/
    - Écrase les anciens résultats à chaque exécution
    - Log global dans logs/vcf_to_DF.log
    """

    # --- Préparation des chemins ---
    input_dir = Path(input_dir).resolve()
    output_dir = Path(output_dir).resolve()
    logs_dir = Path(logs_dir).resolve()
    script_path = Path(script_path).resolve()

    logs_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    log_file = logs_dir / "vcf_to_DF.log"

    # --- Détection des échantillons ---
    sample_ids = sorted({
    re.sub(r"_(samtools|freebayes|HaplotypeCaller|vardict)_vartype\.vcf$", "", f.name)
    for f in input_dir.glob("*_vartype.vcf")})


    if not sample_ids:
        print(f"[ERROR] Aucun échantillon détecté dans {input_dir}")
        return {}

    print(f"[PIPELINE] Étape vcf_to_DF : {len(sample_ids)} échantillons détectés\n")

    results = {}
    global_start = time.time()

    with open(log_file, "w") as log:
        log.write(f"===== Étape vcf_to_DF démarrée : {time.ctime()} =====\n")
        log.write(f"Échantillons : {', '.join(sample_ids)}\n\n")

        for sample_id in tqdm(sample_ids, desc="vcf_to_DF en cours", unit="échantillon"):
            start_sample = time.time()
            sample_output = output_dir / sample_id

            # Écrasement des anciens résultats
            if sample_output.exists():
                shutil.rmtree(sample_output)
            sample_output.mkdir(parents=True, exist_ok=True)

            # Vérification des fichiers VCF nécessaires
            vcf_files = {
                "samtools": input_dir / f"{sample_id}_samtools_vartype.vcf",
                "freeBayes": input_dir / f"{sample_id}_freebayes_vartype.vcf",
                "HaplotypeCaller": input_dir / f"{sample_id}_HaplotypeCaller_vartype.vcf",
                "vardict": input_dir / f"{sample_id}_vardict_vartype.vcf",
            }

            missing = [tool for tool, path in vcf_files.items() if not path.exists()]
            if missing:
                msg = f"[WARNING] {sample_id} : fichiers manquants {missing}\n"
                print(msg.strip())
                log.write(msg)
                continue

            # --- Exécution du script pour chaque outil ---
            for tool, vcf in vcf_files.items():
                cmd = ["python3", str(script_path), "-n", str(vcf), "-o", str(sample_output)]
                try:
                    subprocess.run(cmd, cwd=sample_output, stdout=log, stderr=log, text=True, check=True)
                except subprocess.CalledProcessError as e:
                    log.write(f"[ERREUR] {sample_id} - {tool} a échoué : {e}\n")
                    continue

            elapsed = time.time() - start_sample
            log.write(f"[SUCCESS] {sample_id} terminé en {elapsed:.2f} sec\n")
            results[sample_id] = sample_output

        total_elapsed = time.time() - global_start
        log.write(f"\n===== Étape vcf_to_DF terminée en {total_elapsed:.2f} sec =====\n")

    print(f"\n Étape vcf_to_DF terminée avec succès en {total_elapsed:.2f} sec.")
    print(f"[LOG] Détails : {log_file}\n")

    return results




def run_csv_merge_with_report(input_dir="output/vcf_to_DF", output_dir="output/CSV_merge", script_path="bin/csv_merge.py", logs_dir="logs"):
    """
    Étape CSV_merge améliorée :
    - Fusionne les CSV (samtools, freebayes, HC, vardict)
    - Génère un rapport de complétude pour chaque échantillon
    - Continue même si certains fichiers sont manquants
    """

    input_dir = Path(input_dir).resolve()
    output_dir = Path(output_dir).resolve()
    logs_dir = Path(logs_dir).resolve()
    script_path = Path(script_path).resolve()

    output_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    log_file = logs_dir / "CSV_merge.log"
    report_file = output_dir / "csv_merge_report.csv"

    # Détection des échantillons
    sample_dirs = [d for d in input_dir.iterdir() if d.is_dir()]
    sample_ids = [d.name for d in sample_dirs]

    if not sample_ids:
        print(f"[ERROR] Aucun dossier d'échantillon trouvé dans {input_dir}")
        return {}

    results = {}
    report_lines = []
    global_start = time.time()

    print(f"[PIPELINE] Étape CSV_merge : {len(sample_ids)} échantillons détectés\n")

    with open(log_file, "w") as log:
        log.write(f"\n===== Étape CSV_merge démarrée : {time.ctime()} =====\n")
        log.write(f"Échantillons : {', '.join(sample_ids)}\n\n")

        for sample_id in tqdm(sample_ids, desc="CSV_merge en cours", unit="échantillon"):
            sample_start = time.time()
            sample_input = input_dir / sample_id
            sample_output = output_dir / sample_id
            sample_output.mkdir(parents=True, exist_ok=True)

            # fichiers attendus
            csv_files = {
                "samtools": sample_input / f"{sample_id}_samtools.csv",
                "freebayes": sample_input / f"{sample_id}_freebayes.csv",
                "HaplotypeCaller": sample_input / f"{sample_id}_HaplotypeCaller.csv",
                "vardict": sample_input / f"{sample_id}_vardict.csv"
            }

            missing = [tool for tool, path in csv_files.items() if not path.exists()]
            status = "OK" if not missing else "INCOMPLETE"
            report_lines.append(f"{sample_id},{'|'.join(missing) if missing else 'None'},{status}")

            if missing:
                msg = f"[WARNING] Fichiers CSV manquants pour {sample_id} : {missing}\n"
                print(msg.strip())
                log.write(msg)

            # Commande pour les fichiers existants uniquement
            cmd = ["python3", str(script_path)]
            if csv_files["samtools"].exists(): cmd += ["-v1", str(csv_files["samtools"])]
            if csv_files["freebayes"].exists(): cmd += ["-v2", str(csv_files["freebayes"])]
            if csv_files["HaplotypeCaller"].exists(): cmd += ["-v3", str(csv_files["HaplotypeCaller"])]
            if csv_files["vardict"].exists(): cmd += ["-v4", str(csv_files["vardict"])]
            cmd += ["-o", str(sample_output)]

            # Exécution du merge
            subprocess.run(cmd, stdout=log, stderr=log, text=True, check=False)

            merge_csv = sample_output / f"{sample_id}_merge.csv"
            introns_csv = sample_output / f"{sample_id}_introns.csv"

            results[sample_id] = {"merge": merge_csv, "introns": introns_csv}

            elapsed = time.time() - sample_start
            log.write(f"[DONE] {sample_id} terminé en {elapsed:.2f} sec\n")

        total_elapsed = time.time() - global_start
        log.write(f"\n===== Étape CSV_merge terminée en {total_elapsed:.2f} sec =====\n")

    # Génération du rapport de complétude
    with open(report_file, "w") as rep:
        rep.write("Sample_ID,Missing_Files,Status\n")
        rep.write("\n".join(report_lines))
    print(f"\n[REPORT] Rapport de complétude généré : {report_file}")

    # Résumé global
    nb_ok = sum(1 for line in report_lines if line.endswith("OK"))
    nb_incomplete = len(sample_ids) - nb_ok
    print(f"[SUMMARY] {nb_ok} échantillons complets, {nb_incomplete} incomplets sur {len(sample_ids)}.")

    print("\n[PIPELINE] Étape CSV_merge terminée avec succès.")
    print(f"[LOG] Détails : {log_file}")

    return results


################################################################################################################
# Calcule des Sample_VAF à partir des fichiers merge générés précédemment
################################################################################################################

# ======================================================
# Identifier la taille du pool à partir du Sample_name
# ======================================================
def get_pool_size(sample_name):
    match = re.search(r'P(\d+)', sample_name)
    if match:
        return int(match.group(1))
    return 1


# ======================================================
# Extraire les métadonnées depuis le Sample_name
# ======================================================
def parse_sample_metadata(sample_name):
    return {
        "YEAR": f"20{sample_name[0:2]}",
        "COUNTRY": sample_name[2:4],
        "SITE": sample_name[4:6],
        "TREATMENT_DAY": sample_name[6:8]
    }


# ======================================================
# Calcul du Sample_VAF pour tous les fichiers merge
# ======================================================
def compute_sample_vaf(input_dir="output/CSV_merge", output_dir="output/Sample_VAF"):
       
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    merge_files = glob.glob(str(input_dir / "*" / "*_merge.csv"))

    if not merge_files:
        raise FileNotFoundError("Aucun fichier *_merge.csv trouvé")

    for file in merge_files:
        print(f"[PROCESS] {file}")

        df = pd.read_csv(file)

        sample_name = Path(file).stem.replace("_merge", "")
        pool_size = get_pool_size(sample_name)
        metadata = parse_sample_metadata(sample_name)

        # POOL status et poids
        pool_status = "pooled" if pool_size > 1 else "individual"
        weight = pool_size

        # Calcul du Sample_VAF
        if pool_size > 1:
            df['SVAF'] = df['AVG_VAF'] / pool_size
        else:
            df['SVAF'] = df['AVG_VAF']

        # Sélection et ordre strict des colonnes
        final_df = df[[
            'CHROM',
            'POS',
            'REF',
            'ALT',
            'AA_change',
            'codon',
            'Confidence',
            'Annotation',
            'AVG_COV'
        ]].copy()

        # Ajout des métadonnées
        final_df.insert(0, 'Sample_name', sample_name)
        final_df.insert(1, 'COUNTRY', metadata['COUNTRY'])
        final_df.insert(2, 'SITE', metadata['SITE'])
        final_df.insert(3, 'YEAR', metadata['YEAR'])
        final_df.insert(4, 'TREATMENT_DAY', metadata['TREATMENT_DAY'])
        final_df.insert(5, 'POOL', pool_status)
        final_df.insert(6, 'Weight', weight)

        # Ajout du Sample_VAF
        final_df['SVAF'] = df['SVAF'].round(4)

        # Sauvegarde
        out_file = output_dir / f"{sample_name}_SVAF.csv"
        final_df.to_csv(out_file, index=False)

        print(f"[OK] Sauvegardé : {out_file}")

    print("\n[PIPELINE] Calcul des Sample_VAF terminé.\n")





# ======================================================
# 1. Charger tous les fichiers Sample_VAF
# ======================================================
def load_all_sample_vaf(sample_vaf_dir):
    files = glob.glob(str(Path(sample_vaf_dir) / "*_SVAF.csv"))
    if not files:
        raise FileNotFoundError("Aucun fichier *_SVAF.csv trouvé")

    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    return df

# ======================================================
# 2. Préparer le dataframe sans calcul du VAF
# ======================================================
def prepare_vaf_data(df):
    df = df.copy()
    
    # S'assurer que la colonne Weight existe (sinon par défaut = 1)
    if 'Weight' not in df.columns:
        df['Weight'] = 1

    # Recodage Type
    if 'Type' in df.columns:
        df['Type'] = df['Type'].replace({
            'missense_variant': 'mutant',
            'synonymous_variant': 'wild-type'
        })
    
    # Supprimer les colonnes liées au calcul du VAF si elles existent
    cols_to_drop = [
        'SVAF_weighted', 'VAF_FINAL_SNP', '%VAF_FINAL_SNP', 'TOTAL_WEIGHT', 'num_mut', 'num_mix', 
        'den_mut', 'den_mix'
    ]
    for col in cols_to_drop:
        if col in df.columns:
            df.drop(columns=[col], inplace=True)
    
    return df

# ======================================================
# 3. Script principal simplifié
# ======================================================
def run_SVAF_merge(input_sample_vaf_dir="output/Sample_VAF",
                        output_file="output/SVAF_merge/SVAF_merge.csv"):
    
    print("[PIPELINE] Chargement des Sample_VAF...")
    df = load_all_sample_vaf(input_sample_vaf_dir)
    
    print("[PIPELINE] Préparation des données sans VAF...")
    df_prepared = prepare_vaf_data(df)
    
    # Création du dossier de sortie si nécessaire
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Sauvegarde
    df_prepared.to_csv(output_path, index=False)
    print(f"[OK] Fichier préparé sauvegardé : {output_path}")
    print("\n[PIPELINE] Traitement terminé.\n")





def run_snpfilter(voi_path="pf_3D7_Ref/voinew3.csv", merge_dir="output/CSV_merge", coverage_dir="output/WT_cov",
                    output_dir="output/Snpfilter", script_path="bin/final_snpfilter.py", logs_dir="logs"):
    """
    Exécute final_snpfilter.py pour chaque échantillon.
    - VOI : fichier commun
    - merge.csv : output/CSV_merge/<sample_id>/<sample_id>_merge.csv
    - coverage.csv : output/WT_cov/<sample_id>_coverage.csv
    """

    voi_path = Path(voi_path).resolve()
    merge_dir = Path(merge_dir).resolve()
    coverage_dir = Path(coverage_dir).resolve()
    output_dir = Path(output_dir).resolve()
    logs_dir = Path(logs_dir).resolve()
    script_path = Path(script_path).resolve()

    output_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    # Écrase le log précédent
    log_file = logs_dir / "Snpfilter.log"
    with open(log_file, "w") as log:
        log.write(f"===== Étape Snpfilter démarrée : {time.ctime()} =====\n")

    sample_ids = [d.name for d in merge_dir.iterdir() if d.is_dir()]
    if not sample_ids:
        print(f"[ERROR] Aucun échantillon trouvé dans {merge_dir}")
        return {}

    results = {}
    global_start = time.time()

    print(f"[PIPELINE] Étape Snpfilter : {len(sample_ids)} échantillons détectés\n")

    with open(log_file, "a") as log:
        for sample_id in tqdm(sample_ids, desc="Snpfilter en cours", unit="échantillon"):
            sample_start = time.time()

            merge_csv = merge_dir / sample_id / f"{sample_id}_merge.csv"
            coverage_csv = coverage_dir / f"{sample_id}_coverage.csv"

            # Vérification des fichiers d'entrée
            if not merge_csv.exists():
                log.write(f"[WARNING] {sample_id} : Fichier merge absent : {merge_csv}\n")
                continue
            if not coverage_csv.exists():
                log.write(f"[WARNING] {sample_id} : Fichier coverage absent : {coverage_csv}\n")
                continue

            # Écrase les anciens résultats
            sample_output = output_dir / sample_id
            if sample_output.exists():
                shutil.rmtree(sample_output)
            sample_output.mkdir(parents=True, exist_ok=True)

            final_snp = sample_output / f"{sample_id}_final_snp.csv"

            cmd = [
                "python3", str(script_path),
                "-A", str(merge_csv),
                "-C", str(coverage_csv),
                "-V", str(voi_path)
            ]

            # Exécution silencieuse, erreurs dans le log uniquement
            try:
                with open(final_snp, "w") as out, open(log_file, "a") as log2:
                    subprocess.run(
                        cmd,
                        stdout=out,
                        stderr=log2,
                        text=True,
                        cwd=sample_output,
                        check=True
                    )
                elapsed = time.time() - sample_start
                log.write(f"[SUCCESS] {sample_id} terminé en {elapsed:.2f} sec\n")
                results[sample_id] = final_snp

            except subprocess.CalledProcessError as e:
                log.write(f"[ERROR] {sample_id} a échoué (code {e.returncode}). Voir détails ci-dessus.\n")
                continue

        total_elapsed = time.time() - global_start
        log.write(f"\n===== Étape Snpfilter terminée en {total_elapsed:.2f} sec =====\n")

    print("\n [PIPELINE] Étape Snpfilter terminée.")
    print(f"[LOG] Détails : {log_file}")

    return results





def run_summary_merge(input_dir="output/Snpfilter", output_dir="output/Summary_merge", logs_dir="logs"):
    """
    Fusionne tous les fichiers *_final_snp.csv issus de l'étape Snpfilter
    en un seul fichier merged_final_snp.csv (en supprimant les doublons d'en-tête).
    """

    input_dir = Path(input_dir).resolve()
    output_dir = Path(output_dir).resolve()
    logs_dir = Path(logs_dir).resolve()

    output_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    merged_file = output_dir / "merged_final_snp.csv"
    log_file = logs_dir / "Summary_merge.log"

    with open(log_file, "w") as log:
        log.write(f"===== Étape Summary_merge démarrée : {time.ctime()} =====\n")

    csv_files = sorted(input_dir.glob("*/**/*_final_snp.csv"))

    if not csv_files:
        print(f"[ERROR] Aucun fichier *_final_snp.csv trouvé dans {input_dir}")
        return None

    print(f"[PIPELINE] Étape Summary_merge : {len(csv_files)} fichiers à fusionner\n")

    merged_df = None
    global_start = time.time()

    for csv_file in tqdm(csv_files, desc="Fusion des CSV", unit="fichier"):
        try:
            df = pd.read_csv(csv_file)
            if merged_df is None:
                merged_df = df
            else:
                merged_df = pd.concat([merged_df, df.iloc[1:]], ignore_index=True)
        except Exception as e:
            with open(log_file, "a") as log:
                log.write(f"[ERROR] Impossible de lire {csv_file.name} : {e}\n")

    if merged_df is not None:
        merged_df.to_csv(merged_file, index=False)
        elapsed = time.time() - global_start
        with open(log_file, "a") as log:
            log.write(f"[SUCCESS] Fichier fusionné généré : {merged_file}\n")
            log.write(f"Durée totale : {elapsed:.2f} sec\n")
        print(f"\n Fusion terminée : {merged_file}")
    else:
        print("[ERROR] Aucun fichier valide fusionné.")

    print(f"[LOG] Détails : {log_file}")
    return merged_file



def run_introns_merge(input_dir="output/CSV_merge", output_dir="output/Summary", logs_dir="logs"):
    """
    Fusionne tous les fichiers *_introns.csv des sous-dossiers dans output/CSV_merge/
    et crée un fichier unique : Introns_final_snp.csv
    """
    input_dir = Path(input_dir).resolve()
    output_dir = Path(output_dir).resolve()
    logs_dir = Path(logs_dir).resolve()

    output_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / "Introns_final_snp.csv"
    log_file = logs_dir / "Introns_merge.log"

    start = time.time()

    # Initialiser le log (en écrasant le précédent)
    with open(log_file, "w") as log:
        log.write(f"===== Étape Introns_merge démarrée : {time.ctime()} =====\n")

    introns_files = sorted(Path(f).resolve() for f in input_dir.glob("*/**/*_introns.csv"))

    if not introns_files:
        print(f"[ERROR] Aucun fichier *_introns.csv trouvé dans {input_dir}")
        return None

    merged_df = pd.DataFrame()

    # Capture des warnings pandas et envoi vers le fichier log
    with warnings.catch_warnings(record=True) as caught_warnings:
        warnings.simplefilter("always")

        for f in tqdm(introns_files, desc="Fusion des fichiers introns", unit="fichier"):
            try:
                df = pd.read_csv(f)
                if not df.empty and not df.dropna(how="all").empty:
                    if merged_df.empty:
                        merged_df = df
                    else:
                        merged_df = pd.concat([merged_df, df.iloc[1:]], ignore_index=True)
            except Exception as e:
                with open(log_file, "a") as log:
                    log.write(f"[ERROR] Impossible de lire {f.name} : {e}\n")

        # Écrire les warnings pandas dans le log
        with open(log_file, "a") as log:
            for w in caught_warnings:
                log.write(f"[WARNING] {w.message}\n")

    # Écriture du fichier final fusionné
    merged_df.to_csv(output_file, index=False)

    elapsed = time.time() - start
    with open(log_file, "a") as log:
        log.write(f"\n Fusion terminée : {output_file}\n")
        log.write(f"Durée totale : {elapsed:.2f} sec\n")
        log.write("===== Étape Introns_merge terminée =====\n")

    print(f"\n Fusion terminée : {output_file}")
    print(f"[LOG] Détails : {log_file}")

    return output_file




def run_summary(merged_final_snp="output/Summary_merge/merged_final_snp.csv", output_dir="output/Summary",
                script_path="bin/summary.py", logs_dir="logs"):
    """
    Équivalent Python du process Nextflow 'Summary'.
    Exécute summary.py sur merged_final_snp.csv et publie les fichiers :
    - All_final_snp.csv
    - Reportable_snps.csv
    - Novel_snps.csv
    Tous les fichiers générés seront placés dans output_dir.
    """

    # Préparation des chemins
    merged_final_snp = Path(merged_final_snp).resolve()
    output_dir = Path(output_dir).resolve()
    logs_dir = Path(logs_dir).resolve()
    script_path = Path(script_path).resolve()

    output_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    log_file = logs_dir / "Summary.log"
    start_time = time.time()

    with open(log_file, "w") as log:
        log.write(f"===== Étape Summary démarrée : {time.ctime()} =====\n")

        if not merged_final_snp.exists():
            msg = f"[ERREUR] Fichier introuvable : {merged_final_snp}\n"
            print(msg)
            log.write(msg)
            return

        print(f"[PIPELINE] Étape Summary sur : {merged_final_snp.name}\n")

        cmd = ["python3", str(script_path), "-f", str(merged_final_snp)]

        try:
            # Exécution dans le dossier de sortie
            result = subprocess.run(cmd, cwd=output_dir, capture_output=True, text=True, check=True)

            log.write("[STDOUT]\n" + result.stdout + "\n")
            log.write("[STDERR]\n" + result.stderr + "\n")

            # Vérification et rapport des fichiers générés
            expected_files = [
                "All_final_snp.csv",
                "Reportable_snps.csv",
                "Novel_snps.csv"
            ]

            for f in tqdm(expected_files):
                f_path = output_dir / f
                if f_path.exists():
                    log.write(f"[SUCCESS] Fichier généré : {f_path}\n")
                else:
                    log.write(f"[WARNING] Fichier manquant : {f_path}\n")

            elapsed = time.time() - start_time
            log.write(f"\n===== Étape Summary terminée en {elapsed:.2f} sec =====\n")

            print(" Étape Summary terminée avec succès.")
            print(f"[LOG] Détails : {log_file}")

        except subprocess.CalledProcessError as e:
            log.write(f"\n[ERREUR] L'exécution de summary.py a échoué : {e}\n")
            log.write(f"STDOUT:\n{e.stdout}\nSTDERR:\n{e.stderr}\n")
            print(" Erreur lors de l'exécution de summary.py. Voir le fichier log pour les détails.")
            
            


def run_dataviz_reportable_snps(input_dir="output/Summary", out_dir="output/Dataviz_Reportable_snps",
                                script_path="bin/Dataviz_Reportable_snps.py", logs_dir="logs"):
    """
    Reproduction du process Nextflow 'Dataviz_Reportable_snps' en Python.
    - Parcourt les fichiers CSV dans input_dir.
    - Exécute 'Dataviz_Reportable_snps.py -r <fichier>'.
    - Sauvegarde tout dans logs/Dataviz_Reportable_snps.log.
    - Efface et recrée le dossier de sortie à chaque exécution.
    - N'affiche que la progression, le timer et les prints.
    """

    # --- Préparation des dossiers ---
    input_dir = Path(input_dir).resolve()
    out_dir = Path(out_dir).resolve()
    logs_dir = Path(logs_dir).resolve()
    script_path = Path(script_path).resolve()

    logs_dir.mkdir(parents=True, exist_ok=True)

    log_file = logs_dir / "Dataviz_Reportable_snps.log"

    # Efface le dossier de sortie avant chaque exécution 
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Efface l'ancien log
    if log_file.exists():
        log_file.unlink()

    # --- Recherche des fichiers d'entrée ---
    input_files = sorted(input_dir.glob("*.csv"))
    if not input_files:
        print(f"[ERREUR] Aucun fichier CSV trouvé dans {input_dir}")
        return

    print(f"\n Étape Dataviz_Reportable_snps : {len(input_files)} fichier(s) détecté(s)")
    start_time = time.time()

    with open(log_file, "w") as log:
        log.write(f"===== Début Dataviz_Reportable_snps : {time.ctime()} =====\n")
        log.write(f"Input_dir : {input_dir}\n")
        log.write(f"Output_dir : {out_dir}\n")
        log.write(f"Script : {script_path}\n\n")

        for input_file in tqdm(input_files, desc="Dataviz_Reportable_snps", unit="fichier"):
            sample_start = time.time()
            log.write(f"--- Traitement du fichier : {input_file.name} ---\n")

            cmd = ["python3", str(script_path), "-r", str(input_file)]

            try:
                process = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    check=False,
                )

                log.write("=== STDOUT ===\n" + (process.stdout or "") + "\n")
                log.write("=== STDERR ===\n" + (process.stderr or "") + "\n")

                if process.returncode != 0:
                    log.write(f"[WARNING] {input_file.name} : erreur ignorée (code {process.returncode})\n")

                # 🔹 Récupération automatique des fichiers générés
                expected_outputs = [
                    "DMS_EPI_report.csv",
                    "Reportable_Per_SNP_depth.pdf",
                    "SNPs-Reportable.pdf",
                ]

                for file_name in expected_outputs:
                    found_file = Path(file_name)
                    if found_file.exists():
                        dest = out_dir / f"{input_file.stem}_{file_name}"
                        shutil.move(str(found_file), dest)
                        log.write(f"[OK] Fichier déplacé : {dest}\n")

            except Exception as e:
                log.write(f"[ERREUR] Exception pendant {input_file.name} : {e}\n")

            elapsed = time.time() - sample_start
            log.write(f"[OK] {input_file.name} terminé en {elapsed:.2f} sec\n\n")

        total_time = time.time() - start_time
        log.write(f"\n===== Étape terminée en {total_time:.2f} secondes =====\n")

    print(f"\n Étape Dataviz_Reportable_snps terminée en {total_time:.2f} sec.")
    print(f" Log disponible dans : {log_file}")
    print(f" Résultats enregistrés dans : {out_dir}")




def run_DataViz_Novel_snps(input_dir="output/Summary/Novel_snps.csv", out_dir="output/Dataviz_Novel_snps",
                            script_path="bin/Dataviz_Novel_snps.py", logs_dir="logs"):
    """
    Exécute le script Dataviz_Novel_snps.py sur le fichier spécifié,
    sauvegarde les logs dans logs/Dataviz_Novel_snps.log,P1	Pfdhps	1486	5547	71.0%	A437G	SNP	Mutant

    et copie les fichiers générés dans output/Novel_snps.
    """

    # Initialisation des chemins
    input_path = Path(input_dir)
    output_path = Path(out_dir)
    logs_path = Path(logs_dir)
    log_file = logs_path / "Dataviz_Novel_snps.log"

    # Création des dossiers nécessaires
    output_path.mkdir(parents=True, exist_ok=True)
    logs_path.mkdir(parents=True, exist_ok=True)

    # Début du timer
    start_time = time.time()
    print(f"\n Exécution du processus DataViz_Novel_snps sur {input_path.name}")
    print("-----------------------------------------------------------")

    # Préparation du fichier log (mode écrasement)
    with open(log_file, "w") as log:
        log.write(f"=== DataViz_Novel_snps Execution Log ===\n")
        log.write(f"Date : {datetime.now()}\n")
        log.write(f"Script : {script_path}\n")
        log.write(f"Entrée : {input_path}\n")
        log.write(f"Sortie : {output_path}\n\n")

        # Commande
        cmd = ["python3", script_path, "-n", str(input_path)]

        try:
            # Exécution du script avec redirection des flux vers le log
            process = subprocess.Popen(
                cmd,
                stdout=log,
                stderr=log
            )

            # Barre de progression simple pendant l’exécution
            while process.poll() is None:
                print(" En cours d'exécution...", end="\r")
                time.sleep(1)

            process.wait()
            print(" Exécution terminée.")

        except Exception as e:
            log.write(f"\n[ERREUR] {str(e)}\n")
            print(f" Une erreur est survenue : {e}")

    # Copie des résultats vers le dossier de sortie
    expected_outputs = [
        "SNPs-Novel-missense.pdf",
        "SNPs-Novel-synonymous.pdf"
    ]

    for filename in expected_outputs:
        file_path = Path(filename)
        dest_file = output_path / filename

        if file_path.exists():
            shutil.move(str(file_path), str(dest_file))
            print(f" Résultat copié : {dest_file}")
        else:
            with open(log_file, "a") as log:
                log.write(f"[AVERTISSEMENT] Fichier manquant : {filename}\n")

    # Temps total
    elapsed = time.time() - start_time
    print(f" Durée totale : {elapsed:.2f} secondes")
    print("-----------------------------------------------------------")
    print(f" Log disponible ici : {log_file}\n")





def run_replace_mutations(input_file="output/Dataviz_Reportable_snps/Reportable_snps_DMS_EPI_report.csv",
                            output_dir="output/haplotypes", logs_dir="logs"):
    """
    Parcourt le fichier CSV d'entrée et remplace les valeurs 'MT' ou 'MIX'
    par la dernière lettre du nom de la colonne correspondante et 'WT' par la première lettre.
    Le fichier original n'est pas modifié.
    """

    # === Création des dossiers ===
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)

    log_file = os.path.join(logs_dir, f"replace_mutations.log")

    with open(log_file, "w") as log:
        log.write(f"=== Replace Mutations Log ===\nDate : {datetime.now()}\n")
        log.write(f"Input : {input_file}\nOutput Dir : {output_dir}\n\n")

        try:
            # === Lecture du fichier CSV ===
            df = pd.read_csv(input_file, sep=",", dtype=str).fillna("NA")
            log.write(f"[INFO] Fichier chargé avec {df.shape[0]} lignes et {df.shape[1]} colonnes.\n")

            df_copy = df.copy()  # Ne pas modifier l'original

            # === Parcours des colonnes (sauf la première si elle contient des IDs) ===
            cols_to_process = [col for col in df_copy.columns if col != df_copy.columns[0]]

            for col in tqdm(cols_to_process, desc="Remplacement des mutations"):
                if not any(char.isdigit() for char in col):
                    continue  # ignorer les colonnes non génétiques (comme les # mutations)

                first_letter = col[0]
                last_letter = col[-1]

                # travailler sur une copie propre des valeurs textuelles
                s = df_copy[col].astype(str).str.strip()
                s_upper = s.str.upper()

                # remplacements insensibles à la casse : MT/MIX -> dernière lettre, WT -> première lettre
                s = s.mask(s_upper == "MT", last_letter)
                s = s.mask(s_upper == "MIX", last_letter)
                s = s.mask(s_upper == "WT", first_letter)

                df_copy[col] = s

            # === Sauvegarde dans un nouveau fichier ===
            output_file = os.path.join(output_dir, "Reportable_snps_replaced.csv")
            df_copy.to_csv(output_file, index=False)

            log.write("[INFO] Remplacements terminés avec succès.\n")
            log.write(f"Fichier généré : {output_file}\n")

        except Exception as e:
            log.write(f"[ERREUR] {str(e)}\n")
            raise

    print(f"\n Fichier modifié généré dans : {output_file}")
    print(f" Log enregistré dans : {log_file}")
    return output_file



    
    
def filtered_summary_merge(input_file="output/Summary_merge/merged_final_snp.csv", output_dir="output/Summary_merge_filtered", logs_dir="logs"):
    """
    Extrait et formate certaines colonnes du fichier merged_final_snp.csv :
      - Sélectionne les colonnes d’intérêt
      - Renomme les colonnes pour la cohérence du pipeline
      - Transforme les types de variants en 'Mutant' ou 'WildType'
      - Sauvegarde le résultat dans un fichier CSV filtré
    """

    # Création des dossiers si besoin
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    Path(logs_dir).mkdir(parents=True, exist_ok=True)
    log_file = Path(logs_dir) / "Summary_merge_filter.log"

    try:
        df = pd.read_csv(input_file)
    except Exception as e:
        raise FileNotFoundError(f"Erreur de lecture du fichier : {input_file}\n{e}")

    with open(log_file, "w") as log:
        log.write("\n=== Summary Filter Log ===\n")
        log.write(f"Entrée : {input_file}\n")
        log.write(f"Sortie : {output_dir}\n")

        # Colonnes nécessaires
        required_cols = ["Sample_name", "CHROM", "POS", "AVG_COV", "AVG_VAF", "VOI", "VARTYPE", "Annotation"]
        missing = [c for c in required_cols if c not in df.columns]
        if missing:
            log.write(f" Colonnes manquantes dans le fichier : {missing}\n")
            raise KeyError(f"Colonnes manquantes : {missing}")

        # ✅ Sélectionner uniquement les colonnes d’intérêt
        df_filtered = df[required_cols].copy()

        # ✅ Renommer les colonnes pour standardiser la sortie
        df_filtered.rename(columns={
            "Sample_name": "ID",
            "CHROM": "Gene",
            "AVG_COV": "Average_Coverage",
            "AVG_VAF": "Average_VAF(%)",
            "VOI": "Gene_Annotation",
            "Annotation": "Type"
        }, inplace=True)

        # ✅ Remplacer les valeurs de la colonne "Type"
        df_filtered["Type"] = df_filtered["Type"].replace({
            "missense_variant": "Mutant",
            "synonymous_variant": "WildType"
        })

        # ✅ Sauvegarde du résultat
        out_path = Path(output_dir) / "filtered_summary_merge.csv"
        df_filtered.to_csv(out_path, index=False)

        # ✅ Logs
        log.write(f" Fichier filtré généré : {out_path}\n")
        log.write(f" Nombre de lignes retenues : {len(df_filtered)}\n")
        log.write("=== Fin du traitement ===\n")

    print(f"✅ Résultat sauvegardé dans : {out_path}")
    print(f"🧾 Log disponible ici : {log_file}")

    


def filter_empty_pos(input_file = "output/Summary_merge_filtered/filtered_summary_merge.csv", 
                        output_file="output/Summary_merge_filtered/filtered_summary_clean.csv", logs_dir="logs"):

    """
    Copie le fichier d'entrée et supprime toutes les lignes
    où la colonne 'POS' est vide ou NaN.
    Le fichier d'entrée n'est jamais modifié.
    """

    # === Création des dossiers ===
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    Path(logs_dir).mkdir(parents=True, exist_ok=True)

    log_path = Path(logs_dir) / "filter_POS.log"

    # === Lecture du fichier source en toute sécurité ===
    df = pd.read_csv(input_file)

    # Comptage initial
    initial_rows = len(df)

    # === Filtrage : POS non vide ===
    df_filtered = df[df["POS"].notna() & (df["POS"].astype(str).str.strip() != "")]

    final_rows = len(df_filtered)

    # === Sauvegarde du fichier filtré ===
    df_filtered.to_csv(output_file, index=False)

    # === Log ===
    with open(log_path, "a") as log:
        log.write("\n=== Filter POS Log ===\n")
        # log.write(f"Date: {datetime.datetime.now()}\n")
        log.write(f"Input: {input_file}\n")
        log.write(f"Output: {output_file}\n")
        log.write(f"Total lignes initiales : {initial_rows}\n")
        log.write(f"Lignes conservées      : {final_rows}\n")
        log.write(f"Lignes supprimées      : {initial_rows - final_rows}\n")
        log.write("=== Fin ===\n")

    print(f"Fichier filtré généré : {output_file}")
    print(f"{initial_rows - final_rows} lignes supprimées (POS vide).")



def run_haplotypes_1(input_file="output/Dataviz_Reportable_snps/Reportable_snps_DMS_EPI_report.csv",
                    output_dir="output/haplotypes", logs_dir="logs"):
    """
    Construit des haplotypes pour chaque échantillon.
    Si une position est WT, la première lettre du SNP est utilisée.
    Si une mutation est présente, la dernière lettre du SNP est utilisée.
    """

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    Path(logs_dir).mkdir(parents=True, exist_ok=True)
    log_file = Path(logs_dir) / "Haplotypes.log"

    with open(log_file, "w") as log:
        log.write(f"\n=== Compute Haplotypes Log ===\nDate : {time.ctime()}\n")
        log.write(f"Entrée : {input_file}\nSortie : {output_dir}\n")

        df = pd.read_csv(input_file)
        log.write(f"Fichier lu avec {len(df)} lignes.\n")

        df['Sample'] = df['LSDB_Sequence_ID'].astype(str)
        
        
        # === Détection automatique des gènes selon les préfixes de colonnes ===
        gene_groups = {
            "DHPS": [c for c in df.columns if re.search(r'I431V|S436A|A437G|K540E|A581G|A613S|A613T', c)],
            "DHFR": [c for c in df.columns if re.search(r'N51I|C59R|S108N', c)],
            "CRT":  [c for c in df.columns if re.search(r'C72S|V73V|M74I|N75E|K76T|A220S|Q271E|N326S|C350R|R371I|I356T', c)],
            "MDR":  [c for c in df.columns if re.search(r'N86Y|Y184F|D1246Y|N1042D|S1034C', c)],
            "CytB": [c for c in df.columns if re.search(r'I258M|Y268C|Y268S', c)],
            "K13":  [c for c in df.columns if re.search(r'A481V|A578S|A675V|C469Y|C580Y|D584V|F446I|G449A|G538V|I543T|M476I|N458Y|N537I|P441L|P553L|P574L|R539T|R561H|V568G|Y493H', c)]
        }

        haplo_data = []

        for sample, group in tqdm(df.groupby("Sample"), desc="Construction des haplotypes"):
            sample_haplo = {"Sample": sample}

            for gene, positions in gene_groups.items():
                alleles = []

                for snp in positions:
                    # colonne correspondant au SNP exact
                    snp_cols = [col for col in df.columns if snp in col]
                    if not snp_cols:
                        continue

                    col = snp_cols[0]
                    values = group[col].dropna().unique()

                    if len(values) == 0:
                        alleles.append("")  # pas de donnée
                        continue

                    val = values[0]

                    if val == "WT":
                        # prendre la première lettre du SNP
                        alleles.append(snp[0])
                    elif val == "NA" or val == "MIX" :
                        alleles.append("")
                    else:
                        # mutation : dernière lettre
                        alleles.append(snp[-1])

                # construire haplotype par gène
                haplo = "".join(alleles)
                sample_haplo[gene] = haplo if haplo.strip("-") != "" else "Null(NAN)"

            haplo_data.append(sample_haplo)

        haplo_df = pd.DataFrame(haplo_data)
        out_path = Path(output_dir) / "haplotypes_summary.csv"
        haplo_df.to_csv(out_path, index=False)

        log.write(f"Haplotype généré : {out_path}\n")
        log.write("=== Fin ===\n")

    print(f"\n Haplotype généré avec succès : {out_path}")


def run_haplotypes(input_file="output/Dataviz_Reportable_snps/Reportable_snps_DMS_EPI_report.csv", 
                   output_dir="output/haplotypes", logs_dir="logs"):
    """
    Construit des haplotypes pour chaque échantillon.
    Règle stricte : si NA/MIX apparaît dans les X premières positions d'un gène → "Null".
    Lecture CSV forcée en str pour éviter que "NA" soit interprété comme NaN.
    """

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    Path(logs_dir).mkdir(parents=True, exist_ok=True)
    log_file = Path(logs_dir) / "Haplotypes.log"

    with open(log_file, "w") as log:
        log.write(f"\n=== Compute Haplotypes Log ===\nDate : {time.ctime()}\n")
        log.write(f"Entrée : {input_file}\nSortie : {output_dir}\n")

        # >>> IMPORTANT : forcer lecture en str et désactiver la conversion "NA" -> NaN
        df = pd.read_csv(input_file, dtype=str, keep_default_na=False, na_filter=False)
        log.write(f"Fichier lu avec {len(df)} lignes.\n")

        df['Sample'] = df['LSDB_Sequence_ID'].astype(str)

        gene_groups = {
            "DHPS": [c for c in df.columns if re.search(r'I431V|S436A|A437G|K540E|A581G|A613S|A613T', c)],
            "DHFR": [c for c in df.columns if re.search(r'N51I|C59R|S108N', c)],
            "CRT":  [c for c in df.columns if re.search(r'C72S|V73V|M74I|N75E|K76T|A220S|Q271E|N326S|C350R|R371I|I356T', c)],
            "MDR":  [c for c in df.columns if re.search(r'N86Y|Y184F|D1246Y|N1042D|S1034C', c)],
            "CytB": [c for c in df.columns if re.search(r'I258M|Y268C|Y268S', c)],
            "K13":  [c for c in df.columns if re.search(r'A481V|A578S|A675V|C469Y|C580Y|D584V|F446I|G449A|G538V|I543T|M476I|N458Y|N537I|P441L|P553L|P574L|R539T|R561H|V568G|Y493H', c)]
        }

        # règles NA/MIX : nombre de premières positions à vérifier
        NA_rules = {"DHPS": 5, "DHFR": 3, "CRT": 5, "MDR": 3}

        haplo_data = []

        for sample, group in tqdm(df.groupby("Sample"), desc="Construction des haplotypes"):
            sample_haplo = {"Sample": sample}

            for gene, positions in gene_groups.items():

                # ======== Étape 1 : Vérification NA/MIX sur les X premières colonnes ===========
                if gene in NA_rules and len(positions) > 0:
                    limit = NA_rules[gene]
                    first_cols = positions[:limit]

                    contains_NA_MIX = False
                    triggered = None

                    for col in first_cols:
                        if col not in group.columns:
                            continue

                        # Normaliser toutes les valeurs de la colonne : strip + upper
                        raw_vals = [str(x).strip().upper() for x in group[col].unique()]
                        # retirer vides éventuels
                        raw_vals = [v for v in raw_vals if v != ""]

                        if not raw_vals:
                            # si la colonne est vide après normalisation -> on considère qu'il n'y a pas de "NA" textuel
                            continue

                        # détecter NA/MIX même dans "MIX/WT" ou "NA "
                        for v in raw_vals:
                            #if v == "NA" or v == "MIX" or v.startswith("NA") or "MIX" in v:
                            if v == "NA" or v.startswith("NA"):
                                
                                contains_NA_MIX = True
                                triggered = (col, raw_vals)
                                break
                        if contains_NA_MIX:
                            break

                    if contains_NA_MIX:
                        sample_haplo[gene] = "Null"
                        log.write(f"{sample} | {gene} = Null (detecté: {triggered})\n")
                        continue  # pas de concaténation pour ce gène

                # ======== Étape 2 : Construction classique si pas de NA/MIX détecté ===========
                alleles = []

                for snp in positions:
                    snp_cols = [col for col in df.columns if snp in col]
                    if not snp_cols:
                        continue
                    col = snp_cols[0]

                    # Prenons la valeur telle quelle (déjà lue en str)
                    vals = [str(x).strip().upper() for x in group[col].unique()]
                    vals = [v for v in vals if v != ""]  # enlever vides

                    if not vals:
                        continue

                    v = vals[0]

                
                    if v == "WT":
                        alleles.append(snp[0])

                    elif v == "NA" or v.startswith("NA"):
                        # position non interprétable
                        continue

                    elif "MIX" in v:
                        # MIX considéré comme mutant
                        alleles.append(snp[-1])

                    else:
                        # mutant pur
                        alleles.append(snp[-1])

                haplo = "".join(alleles)
                sample_haplo[gene] = haplo if haplo else "Null"

            haplo_data.append(sample_haplo)

        haplo_df = pd.DataFrame(haplo_data)
        out_path = Path(output_dir) / "haplotypes_summary.csv"
        haplo_df.to_csv(out_path, index=False)

        log.write(f"Haplotype généré : {out_path}\n=== Fin ===\n")

    print(f"\n Haplotype généré avec succès : {out_path}")







def filter_haplotypes(input_file="output/haplotypes/haplotypes_summary.csv",
                      output_file="output/haplotypes/filtered_haplotypes_summary.csv"):
    """
    Charge haplotypes_summary.csv, supprime CytB & K13,
    identifie automatiquement les colonnes haplotype,
    supprime les lignes où toutes les valeurs haplotypes commencent par 'Nul*',
    puis sauvegarde un fichier nettoyé.
    """

    # Charger CSV
    df = pd.read_csv(input_file)

    # Colonnes à supprimer
    df = df.drop(columns=[c for c in ["CytB", "K13"] if c in df.columns], errors="ignore")

    # Motif null-like très robuste
    null_like = re.compile(r'^\s*(nul|null)[\w\(\)\-\/]*\s*$', flags=re.IGNORECASE)

    # Déterminer les colonnes haplotype (toutes sauf celles d'identification)
    non_haplo_cols = {"Sample", "SampleID", "Sample_name", "Sample_ID",
                      "Run", "Barcode", "ReadDepth", "Depth"}
    haplo_cols = [c for c in df.columns if c not in non_haplo_cols]

    # Fonction : vérifie si TOUTES les colonnes haplotype d’une ligne commencent par "Nul"
    def row_all_null_like(row):
        for col in haplo_cols:
            val = str(row[col]).strip()
            if not null_like.match(val):
                return False
        return True

    # Filtrage
    mask = ~df.apply(row_all_null_like, axis=1)
    df_filtered = df[mask].copy()

    # Sauvegarde
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    df_filtered.to_csv(output_file, index=False)

    print("✔ Colonnes haplotype utilisées :", haplo_cols)
    print(f"✔ Fichier filtré enregistré dans : {output_file}")





def run_combined_haplotypes(input_file="output/haplotypes/filtered_haplotypes_summary.csv", out_dir="output/haplotypes",
                                        output_name="Combined_Haplotypes.csv"):

    # Création du dossier de sortie
    os.makedirs(out_dir, exist_ok=True)

    # Lecture du fichier
    df = pd.read_csv(input_file)

    # Construction des nouvelles colonnes
    df["DHFR and DHPS"] = df["DHFR"].astype(str) + "/" + df["DHPS"].astype(str)
    df["CRT"] = df["CRT"].astype(str).str[:5]
    df["MDR"] = df["MDR"].astype(str).str[:3]

    # Sélection des colonnes finales
    df_out = df[["Sample", "DHFR and DHPS", "CRT", "MDR"]]

    # Enregistrement CSV
    output_path = os.path.join(out_dir, output_name)
    df_out.to_csv(output_path, index=False)

    print(f"✅ Fichier généré : {output_path}")







#======================================================================
# Generation du rapport final par site
#======================================================================



def generate_final_report_by_site_1(
    reportable_file="output/Dataviz_Reportable_snps/Reportable_snps_DMS_EPI_report.csv",
    combined_hap_file="output/haplotypes/Combined_Haplotypes.csv",
    vaf_file="output/SVAF_merge/SVAF_merge.csv",
    out_dir="output/haplotypes/Report",
    pdf_name="Final_Report.pdf"
):
    from pathlib import Path
    import pandas as pd
    import re
    import os
    import matplotlib.pyplot as plt
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4 
    from reportlab.lib.styles import getSampleStyleSheet
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, PageTemplate, Frame

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    img_dir = out_dir / "tmp_images"
    img_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = out_dir / pdf_name

    # --- Read inputs ---
    df = pd.read_csv(reportable_file, dtype=str).fillna("")
    dh = pd.read_csv(combined_hap_file, dtype=str).fillna("")
    df_vaf = pd.read_csv(vaf_file, dtype=str).fillna("")

    # --- Add SITE columns ---
    df["SITE"] = df["LSDB_Sequence_ID"].str[4:6]
    dh["SITE"] = dh["Sample"].str[4:6]
    df_vaf["SITE"] = df_vaf["Sample_name"].str[4:6]
    sites = df["SITE"].unique()

    # --- helper to get sample weight ---
    def get_sample_weight(sample_name):
        m = re.search(r'P(\d+)', str(sample_name))
        if m:
            return int(m.group(1))
        return 1

    # --- helper to count MT occurrences ---
    def count_mt(col_name, data):
        ser = data[col_name].astype(str).str.strip().str.upper()
        mask = (ser == "MT")
        return int(data.loc[mask, "__WEIGHT__"].sum())

    # --- robust null-like tester ---
    null_like_re = re.compile(r'^\s*(Null.*|Nul.*|na|n/?a|nan|-+)?\s*$', flags=re.IGNORECASE)
    def is_null_like(series: pd.Series):
        s = series.fillna("").astype(str).str.strip()
        return s.str.match(null_like_re)

    # --- Prepare PDF ---
    styles = getSampleStyleSheet()
    story = []

    # --- CIGASS logo ---
    logo_path = out_dir / "logoCIGASS.png"
    if logo_path.exists():
        story.append(Image(str(logo_path), width=370, height=120))
        story.append(Spacer(1, 10))

    story.append(Paragraph("<b>Final Report SNP & Haplotype</b>", styles["Title"]))
    story.append(Spacer(1, 14))

    # --- Loop over sites ---
    for site in sites:
        story.append(Paragraph(f"<b>Site : {site}</b>", styles["Heading1"]))
        story.append(Spacer(1, 10))

        # Filter data for the site
        site_df = df[df["SITE"] == site].copy()
        site_dh = dh[dh["SITE"] == site].copy()
        site_vaf = df_vaf[df_vaf["SITE"] == site].copy()

        # Sample weights
        site_df["__WEIGHT__"] = site_df["LSDB_Sequence_ID"].apply(get_sample_weight)
        total_samples = site_df["__WEIGHT__"].sum()
        story.append(Paragraph(f"Total samples : <b>{total_samples}</b>", styles["Normal"]))
        story.append(Spacer(1, 10))

        # --- Build VAF dictionary ---
        vaf_pop_dict = {}
        for _, r in site_vaf.iterrows():
            snp = r["Gene_ANNOTATION"] if "Gene_ANNOTATION" in r else r.get("AA_change", "")
            vaf = r.get("%VAF_FINAL_SNP", "")
            if snp:
                vaf_pop_dict[snp] = vaf

        # --- Identify gene blocks ---
        gene_blocks = {}
        current_gene = None
        for col in site_df.columns:
            if ": # Drug resistant mutations" in col:
                current_gene = col.split(":")[0].strip()
                gene_blocks[current_gene] = []
            else:
                if current_gene is not None:
                    gene_blocks[current_gene].append(col)
        for g in list(gene_blocks.keys()):
            cols = [c for c in gene_blocks[g] if re.search(r"\d", c)]
            gene_blocks[g] = cols
            if len(cols) == 0:
                gene_blocks.pop(g, None)
        gene_blocks.pop("CytoB", None)

        # --- Define order ---
        order = []
        if "DHFR" in gene_blocks or "DHPS" in gene_blocks:
            order.append(("DHFR and DHPS", ["DHFR","DHPS"]))
        for g in ["CRT","MDR"]:
            if g in gene_blocks:
                order.append((g, [g]))
        for g in gene_blocks:
            if g not in ("DHFR","DHPS","CRT","MDR","CytB"):
                order.append((g,[g]))

        # --- Loop over gene sections ---
        for section_name, genes_in_section in order:
            story.append(Paragraph(f"<b>{section_name}</b>", styles["Heading2"]))
            story.append(Spacer(1, 8))

            for gene in genes_in_section:
                if gene not in gene_blocks:
                    continue
                snp_cols = gene_blocks[gene]
                if len(snp_cols) == 0:
                    continue

                story.append(Paragraph(f"<b>{gene} — SNP</b>", styles["Heading3"]))
                story.append(Spacer(1, 6))

                table_data = [["SNP", "N of samples", "Percent (%)", "%VAF_pop"]]
                counts = {}
                for col in snp_cols:
                    c = count_mt(col, site_df)
                    pct = (c / total_samples * 100) if total_samples > 0 else 0.0
                    counts[col] = c
                    vaf_pop = "NA" if c == 0 else vaf_pop_dict.get(col, "NA")
                    table_data.append([col, str(c), f"{pct:.1f}%", vaf_pop])

                table = Table(table_data, colWidths=[140, 80, 80, 80])
                table.setStyle(TableStyle([
                    ("BACKGROUND", (0,0), (-1,0), colors.lightblue),
                    ("GRID", (0,0), (-1,-1), 0.4, colors.grey),
                    ("ALIGN", (1,1), (-1,-1), "CENTER")
                ]))
                story.append(table)
                story.append(Spacer(1, 8))

                if sum(counts.values()) > 0:
                    img_path = img_dir / f"{site}_{gene}_snp_bar.png"
                    plt.figure(figsize=(6,3.5))
                    plt.bar(list(counts.keys()), list(counts.values()))
                    plt.xticks(rotation=45, ha="right")
                    plt.ylabel("N of samples")
                    plt.title(f"{gene} — SNP Mutation Distribution")
                    plt.tight_layout()
                    plt.savefig(img_path, dpi=150)
                    plt.close()
                    story.append(Image(str(img_path), width=440, height=220))
                    story.append(Spacer(1, 12))
                else:
                    story.append(Paragraph("Aucune mutation détectée pour ce gène.", styles["Normal"]))
                    story.append(Spacer(1, 8))

            # --- Haplotypes ---
            for gene in genes_in_section:
                if gene not in site_dh.columns:
                    continue
                series = site_dh[gene].astype(str).str.strip()
                valid_mask = ~series.str.contains(r'Nul', case=False)
                df_h_valid = site_dh.loc[valid_mask].copy()
                if df_h_valid.empty:
                    story.append(Paragraph(f"Aucun haplotype valide trouvé pour {gene}.", styles["Normal"]))
                    story.append(Spacer(1,8))
                    continue

                df_h_valid["__WEIGHT__"] = df_h_valid["Sample"].apply(get_sample_weight)
                freq = df_h_valid.groupby(gene)["__WEIGHT__"].sum().reset_index()
                freq.columns = ["Haplotype", "Count"]
                total_h = freq["Count"].sum()
                freq["Percent"] = (freq["Count"] / total_h * 100).round(1)

                story.append(Paragraph(f"<b>Haplotypes {gene}</b>", styles["Heading3"]))
                story.append(Spacer(1,6))
                story.append(Paragraph(f"<b>Total haplotypes detected ({gene}) : {total_h}</b>", styles["Normal"]))
                story.append(Spacer(1, 6))

                table_data = [["Haplotype", "N of samples", "Percent"]]
                for _, r in freq.iterrows():
                    table_data.append([r["Haplotype"], int(r["Count"]), f"{r['Percent']}%"])

                t = Table(table_data, colWidths=[220, 80, 80])
                t.setStyle(TableStyle([
                    ("BACKGROUND", (0,0), (-1,0), colors.lightgrey),
                    ("GRID", (0,0), (-1,-1), 0.4, colors.grey),
                    ("ALIGN", (1,1), (-1,-1), "CENTER")
                ]))
                story.append(t)
                story.append(Spacer(1,8))

                img_path = img_dir / f"{site}_{gene}_haplo.png"
                plt.figure(figsize=(6,3.5))
                plt.bar(freq["Haplotype"].astype(str), freq["Count"])
                plt.xticks(rotation=45, ha="right")
                plt.title(f"{gene} haplotype distribution")
                plt.tight_layout()
                plt.savefig(img_path, dpi=150)
                plt.close()
                story.append(Image(str(img_path), width=440, height=220))
                story.append(Spacer(1, 12))
    
        story.append(Spacer(1, 18))

    # --- Build PDF ---
    from reportlab.lib.units import mm
    from reportlab.platypus import PageTemplate, Frame

    def draw_header(canvas, doc, logo_path):
        canvas.saveState()
        canvas.setFillColor(colors.lightblue)
        canvas.rect(0, doc.height + doc.topMargin + 10, doc.width + 2*doc.leftMargin, 25, fill=1)
        if os.path.exists(logo_path):
            canvas.drawImage(logo_path, x=doc.leftMargin, y=doc.height + doc.topMargin + 12, width=55, height=18, preserveAspectRatio=True, mask='auto')
        canvas.setFillColor(colors.black)
        canvas.setFont("Helvetica-Bold", 12)
        canvas.drawString(doc.leftMargin + 55, doc.height + doc.topMargin + 18, "CIGASS — UCAD — Sénégal")
        canvas.restoreState()

    def draw_footer(canvas, doc, logo_path):
        canvas.saveState()
        footer_y = 8 * mm
        if os.path.exists(logo_path):
            canvas.drawImage(logo_path, x=doc.leftMargin, y=footer_y, width=40, height=15, preserveAspectRatio=True, mask='auto')
        canvas.setFont("Helvetica", 10)
        canvas.drawRightString(doc.width + doc.leftMargin, footer_y + 3, f"Page {doc.page}")
        canvas.restoreState()

    def decorate_page(canvas, doc):
        draw_header(canvas, doc, str(logo_path))
        draw_footer(canvas, doc, str(logo_path))

    doc = SimpleDocTemplate(str(pdf_path), pagesize=A4)
    frame = Frame(doc.leftMargin, doc.bottomMargin+30, doc.width, doc.height-60, id="normal")
    page_template = PageTemplate(id="decorated", frames=[frame], onPage=decorate_page)
    doc.addPageTemplates([page_template])
    doc.build(story)

    # cleanup images
    for f in img_dir.glob("*"):
        try: f.unlink()
        except: pass

    print(f"✅ Rapport combiné généré : {pdf_path}")
    return str(pdf_path)





from pathlib import Path
import pandas as pd
import re
import os
import matplotlib.pyplot as plt
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, PageTemplate, Frame
from reportlab.lib.units import mm
from reportlab.lib.pagesizes import landscape, A4
from reportlab.platypus import Paragraph
from reportlab.lib.styles import ParagraphStyle


def generate_final_report_by_site(
    reportable_file="output/Dataviz_Reportable_snps/Reportable_snps_DMS_EPI_report.csv",
    combined_hap_file="output/haplotypes/Combined_Haplotypes.csv",
    vaf_file="output/SVAF_merge/SVAF_merge.csv",
    out_dir="output/haplotypes/Report",
    pdf_name="Final_Report.pdf"
):

    # --- Créer les dossiers ---
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    img_dir = out_dir / "tmp_images"
    img_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = out_dir / pdf_name

    # --- Lire les fichiers ---
    df = pd.read_csv(reportable_file, dtype=str).fillna("")
    dh = pd.read_csv(combined_hap_file, dtype=str).fillna("")
    df_vaf = pd.read_csv(vaf_file, dtype=str).fillna("")

    # --- Extraire le site depuis le nom d'échantillon ---
    df["SITE"] = df["LSDB_Sequence_ID"].str[4:6]
    dh["SITE"] = dh["Sample"].str[4:6]
    df_vaf["SITE"] = df_vaf["Sample_name"].str[4:6]
    sites = df["SITE"].unique()

    # --- Fonctions utilitaires ---
    def get_sample_weight(sample_name):
        m = re.search(r'P(\d+)', str(sample_name))
        return int(m.group(1)) if m else 1

    def load_svaf(sample_name, snp, sample_vaf_dir="output/Sample_VAF"):
        pattern = f"{sample_name}*_SVAF.csv"
        files = list(Path(sample_vaf_dir).glob(pattern))
        if not files:
            return 0.0
        try:
            df_svaf = pd.read_csv(files[0])
        except Exception:
            return 0.0
        if "AA_change" not in df_svaf.columns or "SVAF" not in df_svaf.columns:
            return 0.0
        row = df_svaf.loc[df_svaf["AA_change"] == snp]
        if row.empty:
            return 0.0
        try:
            return float(row.iloc[0]["SVAF"])
        except Exception:
            return 0.0

    # --- Préparer le PDF ---
    styles = getSampleStyleSheet()
    story = []

    # Logo
    logo_path = Path(out_dir).parent.parent.parent / "images" / "logoCIGASS.png"
    if logo_path.exists():
        story.append(Image(str(logo_path), width=370, height=120))
        story.append(Spacer(1, 10))

    story.append(Paragraph("<b>Final Report SNP & Haplotype</b>", styles["Title"]))
    story.append(Spacer(1, 14))


    # ======================================================
    # Résumé global au début du rapport (corrigé pour J0/JE)
    # ======================================================

    def get_patient_id(sample_id):
        """
        Retourne l'ID principal du patient en ignorant le jour (day)
        Exemple : "23SNKD00I0009PfB4721" et "23SNKD28I0009PfB4721" -> même patient_id
        """
        sid = str(sample_id).split("_")[0]  # tout avant "_"
        # retirer le jour (positions 6-7)
        return sid[:6] + sid[8:]

    def get_day(sample_id):
        
        #Day extrait de l'ID principal (positions 6-7)
        
        sid = str(sample_id).split("_")[0]
        return sid[6:8] if len(sid) >= 8 else None

    metrics = [
    "Number of patients (J0/JE pairs)",
    "Number of individual samples",     # NON poolés uniquement
    "Day 0 samples",
    "Day Failure samples",              # NOUVELLE LIGNE
    "Number of Pools",
    "Number of pooled Day 0 samples"
    ]

    summary_data = [["Metric"] + list(sites) + ["Total"]]

    for metric in metrics:
        row = [metric]
        total_value = 0

        for site in sites:
            site_df = df[df["SITE"] == site].copy()
            site_df["__WEIGHT__"] = site_df["LSDB_Sequence_ID"].apply(get_sample_weight)
            site_df["PATIENT_ID"] = site_df["LSDB_Sequence_ID"].apply(get_patient_id)
            site_df["DAY"] = site_df["LSDB_Sequence_ID"].apply(get_day)

            if metric == "Number of patients (J0/JE pairs)":
                n_patients = 0
                for pid, g in site_df.groupby("PATIENT_ID"):
                    has_j0 = (g["DAY"] == "00").any()
                    has_je = (g["DAY"] != "00").any()
                    if has_j0 and has_je:
                        n_patients += 1
                value = n_patients

            elif metric == "Number of individual samples":
                value = site_df[site_df["__WEIGHT__"] == 1].shape[0]

            elif metric == "Day 0 samples":
                value = site_df[
                    (site_df["DAY"] == "00") &
                    (site_df["__WEIGHT__"] == 1)
                ].shape[0]

            elif metric == "Day Failure samples":
                value = site_df[
                    (site_df["DAY"] != "00") &
                    (site_df["__WEIGHT__"] == 1)
                ].shape[0]

            elif metric == "Number of Pools":
                value = (site_df["__WEIGHT__"] > 1).sum()

            elif metric == "Number of pooled Day 0 samples":
                mask = (site_df["DAY"] == "00") & (site_df["__WEIGHT__"] > 1)
                value = site_df.loc[mask, "__WEIGHT__"].sum()
                
            else:
                value = 0

            row.append(int(value))
            total_value += value

        row.append(int(total_value))
        summary_data.append(row)

    
    # ============================
    # Affichage dans le PDF
    # ============================

    story.append(Paragraph("<b>Global Summary by Site</b>", styles["Heading2"]))
    story.append(Spacer(1, 6))
    
    from reportlab.lib.pagesizes import A4

    # --- Calcul dynamique des largeurs ---
    page_width = A4[0] - 2*30  # largeur A4 moins marges
    n_cols = len(sites) + 2    # "Metric" + sites + "Total"

    metric_col_width = page_width * 0.25            # 25% pour la colonne "Metric"
    other_col_width = (page_width - metric_col_width) / (n_cols - 1)
    col_widths = [metric_col_width] + [other_col_width]*(n_cols - 1)

    # --- Création du tableau ---
    table_summary = Table(summary_data, colWidths=col_widths)

    # --- Style du tableau ---
    table_style = TableStyle([
        ("BACKGROUND", (0,0), (-1,0), colors.lightblue),
        ("GRID", (0,0), (-1,-1), 0.4, colors.grey),
        ("ALIGN", (1,1), (-1,-1), "CENTER"),
        ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
        ("FONTSIZE", (0,0), (-1,-1), 8),      # taille police réduite
        ("ROTATE", (1,0), (-2,0), 90),        # rotation des en-têtes sites
        ("ALIGN", (0,0), (0,-1), "LEFT")      # colonne Metric alignée à gauche
    ])

    table_summary.setStyle(table_style)

    # --- Ajout au PDF ---
    story.append(Paragraph("<b>Global Summary by Site</b>", styles["Heading2"]))
    story.append(Spacer(1, 6))
    story.append(table_summary)
    story.append(Spacer(1,12))

        

    
    # --- Boucle sur les sites ---
    for site in sites:
        story.append(Paragraph(f"<b>Site : {site}</b>", styles["Heading1"]))
        story.append(Spacer(1, 10))

        site_df = df[df["SITE"] == site].copy()
        site_dh = dh[dh["SITE"] == site].copy()
        site_vaf = df_vaf[df_vaf["SITE"] == site].copy()

        site_df["__WEIGHT__"] = site_df["LSDB_Sequence_ID"].apply(get_sample_weight)
        total_samples = site_df["__WEIGHT__"].sum()
        story.append(Paragraph(f"Total samples : <b>{total_samples}</b>", styles["Normal"]))
        story.append(Spacer(1, 10))

        # --- Identifier les blocs de gènes ---
        gene_blocks = {}
        current_gene = None
        for col in site_df.columns:
            if ": # Drug resistant mutations" in col:
                current_gene = col.split(":")[0].strip()
                gene_blocks[current_gene] = []
            elif current_gene is not None:
                gene_blocks[current_gene].append(col)
        for g in list(gene_blocks.keys()):
            cols = [c for c in gene_blocks[g] if re.search(r"\d", c)]
            gene_blocks[g] = cols
            if len(cols) == 0:
                gene_blocks.pop(g, None)
        gene_blocks.pop("CytoB", None)

        # --- Définir l'ordre des sections ---
        order = []
        if "DHFR" in gene_blocks or "DHPS" in gene_blocks:
            order.append(("DHFR and DHPS", ["DHFR","DHPS"]))
        for g in ["CRT","MDR"]:
            if g in gene_blocks:
                order.append((g, [g]))
        for g in gene_blocks:
            if g not in ("DHFR","DHPS","CRT","MDR","CytoB"):
                order.append((g,[g]))

        # --- Boucle sur les sections de gènes ---
        for section_name, genes_in_section in order:
            story.append(Paragraph(f"<b>{section_name}</b>", styles["Heading2"]))
            story.append(Spacer(1, 8))

            processed_haplo_sections = set()  # pour ne pas répéter DHFR/DHPS

            for gene in genes_in_section:
                snps = gene_blocks.get(gene, [])
                if not snps:
                    continue
                
                    
                # --- Tableau SNPs pour le gène ---


                # --- Titre du gène avant le tableau SNP ---
                story.append(Spacer(1, 6))
                story.append(Paragraph(f"<b>{gene}</b>", styles["Heading3"]))
                story.append(Spacer(1, 4))

                # --- Styles pour le tableau ---
                style_cell = ParagraphStyle(name='cell', alignment=1, fontSize=8)  # centré, petite police
                style_header = ParagraphStyle(name='header', alignment=1, fontSize=9, leading=10)

                # --- Préparer les données avec wrapping ---
                table_data_wrapped = []
                # Ajouter l'entête avec style_header
                table_data_wrapped.append([Paragraph(str(col), style_header) for col in ["SNP", "N of sample", "Wild-type (N, %)", 
                                                                                        "Mutant (N, %)", "%VAF_Mutant", 
                                                                                        "Mix (N, %)", "%VAF_Mix"]])

                snp_values_for_plot = []

                for snp in snps:
                    df_snp = site_df[[snp, "LSDB_Sequence_ID", "__WEIGHT__"]].copy()
                    df_snp["SVAF"] = df_snp["LSDB_Sequence_ID"].apply(lambda x: load_svaf(x, snp))

                    # Masques
                    mask_wt  = df_snp[snp] == "WT"
                    mask_mut = df_snp[snp] == "MT"
                    mask_mix = df_snp[snp] == "MIX"
                    mask_valid = df_snp[snp].isin(["WT", "MT", "MIX"])

                    # Total pondéré
                    n_total = df_snp.loc[mask_valid, "__WEIGHT__"].sum()

                    # Comptage pondéré
                    n_wt  = df_snp.loc[mask_wt, "__WEIGHT__"].sum()
                    n_mut = df_snp.loc[mask_mut, "__WEIGHT__"].sum()
                    n_mix = df_snp.loc[mask_mix, "__WEIGHT__"].sum()

                    # Pourcentages
                    pct_wt  = (n_wt / n_total * 100) if n_total > 0 else 0
                    pct_mut = (n_mut / n_total * 100) if n_total > 0 else 0
                    pct_mix = (n_mix / n_total * 100) if n_total > 0 else 0

                    # Numérateurs VAF
                    num_mut = (df_snp.loc[mask_mut, "SVAF"] * df_snp.loc[mask_mut, "__WEIGHT__"]).sum()
                    num_mix = (df_snp.loc[mask_mix, "SVAF"] * df_snp.loc[mask_mix, "__WEIGHT__"]).sum()

                    # Dénominateurs VAF
                    den_mut = df_snp.loc[mask_mut | mask_wt, "__WEIGHT__"].sum()
                    den_mix = df_snp.loc[mask_mix | mask_wt, "__WEIGHT__"].sum()

                    # %VAF
                    vaf_mut = (num_mut / den_mut * 100) if den_mut > 0 else 0
                    vaf_mix = (num_mix / den_mix * 100) if den_mix > 0 else 0

                    # Ajouter la ligne avec style_cell
                    table_data_wrapped.append([
                        Paragraph(str(snp), style_cell),
                        Paragraph(str(int(n_total)), style_cell),
                        Paragraph(f"({int(n_wt)}, {pct_wt:.1f}%)", style_cell),
                        Paragraph(f"({int(n_mut)}, {pct_mut:.1f}%)", style_cell),
                        Paragraph(f"{vaf_mut:.1f}%", style_cell),
                        Paragraph(f"({int(n_mix)}, {pct_mix:.1f}%)", style_cell),
                        Paragraph(f"{vaf_mix:.1f}%", style_cell)
                    ])

                    # Pour le graphe
                    #snp_values_for_plot.append((snp, int(n_mut), int(n_mix)))
                    snp_values_for_plot.append((snp, int(n_wt), int(n_mut), int(n_mix)))

                # --- Créer le tableau ---
                col_widths = [50,50,70,70,50,70,50]  # ajusté pour paysage A4
                table = Table(table_data_wrapped, colWidths=col_widths, repeatRows=1)
                table.setStyle(TableStyle([
                    ("BACKGROUND",(0,0),(-1,0),colors.lightblue),
                    ("GRID",(0,0),(-1,-1),0.4,colors.grey),
                    ("ALIGN",(0,0),(-1,-1),"CENTER"),
                ]))

                story.append(table)
                story.append(Spacer(1,6))

                # --- Graphe SNP avec WT / Mutant / Mix ---
                if snp_values_for_plot:
                    img_path = img_dir / f"{site}_{gene}_SNP.png"
                    plt.figure(figsize=(6.5,3.8))

                    labels = [x[0] for x in snp_values_for_plot]
                    wt_counts  = [x[1] for x in snp_values_for_plot]
                    mut_counts = [x[2] for x in snp_values_for_plot]
                    mix_counts = [x[3] for x in snp_values_for_plot]

                    width = 0.25
                    x = range(len(labels))

                    plt.bar(x, wt_counts, width, label="Wild-type")
                    plt.bar([i + width for i in x], mut_counts, width, label="Mutant")
                    plt.bar([i + 2*width for i in x], mix_counts, width, label="Mix")

                    plt.xticks([i + width for i in x], labels, rotation=45, ha="right")
                    plt.ylabel("Nombre d'échantillons")
                    plt.title(f"{gene} SNP distribution — Site {site}")
                    plt.legend()
                    plt.tight_layout()

                    plt.savefig(img_path, dpi=150)
                    plt.close()

                    story.append(Image(str(img_path), width=440, height=220))
                    story.append(Spacer(1,12))


                    
                    
            # --- Haplotypes pour le bloc de gènes ---
            haplo_col = section_name  # "DHFR and DHPS" ou "CRT", etc
            if haplo_col in processed_haplo_sections:
                continue
            processed_haplo_sections.add(haplo_col)

            if haplo_col not in site_dh.columns:
                story.append(Paragraph(f"Aucun haplotype valide trouvé pour {haplo_col}.", styles["Normal"]))
                story.append(Spacer(1,6))
                continue

            df_h_valid = site_dh.copy()
            df_h_valid = df_h_valid[~df_h_valid[haplo_col].astype(str).str.contains(r'Nul',case=False)]
            if df_h_valid.empty:
                story.append(Paragraph(f"Aucun haplotype valide pour {haplo_col} — site {site}.", styles["Normal"]))
                story.append(Spacer(1,6))
                continue

            df_h_valid["__WEIGHT__"] = df_h_valid["Sample"].apply(get_sample_weight)
            freq = df_h_valid.groupby(haplo_col)["__WEIGHT__"].sum().reset_index()
            freq.columns = ["Haplotype","Count"]
            total_h = freq["Count"].sum()
            freq["Percent"] = (freq["Count"]/total_h*100).round(1)

            story.append(Paragraph(f"<b>Haplotypes {haplo_col}</b>", styles["Heading3"]))
            story.append(Spacer(1,6))
            story.append(Paragraph(f"<b>Total haplotypes detected ({haplo_col}) : {total_h}</b>", styles["Normal"]))
            story.append(Spacer(1,6))

            table_data = [["Haplotype","N of samples","Percent (%)"]]
            for _,r in freq.iterrows():
                table_data.append([r["Haplotype"],int(r["Count"]),f"{r['Percent']}%"])
            t = Table(table_data,colWidths=[220,80,80])
            t.setStyle(TableStyle([
                ("BACKGROUND",(0,0),(-1,0),colors.lightgrey),
                ("GRID",(0,0),(-1,-1),0.4,colors.grey),
                ("ALIGN",(1,1),(-1,-1),"CENTER")
            ]))
            story.append(t)
            story.append(Spacer(1,6))

            # Graphe haplotype
            img_path = img_dir / f"{site}_{haplo_col}_haplo.png"
            plt.figure(figsize=(6,3.5))
            plt.bar(freq["Haplotype"].astype(str),freq["Count"])
            plt.xticks(rotation=45,ha="right")
            plt.title(f"{haplo_col} haplotype distribution — Site {site}")
            plt.tight_layout()
            plt.savefig(img_path,dpi=150)
            plt.close()
            story.append(Image(str(img_path),width=440,height=220))
            story.append(Spacer(1,12))


        story.append(Spacer(1,18))

    # --- Header & Footer ---
    def draw_header(canvas, doc, logo_path):
        canvas.saveState()
        canvas.setFillColor(colors.lightblue)
        canvas.rect(0, doc.height+doc.topMargin+10, doc.width+2*doc.leftMargin, 25, fill=1)
        if os.path.exists(logo_path):
            canvas.drawImage(logo_path, x=doc.leftMargin, y=doc.height+doc.topMargin+12, width=55, height=18, preserveAspectRatio=True, mask='auto')
        canvas.setFillColor(colors.black)
        canvas.setFont("Helvetica-Bold",12)
        canvas.drawString(doc.leftMargin+55, doc.height+doc.topMargin+18,"CIGASS — UCAD — Sénégal")
        canvas.restoreState()

    def draw_footer(canvas, doc, logo_path):
        canvas.saveState()
        footer_y = 8*mm
        if os.path.exists(logo_path):
            canvas.drawImage(logo_path, x=doc.leftMargin, y=footer_y, width=40, height=15, preserveAspectRatio=True, mask='auto')
        canvas.setFont("Helvetica",10)
        canvas.drawRightString(doc.width+doc.leftMargin, footer_y+3,f"Page {doc.page}")
        canvas.restoreState()

    def decorate_page(canvas, doc):
        draw_header(canvas, doc, str(logo_path))
        draw_footer(canvas, doc, str(logo_path))

    # --- Générer le PDF ---
    doc = SimpleDocTemplate(str(pdf_path),pagesize=A4)
    frame = Frame(doc.leftMargin, doc.bottomMargin+30, doc.width, doc.height-60, id="normal")
    page_template = PageTemplate(id="decorated", frames=[frame], onPage=decorate_page)
    doc.addPageTemplates([page_template])
    doc.build(story)

    # --- Nettoyer les images temporaires ---
    for f in img_dir.glob("*"):
        try: f.unlink()
        except: pass

    print(f"✅ Rapport combiné généré : {pdf_path}")
    return str(pdf_path)




from pathlib import Path
import pandas as pd
import re
import os
import matplotlib.pyplot as plt
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, PageTemplate, Frame
from reportlab.lib.units import mm
from reportlab.platypus import Paragraph
from reportlab.lib.styles import ParagraphStyle

def generate_final_report_by_site_0(
    reportable_file="output/Dataviz_Reportable_snps/Reportable_snps_DMS_EPI_report.csv",
    combined_hap_file="output/haplotypes/Combined_Haplotypes.csv",
    vaf_file="output/SVAF_merge/SVAF_merge.csv",
    out_dir="output/haplotypes/Report",
    pdf_name="Final_Report_MT_MIX.pdf"
):

    # --- Créer les dossiers ---
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    img_dir = out_dir / "tmp_images"
    img_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = out_dir / pdf_name

    # --- Lire les fichiers ---
    df = pd.read_csv(reportable_file, dtype=str).fillna("")
    dh = pd.read_csv(combined_hap_file, dtype=str).fillna("")
    df_vaf = pd.read_csv(vaf_file, dtype=str).fillna("")

    # --- Extraire le site depuis le nom d'échantillon ---
    df["SITE"] = df["LSDB_Sequence_ID"].str[4:6]
    dh["SITE"] = dh["Sample"].str[4:6]
    df_vaf["SITE"] = df_vaf["Sample_name"].str[4:6]
    sites = df["SITE"].unique()

    # --- Fonctions utilitaires ---
    def get_sample_weight(sample_name):
        m = re.search(r'P(\d+)', str(sample_name))
        return int(m.group(1)) if m else 1

    def load_svaf(sample_name, snp, sample_vaf_dir="output/Sample_VAF"):
        pattern = f"{sample_name}*_SVAF.csv"
        files = list(Path(sample_vaf_dir).glob(pattern))
        if not files:
            return 0.0
        try:
            df_svaf = pd.read_csv(files[0])
        except Exception:
            return 0.0
        if "AA_change" not in df_svaf.columns or "SVAF" not in df_svaf.columns:
            return 0.0
        row = df_svaf.loc[df_svaf["AA_change"] == snp]
        if row.empty:
            return 0.0
        try:
            return float(row.iloc[0]["SVAF"])
        except Exception:
            return 0.0

    # --- Préparer le PDF ---
    styles = getSampleStyleSheet()
    story = []

    # Logo
    logo_path = Path(out_dir).parent.parent.parent / "images" / "logoCIGASS.png"
    if logo_path.exists():
        story.append(Image(str(logo_path), width=370, height=120))
        story.append(Spacer(1, 10))

    story.append(Paragraph("<b>Final Report SNP & Haplotype (MIX as Mutant)</b>", styles["Title"]))
    story.append(Spacer(1, 14))

    # ======================================================
    # Résumé global au début du rapport (J0/JE)
    # ======================================================

    def get_patient_id(sample_id):
        sid = str(sample_id).split("_")[0]
        return sid[:6] + sid[8:]  # retirer le day pour ID patient

    def get_day(sample_id):
        sid = str(sample_id).split("_")[0]
        return sid[6:8] if len(sid) >= 8 else None
    
    metrics = [
    "Number of patients (J0/JE pairs)",
    "Number of individual samples",     # NON poolés uniquement
    "Day 0 samples",
    "Day Failure samples",              # NOUVELLE LIGNE
    "Number of Pools",
    "Number of pooled Day 0 samples"
    ]


    summary_data = [["Metric"] + list(sites) + ["Total"]]

    for metric in metrics:
        row = [metric]
        total_value = 0

        for site in sites:
            site_df = df[df["SITE"] == site].copy()
            site_df["__WEIGHT__"] = site_df["LSDB_Sequence_ID"].apply(get_sample_weight)
            site_df["PATIENT_ID"] = site_df["LSDB_Sequence_ID"].apply(get_patient_id)
            site_df["DAY"] = site_df["LSDB_Sequence_ID"].apply(get_day)

            if metric == "Number of patients (J0/JE pairs)":
                n_patients = 0
                for pid, g in site_df.groupby("PATIENT_ID"):
                    has_j0 = (g["DAY"] == "00").any()
                    has_je = (g["DAY"] != "00").any()
                    if has_j0 and has_je:
                        n_patients += 1
                value = n_patients

            elif metric == "Number of individual samples":
                value = site_df[site_df["__WEIGHT__"] == 1].shape[0]

            elif metric == "Day 0 samples":
                value = site_df[
                    (site_df["DAY"] == "00") &
                    (site_df["__WEIGHT__"] == 1)
                ].shape[0]

            elif metric == "Day Failure samples":
                value = site_df[
                    (site_df["DAY"] != "00") &
                    (site_df["__WEIGHT__"] == 1)
                ].shape[0]

            elif metric == "Number of Pools":
                value = (site_df["__WEIGHT__"] > 1).sum()

            elif metric == "Number of pooled Day 0 samples":
                mask = (site_df["DAY"] == "00") & (site_df["__WEIGHT__"] > 1)
                value = site_df.loc[mask, "__WEIGHT__"].sum()

            else:
                value = 0

            row.append(int(value))
            total_value += value

        row.append(int(total_value))
        summary_data.append(row)

    
    # --- Affichage dans le PDF ---
    story.append(Paragraph("<b>Global Summary by Site</b>", styles["Heading2"]))
    story.append(Spacer(1, 6))
    """
    table_summary = Table(summary_data, colWidths=[260] + [70]*len(sites) + [70])
    table_summary.setStyle(TableStyle([
        ("BACKGROUND", (0,0), (-1,0), colors.lightblue),
        ("GRID",(0,0),(-1,-1),0.4,colors.grey),
        ("ALIGN",(1,1),(-1,-1),"CENTER"),
        ("VALIGN",(0,0),(-1,-1),"MIDDLE")
    ]))
    story.append(table_summary)
    story.append(Spacer(1,12))
    """
    
    from reportlab.lib.pagesizes import A4

    # --- Calcul dynamique des largeurs ---
    page_width = A4[0] - 2*30  # largeur A4 moins marges
    n_cols = len(sites) + 2    # "Metric" + sites + "Total"

    metric_col_width = page_width * 0.25            # 25% pour la colonne "Metric"
    other_col_width = (page_width - metric_col_width) / (n_cols - 1)
    col_widths = [metric_col_width] + [other_col_width]*(n_cols - 1)

    # --- Création du tableau ---
    table_summary = Table(summary_data, colWidths=col_widths)

    # --- Style du tableau ---
    table_style = TableStyle([
        ("BACKGROUND", (0,0), (-1,0), colors.lightblue),
        ("GRID", (0,0), (-1,-1), 0.4, colors.grey),
        ("ALIGN", (1,1), (-1,-1), "CENTER"),
        ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
        ("FONTSIZE", (0,0), (-1,-1), 8),      # taille police réduite
        ("ROTATE", (1,0), (-2,0), 90),        # rotation des en-têtes sites
        ("ALIGN", (0,0), (0,-1), "LEFT")      # colonne Metric alignée à gauche
    ])

    table_summary.setStyle(table_style)

    # --- Ajout au PDF ---
    story.append(Paragraph("<b>Global Summary by Site</b>", styles["Heading2"]))
    story.append(Spacer(1, 6))
    story.append(table_summary)
    story.append(Spacer(1,12))

    
    # --- Boucle sur les sites ---
    for site in sites:
        story.append(Paragraph(f"<b>Site : {site}</b>", styles["Heading1"]))
        story.append(Spacer(1, 10))

        site_df = df[df["SITE"] == site].copy()
        site_dh = dh[dh["SITE"] == site].copy()
        site_vaf = df_vaf[df_vaf["SITE"] == site].copy()

        site_df["__WEIGHT__"] = site_df["LSDB_Sequence_ID"].apply(get_sample_weight)
        total_samples = site_df["__WEIGHT__"].sum()
        story.append(Paragraph(f"Total samples : <b>{total_samples}</b>", styles["Normal"]))
        story.append(Spacer(1, 10))

        # --- Identifier les blocs de gènes ---
        gene_blocks = {}
        current_gene = None
        for col in site_df.columns:
            if ": # Drug resistant mutations" in col:
                current_gene = col.split(":")[0].strip()
                gene_blocks[current_gene] = []
            elif current_gene is not None:
                gene_blocks[current_gene].append(col)
        for g in list(gene_blocks.keys()):
            cols = [c for c in gene_blocks[g] if re.search(r"\d", c)]
            gene_blocks[g] = cols
            if len(cols) == 0:
                gene_blocks.pop(g, None)
        gene_blocks.pop("CytoB", None)

        # --- Définir l'ordre des sections ---
        order = []
        if "DHFR" in gene_blocks or "DHPS" in gene_blocks:
            order.append(("DHFR and DHPS", ["DHFR","DHPS"]))
        for g in ["CRT","MDR"]:
            if g in gene_blocks:
                order.append((g, [g]))
        for g in gene_blocks:
            if g not in ("DHFR","DHPS","CRT","MDR","CytoB"):
                order.append((g,[g]))

        # --- Boucle sur les sections de gènes ---
        for section_name, genes_in_section in order:
            story.append(Paragraph(f"<b>{section_name}</b>", styles["Heading2"]))
            story.append(Spacer(1, 8))

            processed_haplo_sections = set()

            for gene in genes_in_section:
                snps = gene_blocks.get(gene, [])
                if not snps: continue

                # --- Titre du gène avant le tableau SNP ---
                story.append(Spacer(1, 6))
                story.append(Paragraph(f"<b>{gene}</b>", styles["Heading3"]))
                story.append(Spacer(1, 4))

                style_cell = ParagraphStyle(name='cell', alignment=1, fontSize=8)
                style_header = ParagraphStyle(name='header', alignment=1, fontSize=9, leading=10)

                table_data_wrapped = []
                table_data_wrapped.append([Paragraph(str(col), style_header) for col in ["SNP","N of sample","Wild-type (N, %)","Mutant (N, %)","%VAF_Mutant"]])

                snp_values_for_plot = []

                for snp in snps:
                    df_snp = site_df[[snp, "LSDB_Sequence_ID", "__WEIGHT__"]].copy()
                    df_snp["SVAF"] = df_snp["LSDB_Sequence_ID"].apply(lambda x: load_svaf(x, snp))

                    mask_wt  = df_snp[snp] == "WT"
                    mask_mut = df_snp[snp].isin(["MT","MIX"])  # <-- MIX compté comme MT
                    mask_valid = df_snp[snp].isin(["WT","MT","MIX"])

                    n_total = df_snp.loc[mask_valid, "__WEIGHT__"].sum()
                    n_wt = df_snp.loc[mask_wt, "__WEIGHT__"].sum()
                    n_mut = df_snp.loc[mask_mut, "__WEIGHT__"].sum()

                    pct_wt = (n_wt / n_total * 100) if n_total > 0 else 0
                    pct_mut = (n_mut / n_total * 100) if n_total > 0 else 0

                    num_mut = (df_snp.loc[mask_mut, "SVAF"] * df_snp.loc[mask_mut, "__WEIGHT__"]).sum()
                    den_mut = df_snp.loc[mask_mut | mask_wt, "__WEIGHT__"].sum()
                    vaf_mut = (num_mut / den_mut * 100) if den_mut > 0 else 0

                    table_data_wrapped.append([
                        Paragraph(str(snp), style_cell),
                        Paragraph(str(int(n_total)), style_cell),
                        Paragraph(f"({int(n_wt)}, {pct_wt:.1f}%)", style_cell),
                        Paragraph(f"({int(n_mut)}, {pct_mut:.1f}%)", style_cell),
                        Paragraph(f"{vaf_mut:.1f}%", style_cell)
                    ])

                    snp_values_for_plot.append((snp, int(n_wt), int(n_mut)))

                # --- Créer le tableau ---
                col_widths = [50,50,100,100,60]
                table = Table(table_data_wrapped, colWidths=col_widths, repeatRows=1)
                table.setStyle(TableStyle([
                    ("BACKGROUND",(0,0),(-1,0),colors.lightblue),
                    ("GRID",(0,0),(-1,-1),0.4,colors.grey),
                    ("ALIGN",(0,0),(-1,-1),"CENTER")
                ]))
                story.append(table)
                story.append(Spacer(1,6))

                # --- Graphe SNP WT/Mutant ---
                if snp_values_for_plot:
                    img_path = img_dir / f"{site}_{gene}_SNP.png"
                    plt.figure(figsize=(6.5,3.8))
                    labels = [x[0] for x in snp_values_for_plot]
                    wt_counts  = [x[1] for x in snp_values_for_plot]
                    mut_counts = [x[2] for x in snp_values_for_plot]

                    width = 0.35
                    x = range(len(labels))
                    plt.bar(x, wt_counts, width, label="Wild-type")
                    plt.bar([i + width for i in x], mut_counts, width, label="Mutant")

                    plt.xticks([i + width/2 for i in x], labels, rotation=45, ha="right")
                    plt.ylabel("Nombre d'échantillons")
                    plt.title(f"{gene} SNP distribution — Site {site}")
                    plt.legend()
                    plt.tight_layout()
                    plt.savefig(img_path, dpi=150)
                    plt.close()
                    story.append(Image(str(img_path), width=440, height=220))
                    story.append(Spacer(1,12))

            # --- Haplotypes (inchangé) ---
            haplo_col = section_name
            if haplo_col in processed_haplo_sections: continue
            processed_haplo_sections.add(haplo_col)

            if haplo_col not in site_dh.columns:
                story.append(Paragraph(f"Aucun haplotype valide trouvé pour {haplo_col}.", styles["Normal"]))
                story.append(Spacer(1,6))
                continue

            df_h_valid = site_dh.copy()
            df_h_valid = df_h_valid[~df_h_valid[haplo_col].astype(str).str.contains(r'Nul',case=False)]
            if df_h_valid.empty:
                story.append(Paragraph(f"Aucun haplotype valide pour {haplo_col} — site {site}.", styles["Normal"]))
                story.append(Spacer(1,6))
                continue

            df_h_valid["__WEIGHT__"] = df_h_valid["Sample"].apply(get_sample_weight)
            freq = df_h_valid.groupby(haplo_col)["__WEIGHT__"].sum().reset_index()
            freq.columns = ["Haplotype","Count"]
            total_h = freq["Count"].sum()
            freq["Percent"] = (freq["Count"]/total_h*100).round(1)

            story.append(Paragraph(f"<b>Haplotypes {haplo_col}</b>", styles["Heading3"]))
            story.append(Spacer(1,6))
            story.append(Paragraph(f"<b>Total haplotypes detected ({haplo_col}) : {total_h}</b>", styles["Normal"]))
            story.append(Spacer(1,6))

            table_data = [["Haplotype","N of samples","Percent (%)"]]
            for _,r in freq.iterrows():
                table_data.append([r["Haplotype"],int(r["Count"]),f"{r['Percent']}%"])
            t = Table(table_data,colWidths=[220,80,80])
            t.setStyle(TableStyle([
                ("BACKGROUND",(0,0),(-1,0),colors.lightgrey),
                ("GRID",(0,0),(-1,-1),0.4,colors.grey),
                ("ALIGN",(1,1),(-1,-1),"CENTER")
            ]))
            story.append(t)
            story.append(Spacer(1,6))

            img_path = img_dir / f"{site}_{haplo_col}_haplo.png"
            plt.figure(figsize=(6,3.5))
            plt.bar(freq["Haplotype"].astype(str),freq["Count"])
            plt.xticks(rotation=45,ha="right")
            plt.title(f"{haplo_col} haplotype distribution — Site {site}")
            plt.tight_layout()
            plt.savefig(img_path,dpi=150)
            plt.close()
            story.append(Image(str(img_path),width=440,height=220))
            story.append(Spacer(1,12))

        story.append(Spacer(1,18))

    # --- Header & Footer (inchangé) ---
    def draw_header(canvas, doc, logo_path):
        canvas.saveState()
        canvas.setFillColor(colors.lightblue)
        canvas.rect(0, doc.height+doc.topMargin+10, doc.width+2*doc.leftMargin, 25, fill=1)
        if os.path.exists(logo_path):
            canvas.drawImage(logo_path, x=doc.leftMargin, y=doc.height+doc.topMargin+12, width=55, height=18, preserveAspectRatio=True, mask='auto')
        canvas.setFillColor(colors.black)
        canvas.setFont("Helvetica-Bold",12)
        canvas.drawString(doc.leftMargin+55, doc.height+doc.topMargin+18,"CIGASS — UCAD — Sénégal")
        canvas.restoreState()

    def draw_footer(canvas, doc, logo_path):
        canvas.saveState()
        footer_y = 8*mm
        if os.path.exists(logo_path):
            canvas.drawImage(logo_path, x=doc.leftMargin, y=footer_y, width=40, height=15, preserveAspectRatio=True, mask='auto')
        canvas.setFont("Helvetica",10)
        canvas.drawRightString(doc.width+doc.leftMargin, footer_y+3,f"Page {doc.page}")
        canvas.restoreState()

    def decorate_page(canvas, doc):
        draw_header(canvas, doc, str(logo_path))
        draw_footer(canvas, doc, str(logo_path))

    # --- Générer le PDF ---
    doc = SimpleDocTemplate(str(pdf_path),pagesize=A4)
    frame = Frame(doc.leftMargin, doc.bottomMargin+30, doc.width, doc.height-60, id="normal")
    page_template = PageTemplate(id="decorated", frames=[frame], onPage=decorate_page)
    doc.addPageTemplates([page_template])
    doc.build(story)

    # --- Nettoyer les images temporaires ---
    for f in img_dir.glob("*"):
        try: f.unlink()
        except: pass

    print(f"✅ Rapport combiné généré : {pdf_path}")
    return str(pdf_path)





if __name__ == "__main__":
    
    #======================================================================
    #        Pipeline Main
    #======================================================================

    
    print("[PIPELINE] Début des l'analyses.............")
    global_start = time.time()
    """
    
    # Lancement QC_pre_trimming
    QC_pre_trimming(input_dir="data", output_dir="output/QC_pre_trimming")


    # Lancement MultiQC_pre_trimming
    MultiQC_pre_trimming(output_dir="output/QC_pre_trimming", report_name="MultiQC_pre_trimming_report.html")
    

    # Lancement Trimming avec BBduk
    trimming(input_dir="data", output_dir="output/trimmed_reads", adapter_file="pf_3D7_Ref/adapters.fa")


    # Lancement QC_post_trimming
    QC_post_trimming(input_dir="output/trimmed_reads", output_dir="output/QC_post_trimming")


    # Lancement MultiQC_post_trimming
    MultiQC_post_trimming(output_dir="output/QC_post_trimming", report_name="MultiQC_post_trimming_report.html")

    
    # Lancement bwa_index
    reference = "pf_3D7_Ref/mars_pf_ref.fasta"  # chemin vers la référence
    bwa_index(reference)

    
    # Lancement bwa_align
    bwa_align()


    # Lancement picard_add_readgroups
    picard_add_readgroups()
    

    # lancement get_bed
    get_bed()
    
    
    # Lancement VCF call
    vcf_call_all_samples()

    
    # Lancement la creation de la base de donnee
    build_snpeff_db(db_name="pf_3D7_snpEff_db", 
        snpeff_config="pf_3D7_snpEff_db",
        output_dir="output")
        
    # Lancement l'annotation
    annotated = annotate_all_vcfs()

    
    # Lancement vartype
    BASE_DIR = Path(__file__).resolve().parent
    env_path = BASE_DIR / "miniconda3" / "envs" / "pipeline_env"
    logs_dir = BASE_DIR / "logs"
    
    results = run_vartype_all(
        input_dir=BASE_DIR / "output" / "vcf_annotated",
        output_dir=BASE_DIR / "output" / "vartype",
        env_path=env_path,
        logs_dir=logs_dir)
    
    print("\n[PIPELINE] Résultats Vartype :")
    for sample, files in results.items():
        print(f"  - {sample}: {files}")
    
        
    # Lancement coverage
    get_coverage(bam_dir="output/bam_picard_readgroups",
        output_dir="output/samtools_coverage", logs_dir="logs")


    # Lancement wt_coverage
    ref = "pf_3D7_Ref/mars_pf_ref.fasta"
    gff = "pf_3D7_Ref/mars_pf.gff"
    voi = "pf_3D7_Ref/voinew3.csv"
    results = run_wt_cov(ref, gff, voi)
    print("[PIPELINE] Résultats WT_cov générés :")
    for sample, path in results.items():
        print(f"  - {sample}: {path}")


    # Lancement Trim_Stats
    trim_stats_results = run_trim_stats()
    print(f"Trim_Stats terminé pour {len(trim_stats_results)} échantillons.\n")


    # Lancement run_reads_merge
    merged_file = run_reads_merge()
    print(f"Matrice finale de lecture générée : {merged_file}")
    

    
    # Lancement run_vcf_to_df
    run_vcf_to_df()
    
    
    # Lancement run_csv_merge_with_report
    run_csv_merge_with_report()
    
    
    # Lancement compute_sample_vaf
    compute_sample_vaf()
    
    # Lancement run_SVAF_merge
    run_SVAF_merge()
    
    
    #Lancement run_snpfilter
    run_snpfilter()
    
    
    #Lancement run_summary_merge
    run_summary_merge()


    # Lancement run_introns_merge
    run_introns_merge()
    
    
    # Lancement run_summary
    run_summary()
    
    
    # Lancement run_dataviz_reportable_snps
    run_dataviz_reportable_snps()
    
    
    # Lancement run_dataviz_novel_snps
    run_DataViz_Novel_snps()
    
    
    # Lancement run_replace_mutations
    run_replace_mutations()
    
    
    # Lancement filtered_summary_merge
    filtered_summary_merge()
    
    
    # Lancement filter_empty_pos
    filter_empty_pos()
    
    # Lancement run_haplotypes
    run_haplotypes()
    
    
    # Lancement filter_haplotypes
    filter_haplotypes()
    
    
    # génération du fichier CSV combiné des haplotypes
    run_combined_haplotypes()
    
    """
    # Lancement generate_final_report_by_site
    generate_final_report_by_site()
    
    
    generate_final_report_by_site_0()

    total_elapsed = time.time() - global_start
    print(f"\n[PIPELINE] Analyses terminées avec succès en {total_elapsed:.2f} sec.")
    print("[PIPELINE] Fin du pipeline.")   
    