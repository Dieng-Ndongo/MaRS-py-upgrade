import streamlit as st
import subprocess
import os
from pathlib import Path
import re
import zipfile

# ========================
# CONFIGURATION PAGE
# ========================
st.set_page_config(
    page_title="Pipeline FASTQ → Haplotypes",
    layout="wide",
    page_icon="🧬"
)

# ========================
# STYLE CSS PERSONNALISÉ
# ========================
st.markdown("""
<style>
    .main-title {
        font-size: 40px;
        font-weight: bold;
        color: #4CAF50;
    }
    .subtitle {
        font-size: 18px;
        color: #888;
    }
    .success-box {
        padding: 10px;
        border-radius: 10px;
        background-color: #e8f5e9;
        border-left: 5px solid #4CAF50;
    }
</style>
""", unsafe_allow_html=True)

# ========================
# HEADER
# ========================
st.markdown('<div class="main-title">🧬 MaRS-py-upgrade</div>', unsafe_allow_html=True)
st.markdown('<div class="subtitle">Pipeline FASTQ → Rapport automatisé</div>', unsafe_allow_html=True)
st.divider()

# ========================
# UPLOAD DES FASTQ
# ========================
col1, col2 = st.columns([2, 1])

with col1:
    st.subheader("📤 Import des fichiers FASTQ")

    uploaded_files = st.file_uploader(
        "Glissez-déposez vos fichiers (.fastq.gz)",
        type=["fastq.gz"],
        accept_multiple_files=True
    )

    input_dir = Path.home() / "pipeline" / "data"
    input_dir.mkdir(parents=True, exist_ok=True)

    if uploaded_files:
        st.markdown(f"### ✅ {len(uploaded_files)} fichier(s) importé(s)")

        for file in uploaded_files:
            file_path = input_dir / file.name
            with open(file_path, "wb") as f:
                f.write(file.getbuffer())

            st.markdown(f"✔ **{file.name}**")

with col2:
    st.subheader("⚙️ Pipeline")
    st.info("Cliquez pour lancer l’analyse")

    run_pipeline = st.button("▶ Lancer", use_container_width=True)

# ========================
# EXECUTION PIPELINE
# ========================
if run_pipeline:

    st.divider()
    st.subheader("🚀 Exécution du pipeline")

    info_box = st.info("🧠 Analyse en cours... veuillez patienter")

    progress_bar = st.progress(0)
    status_text = st.empty()
    log_container = st.empty()

    docker_cmd = [
        "docker", "run", "--rm",
        "-v", f"{Path.home()}/pipeline:/app",
        "bioinfo_pipeline"
    ]

    logs = ""
    process = subprocess.Popen(
        docker_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True
    )

    step_pattern = re.compile(r"(?:Step|Etape)\s*(\d+)\s*/\s*(\d+)", re.IGNORECASE)

    for line in iter(process.stdout.readline, ""):
        logs += line

        # Logs stylés
        log_container.code("\n".join(logs.splitlines()[-15:]), language="bash")

        match = step_pattern.search(line)
        if match:
            current = int(match.group(1))
            total = int(match.group(2))

            if total > 0:
                progress = min(current / total, 1.0)

                progress_bar.progress(progress)

                percent = int(progress * 100)
                status_text.markdown(f"""
                ⏳ **Étape {current}/{total}**
                
                🔄 Progression : **{percent}%**
                """)

    process.wait()

    info_box.empty()  # ✅ enlève le message "en cours"

    if process.returncode == 0:
        progress_bar.progress(1.0)
        status_text.success("✅ Pipeline terminé avec succès !")
        st.balloons()  # optionnel 🎉
    else:
        status_text.error("❌ Une erreur est survenue.")

# ========================
# TELECHARGEMENT
# ========================
output_dir = Path.home() / "pipeline" / "output"
zip_path = Path.home() / "pipeline" / "resultats_pipeline.zip"

if output_dir.exists():
    st.divider()
    st.subheader("📥 Résultats")

    with zipfile.ZipFile(zip_path, "w") as zipf:
        for root, dirs, files in os.walk(output_dir):
            for file in files:
                zipf.write(
                    os.path.join(root, file),
                    arcname=os.path.relpath(os.path.join(root, file), output_dir)
                )

    with open(zip_path, "rb") as f:
        st.download_button(
            "⬇ Télécharger les résultats",
            f,
            file_name="resultats_pipeline.zip",
            use_container_width=True
        )

