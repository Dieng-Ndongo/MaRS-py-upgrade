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
# INIT SESSION STATE
# ========================
if "pipeline_done" not in st.session_state:
    st.session_state["pipeline_done"] = False

if "zip_created" not in st.session_state:
    st.session_state["zip_created"] = False

# ========================
# STYLE CSS
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
</style>
""", unsafe_allow_html=True)

# ========================
# HEADER
# ========================
st.markdown('<div class="main-title">🧬 MaRS-py-upgrade</div>', unsafe_allow_html=True)
st.markdown('<div class="subtitle">Pipeline FASTQ → Rapport automatisé</div>', unsafe_allow_html=True)
st.divider()

# ========================
# UPLOAD
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

    st.session_state["pipeline_done"] = False
    st.session_state["zip_created"] = False

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
    info_box.empty()

    if process.returncode == 0:
        progress_bar.progress(1.0)
        status_text.success("✅ Pipeline terminé avec succès !")
        st.session_state["pipeline_done"] = True
        st.balloons()
    else:
        status_text.error("❌ Une erreur est survenue.")

output_dir = Path.home() / "pipeline" / "output"
zip_path = Path.home() / "pipeline" / "resultats_pipeline.zip"

if "download_ready" not in st.session_state:
    st.session_state["download_ready"] = False

if output_dir.exists() and st.session_state.get("pipeline_done"):

    st.divider()
    st.subheader("📥 Résultats")

    # 1️⃣ Créer ZIP automatiquement (UNE SEULE FOIS)
    if not st.session_state.get("zip_created", False):
        with zipfile.ZipFile(zip_path, "w") as zipf:
            for root, dirs, files in os.walk(output_dir):
                for file in files:
                    zipf.write(
                        os.path.join(root, file),
                        arcname=os.path.relpath(os.path.join(root, file), output_dir)
                    )

        st.session_state["zip_created"] = True

    # 2️⃣ Bouton normal (NE déclenche pas download)
    if st.button("⬇ Télécharger les résultats"):
        st.session_state["download_ready"] = True
       
    # 3️⃣ Download seulement après clic utilisateur
    if st.session_state["download_ready"] and zip_path.exists():

        with open(zip_path, "rb") as f:
            downloaded = st.download_button(
                "📦 Cliquer pour télécharger le fichier",
                f,
                file_name="resultats_pipeline.zip",
                use_container_width=True
            )

        # 🔥 RESET après téléchargement
        if downloaded:
            st.session_state["download_ready"] = False