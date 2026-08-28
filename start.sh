#!/bin/bash
set -euo pipefail

# ============================================================
# MaRS CIGASS Streamlit launcher
# Production user: mars
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

VENV_DIR="$SCRIPT_DIR/venv"
REQ_FILE="$SCRIPT_DIR/requirements.txt"

# Runtime state — never store this inside venv
STATE_DIR="$SCRIPT_DIR/.runtime"
REQ_HASH_FILE="$STATE_DIR/.req_hash"
DOCKER_HASH_FILE="$STATE_DIR/.docker_hash"

# Streamlit local configuration
STREAMLIT_DIR="$SCRIPT_DIR/.streamlit"
CREDENTIALS_FILE="$STREAMLIT_DIR/credentials.toml"

# Docker image
IMAGE_NAME="bioinfo_pipeline"

# ============================================================
# Safety: this script MUST run as mars
# ============================================================

if [ "$(id -un)" != "mars" ]; then
    echo "[ERREUR] start.sh doit être exécuté avec l'utilisateur 'mars'." >&2
    echo "[INFO] Utilisez :" >&2
    echo "  sudo -u mars $SCRIPT_DIR/start.sh" >&2
    exit 1
fi

# ============================================================
# Validate required files
# ============================================================

if [ ! -f "$REQ_FILE" ]; then
    echo "[ERREUR] requirements.txt introuvable : $REQ_FILE" >&2
    exit 1
fi

if [ ! -f "$SCRIPT_DIR/Dockerfile" ]; then
    echo "[ERREUR] Dockerfile introuvable : $SCRIPT_DIR/Dockerfile" >&2
    exit 1
fi

if [ ! -f "$SCRIPT_DIR/environment.yml" ]; then
    echo "[ERREUR] environment.yml introuvable : $SCRIPT_DIR/environment.yml" >&2
    exit 1
fi

if [ ! -f "$SCRIPT_DIR/app.py" ]; then
    echo "[ERREUR] app.py introuvable : $SCRIPT_DIR/app.py" >&2
    exit 1
fi

# ============================================================
# Create runtime directories
# ============================================================

mkdir -p "$STATE_DIR"
mkdir -p "$STREAMLIT_DIR"

# Make sure mars owns runtime directories
chown mars:mars "$STATE_DIR" "$STREAMLIT_DIR" 2>/dev/null || true

# ============================================================
# Python / Streamlit virtual environment
# ============================================================

if [ ! -d "$VENV_DIR" ]; then

    if ! command -v python3.12 >/dev/null 2>&1; then
        echo "[ERREUR] python3.12 est requis mais introuvable." >&2
        echo "[INFO] Installez python3.12 et python3.12-venv." >&2
        exit 1
    fi

    echo "[INFO] Création du venv Streamlit..."

    python3.12 -m venv "$VENV_DIR"

    "$VENV_DIR/bin/python" -m pip install --upgrade pip
fi

# ============================================================
# Python dependencies
# ============================================================

CURRENT_REQ_HASH=$(md5sum "$REQ_FILE" | awk '{print $1}')
SAVED_REQ_HASH=""

if [ -f "$REQ_HASH_FILE" ]; then
    SAVED_REQ_HASH=$(cat "$REQ_HASH_FILE")
fi

if [ "$CURRENT_REQ_HASH" != "$SAVED_REQ_HASH" ]; then

    echo "[INFO] requirements.txt modifié — installation des dépendances..."

    "$VENV_DIR/bin/python" -m pip install --upgrade pip
    "$VENV_DIR/bin/python" -m pip install -r "$REQ_FILE"

    printf '%s\n' "$CURRENT_REQ_HASH" > "$REQ_HASH_FILE"

    echo "[INFO] Dépendances Python à jour."
else
    echo "[INFO] Dépendances Python déjà à jour."
fi

# ============================================================
# Streamlit configuration
# ============================================================

if [ ! -f "$CREDENTIALS_FILE" ]; then

    printf '[general]\nemail = ""\n' > "$CREDENTIALS_FILE"

    chmod 600 "$CREDENTIALS_FILE"
fi

# ============================================================
# Streamlit static metadata
# ============================================================

INDEX_HTML=$(find "$VENV_DIR/lib" \
    -path "*/site-packages/streamlit/static/index.html" \
    -type f \
    -print -quit 2>/dev/null || true)

if [ -n "$INDEX_HTML" ]; then

    if ! grep -q "CIGASS MaRS" "$INDEX_HTML"; then

        echo "[INFO] Personnalisation du HTML Streamlit..."

        sed -i \
            -e 's#<title>Streamlit</title>#<title>CIGASS MaRS · Pipeline</title>#' \
            -e 's#</title>#</title>\n    <meta property="og:title" content="CIGASS MaRS · Pipeline">#' \
            "$INDEX_HTML"

    fi
fi

# ============================================================
# Docker availability
# ============================================================

if ! command -v docker >/dev/null 2>&1; then
    echo "[ERREUR] Docker est introuvable." >&2
    exit 1
fi

if ! docker info >/dev/null 2>&1; then
    echo "[ERREUR] L'utilisateur mars ne peut pas accéder à Docker." >&2
    echo "[INFO] Vérifiez :" >&2
    echo "  sudo usermod -aG docker mars" >&2
    exit 1
fi

# ============================================================
# Docker image hash
#
# Include every file/directory used by Docker COPY.
# UID/GID are included because they affect the image.
# ============================================================

DOCKER_INPUT_HASH=$(
    {
        printf '%s\0' \
            "$SCRIPT_DIR/Dockerfile" \
            "$SCRIPT_DIR/environment.yml"

    } |
    sort -z |
    xargs -0 md5sum 2>/dev/null |
    md5sum |
    awk '{print $1}'
)

DOCKER_HASH="${DOCKER_INPUT_HASH}-$(id -u)-$(id -g)"

SAVED_DOCKER_HASH=""

if [ -f "$DOCKER_HASH_FILE" ]; then
    SAVED_DOCKER_HASH=$(cat "$DOCKER_HASH_FILE")
fi

# ============================================================
# Build Docker image when necessary
# ============================================================

if ! docker image inspect "$IMAGE_NAME" >/dev/null 2>&1 || \
   [ "$DOCKER_HASH" != "$SAVED_DOCKER_HASH" ]; then

    echo "[INFO] Construction de l'image Docker $IMAGE_NAME..."

    docker build \
        --build-arg PUID="$(id -u)" \
        --build-arg PGID="$(id -g)" \
        -t "$IMAGE_NAME" \
        "$SCRIPT_DIR"

    printf '%s\n' "$DOCKER_HASH" > "$DOCKER_HASH_FILE"

    echo "[INFO] Image Docker construite."
else
    echo "[INFO] Image Docker déjà à jour."
fi

# ============================================================
# Display runtime information
# ============================================================

echo ""
echo "============================================================"
echo " MaRS CIGASS"
echo "============================================================"
echo " User       : $(id -un)"
echo " UID        : $(id -u)"
echo " GID        : $(id -g)"
echo " Project    : $SCRIPT_DIR"
echo " Docker     : $IMAGE_NAME"
echo "============================================================"
echo ""

# ============================================================
# Start Streamlit
# ============================================================

echo "[INFO] Démarrage de Streamlit..."
echo "[INFO] Application : $SCRIPT_DIR/app.py"

exec "$VENV_DIR/bin/streamlit" run "$SCRIPT_DIR/app.py"
