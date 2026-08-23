#!/bin/bash

# Répertoire du script (fonctionne quel que soit l'endroit où le repo est cloné)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV_DIR="$SCRIPT_DIR/venv"
REQ_FILE="$SCRIPT_DIR/requirements.txt"
REQ_HASH_FILE="$VENV_DIR/.req_hash"

# Créer le venv si absent
if [ ! -d "$VENV_DIR" ]; then
    if ! command -v python3.12 >/dev/null 2>&1; then
        echo "[ERREUR] python3.12 est requis mais introuvable." >&2
        echo "Installez-le avec :" >&2
        echo "  sudo apt update && sudo apt install -y software-properties-common" >&2
        echo "  sudo add-apt-repository -y ppa:deadsnakes/ppa" >&2
        echo "  sudo apt update && sudo apt install -y python3.12 python3.12-venv" >&2
        exit 1
    fi
    echo "[INFO] Création du venv Streamlit..."
    python3.12 -m venv "$VENV_DIR"
    "$VENV_DIR/bin/pip" install --upgrade pip
fi

# Réinstaller les dépendances si requirements.txt a changé
CURRENT_HASH=$(md5sum "$REQ_FILE" | cut -d' ' -f1)
SAVED_HASH=$(cat "$REQ_HASH_FILE" 2>/dev/null || echo "")

if [ "$CURRENT_HASH" != "$SAVED_HASH" ]; then
    echo "[INFO] requirements.txt modifié — mise à jour des dépendances..."
    "$VENV_DIR/bin/pip" install --upgrade pip
    "$VENV_DIR/bin/pip" install -r "$REQ_FILE"
    echo "$CURRENT_HASH" > "$REQ_HASH_FILE"
    echo "[INFO] Dépendances à jour."
fi

# Construire l'image Docker du pipeline si absente ou si les sources ont changé
DOCKER_HASH_FILE="$VENV_DIR/.docker_hash"
DOCKER_HASH=$(cat "$SCRIPT_DIR/Dockerfile" "$SCRIPT_DIR/environment.yml" | md5sum | cut -d' ' -f1)
SAVED_DOCKER_HASH=$(cat "$DOCKER_HASH_FILE" 2>/dev/null || echo "")

if ! docker image inspect bioinfo_pipeline >/dev/null 2>&1 || [ "$DOCKER_HASH" != "$SAVED_DOCKER_HASH" ]; then
    echo "[INFO] Construction de l'image Docker bioinfo_pipeline..."
    docker build -t bioinfo_pipeline "$SCRIPT_DIR"
    echo "$DOCKER_HASH" > "$DOCKER_HASH_FILE"
fi

# Éviter le prompt d'onboarding Streamlit (bloquant en exécution non-interactive)
mkdir -p "$HOME/.streamlit"
if [ ! -f "$HOME/.streamlit/credentials.toml" ]; then
    printf '[general]\nemail = ""\n' > "$HOME/.streamlit/credentials.toml"
fi

source "$VENV_DIR/bin/activate"
exec streamlit run "$SCRIPT_DIR/app.py"
