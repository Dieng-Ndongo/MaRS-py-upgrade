#!/bin/bash

# Répertoire du script (fonctionne quel que soit l'endroit où le repo est cloné)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV_DIR="$SCRIPT_DIR/venv"

# Créer le venv + installer les dépendances si absent
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
    "$VENV_DIR/bin/pip" install -r "$SCRIPT_DIR/requirements.txt"
fi

source "$VENV_DIR/bin/activate"
streamlit run "$SCRIPT_DIR/app.py"
