#!/bin/bash

# Répertoire du script (fonctionne quel que soit l'endroit où le repo est cloné)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV_DIR="$SCRIPT_DIR/venv"

# Créer le venv + installer les dépendances si absent
if [ ! -d "$VENV_DIR" ]; then
    echo "[INFO] Création du venv Streamlit..."
    python3 -m venv "$VENV_DIR"
    "$VENV_DIR/bin/pip" install --upgrade pip
    "$VENV_DIR/bin/pip" install -r "$SCRIPT_DIR/requirements.txt"
fi

source "$VENV_DIR/bin/activate"
streamlit run "$SCRIPT_DIR/app.py"
