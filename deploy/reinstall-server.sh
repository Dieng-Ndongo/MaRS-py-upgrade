#!/bin/bash
# Full clean reinstall of the mars-streamlit app on this server.
#
# Resets the repo to a pristine `production` checkout (discards any drift, stray
# root-owned files, half-built venv, stale Docker image/containers) while explicitly
# preserving runs/ (past analysis results) and .streamlit/secrets.toml (APP_PASSWORD
# etc., not reproducible from git history — both are gitignored, see .gitignore:170,196).
#
# Run ONCE, as root:
#   sudo bash /opt/mars-py-upgrade/deploy/reinstall-server.sh

set -euo pipefail

if [ "$(id -u)" -ne 0 ]; then
    echo "[ERREUR] Exécutez ce script avec sudo (sudo bash deploy/reinstall-server.sh)." >&2
    exit 1
fi

REPO_DIR="/opt/mars-py-upgrade"
BRANCH="production"

if [ ! -d "$REPO_DIR" ]; then
    echo "[ERREUR] $REPO_DIR introuvable." >&2
    exit 1
fi

echo "=== 1. Arrêt du service ==="
systemctl stop mars-streamlit.timer 2>/dev/null || true
systemctl stop mars-streamlit.service 2>/dev/null || true

echo
echo "=== 2. Suppression des conteneurs/image Docker existants ==="
docker ps -aq --filter "ancestor=bioinfo_pipeline" | xargs -r docker rm -f
docker image rm -f bioinfo_pipeline 2>/dev/null || echo "(aucune image à supprimer)"

echo
echo "=== 3. Réinitialisation du dépôt sur '$BRANCH' (conservation de runs/ et des secrets) ==="
cd "$REPO_DIR"
sudo -u mars git fetch origin
sudo -u mars git checkout "$BRANCH"
sudo -u mars git reset --hard "origin/$BRANCH"
sudo -u mars git clean -fdx -e runs -e .streamlit/secrets.toml

echo
echo "=== 4. Propriété correcte sur tout le dépôt ==="
chown -R mars:mars "$REPO_DIR"

echo
echo "=== 5. Réactivation du service (le premier démarrage recrée le venv et reconstruit"
echo "        l'image Docker — patientez quelques minutes) ==="
systemctl daemon-reload
systemctl enable --now mars-streamlit.timer
systemctl start mars-streamlit

echo
echo "=== 6. État actuel ==="
sleep 3
systemctl status mars-streamlit --no-pager || true

echo
echo "Suivez la reconstruction avec : journalctl -u mars-streamlit -f"
echo "Une fois prêt, vérifiez : ls $REPO_DIR/runs   (doit encore contenir vos anciens runs)"
