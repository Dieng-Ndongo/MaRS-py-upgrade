#!/bin/bash
# One-shot fix for "permission denied" errors from the mars-streamlit service.
#
# Root cause: being in the `docker` group is already root-equivalent on the host
# (any docker-group member can trivially become root via a bind-mounted container),
# so `mars` being in that group is not itself a problem. The actual bug is usually
# one of:
#   1. `mars` was added to the `docker` group *after* mars-streamlit was already
#      running — systemd only reads supplementary groups when a unit starts, so the
#      running process still has the old group list until it's restarted.
#   2. The repo tree (or runs/) drifted away from mars:mars ownership (e.g. a file
#      edited as root).
#   3. A stale `bioinfo_pipeline` image was built manually (README step 4) without
#      the right --build-arg PUID/PGID, baking in a UID that doesn't match `mars`.
#
# This script re-applies all three fixes and restarts the service. Run as root:
#   sudo bash deploy/fix-server-permissions.sh

set -euo pipefail

if [ "$(id -u)" -ne 0 ]; then
    echo "[ERREUR] Ce script doit être exécuté en root (sudo bash deploy/fix-server-permissions.sh)." >&2
    exit 1
fi

REPO_DIR="/opt/mars-py-upgrade"

echo "=== État actuel ==="
id mars || { echo "[ERREUR] L'utilisateur 'mars' n'existe pas." >&2; exit 1; }
getent group docker || echo "[ATTENTION] Le groupe 'docker' n'existe pas — Docker est-il installé ?"

echo
echo "=== 1. Ajout de mars au groupe docker (idempotent) ==="
usermod -aG docker mars

echo
echo "=== 2. Réapplication de la propriété mars:mars sur $REPO_DIR ==="
if [ -d "$REPO_DIR" ]; then
    chown -R mars:mars "$REPO_DIR"
else
    echo "[ATTENTION] $REPO_DIR introuvable, étape ignorée."
fi

echo
echo "=== 3. Suppression de l'image bioinfo_pipeline existante (forcer un rebuild propre) ==="
docker image rm -f bioinfo_pipeline 2>/dev/null || echo "(aucune image à supprimer)"
rm -f "$REPO_DIR/venv/.docker_hash" 2>/dev/null || true

echo
echo "=== 4. Redémarrage de mars-streamlit (prise en compte du nouveau groupe) ==="
systemctl restart mars-streamlit

echo
echo "=== 5. Vérification ==="
sleep 2
if sudo -u mars docker ps >/dev/null 2>&1; then
    echo "[OK] mars peut parler au démon Docker sans permission denied."
else
    echo "[ERREUR] mars ne peut toujours pas accéder à Docker — vérifiez manuellement." >&2
fi

echo
echo "--- Dernières lignes du journal mars-streamlit ---"
journalctl -u mars-streamlit -n 20 --no-pager
