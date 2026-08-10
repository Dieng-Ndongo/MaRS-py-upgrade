import re
import subprocess
import zipfile
import os
from datetime import datetime
from pathlib import Path

import streamlit as st

# ─────────────────────────────────────────────
# Cartographie des 37 étapes → phases visuelles
# ─────────────────────────────────────────────
# Format : (step_min, step_max, label_phase)
PHASES = [
    (1,  2,  "QC pré-trim"),
    (3,  3,  "Trimming"),
    (4,  5,  "QC post-trim"),
    (6,  7,  "Alignement"),
    (8,  8,  "Picard RG"),
    (9,  9,  "BED file"),
    (10, 10, "VCF calling"),
    (11, 12, "SnpEff"),
    (13, 13, "VarType"),
    (14, 15, "Coverage"),
    (16, 17, "Reads & Stats"),
    (18, 19, "VCF → DataFrame"),
    (20, 21, "VAF calcul"),
    (22, 25, "Filtres"),
    (26, 28, "DataViz"),
    (29, 31, "Haplotypes"),
    (32, 33, "Haplo filtrés"),
    (34, 37, "Rapports finaux"),
]

# Labels complets des 37 étapes (index 0 = step 1)
STEP_LABELS = [
    "QC pré-trimming (FastQC)",
    "MultiQC pré-trimming",
    "Trimming (BBduk)",
    "QC post-trimming (FastQC)",
    "MultiQC post-trimming",
    "BWA index",
    "BWA align",
    "Picard add readgroups",
    "Get BED",
    "VCF call (freebayes)",
    "Build snpEff database",
    "Annotate VCFs (snpEff)",
    "Run VarType",
    "Coverage (samtools)",
    "WT Coverage",
    "Trim Stats",
    "Reads merge",
    "VCF → DataFrame",
    "CSV merge + report",
    "Compute sample VAF",
    "SVAF merge",
    "SNP filter",
    "Summary merge",
    "Introns merge",
    "Run summary",
    "DataViz reportable SNPs",
    "DataViz novel SNPs",
    "Replace mutations",
    "Filtered summary merge",
    "Filter empty positions",
    "Run haplotypes",
    "Filter haplotypes",
    "Combined haplotypes",
    "Final report by site",
    "Final report by site (MIX→MT)",
    "Pondered report (MIX→MT) — 1",
    "Pondered report (MIX→MT) — 2",
]

TOTAL_STEPS = 37
_CSS = """
<style>
.trk-section-label {
    font-size: 0.78rem;
    font-weight: 700;
    text-transform: uppercase;
    letter-spacing: 0.12em;
    color: #4b5563;
    margin: 0 0 0.9rem;
}
.trk-prog-header {
    display: flex;
    justify-content: space-between;
    align-items: baseline;
    margin-bottom: 7px;
}
.trk-prog-label { font-size: 0.85rem; font-weight: 600; color: #374151; }
.trk-prog-pct   { font-size: 0.95rem; font-weight: 700; color: #185FA5; }
.trk-prog-pct.done { color: #3B6D11; }
.trk-prog-outer {
    background: #e5e7eb;
    border-radius: 999px;
    height: 10px;
    width: 100%;
    overflow: hidden;
    margin-bottom: 1.1rem;
}
.trk-prog-fill {
    height: 100%;
    border-radius: 999px;
    background: #378ADD;
    transition: width 0.4s ease;
}
.trk-prog-fill.done { background: #639922; }
.trk-metrics {
    display: flex;
    border: 1px solid #d1d5db;
    border-radius: 8px;
    overflow: hidden;
    margin-bottom: 1.1rem;
}
.trk-metric {
    flex: 1;
    padding: 12px 16px;
    border-right: 1px solid #d1d5db;
}
.trk-metric:last-child { border-right: none; }
.trk-metric-label {
    font-size: 0.75rem;
    font-weight: 600;
    color: #6b7280;
    margin-bottom: 5px;
    text-transform: uppercase;
    letter-spacing: 0.06em;
}
.trk-metric-value {
    font-size: 1.05rem;
    font-weight: 700;
    color: #111827;
    line-height: 1.2;
}
.trk-step-box {
    display: flex;
    align-items: center;
    gap: 12px;
    background: #f0f6ff;
    border-left: 4px solid #378ADD;
    border-radius: 0 8px 8px 0;
    padding: 12px 16px;
    margin-bottom: 1.2rem;
}
.trk-step-badge {
    background: #dbeafe;
    color: #185FA5;
    border-radius: 50%;
    width: 34px; height: 34px;
    display: flex; align-items: center; justify-content: center;
    font-size: 0.85rem; font-weight: 700;
    flex-shrink: 0;
}
.trk-step-name  { font-size: 0.95rem; font-weight: 700; color: #111827; }
.trk-step-phase { font-size: 0.82rem; font-weight: 500; color: #4b5563; margin-top: 3px; }
.trk-phases-label {
    font-size: 0.78rem;
    font-weight: 700;
    text-transform: uppercase;
    letter-spacing: 0.1em;
    color: #6b7280;
    margin-bottom: 9px;
}
.trk-phases-grid {
    display: grid;
    grid-template-columns: repeat(6, 1fr);
    gap: 6px;
    margin-bottom: 1.1rem;
}
.trk-phase {
    border-radius: 6px;
    padding: 8px 5px;
    font-size: 0.74rem;
    font-weight: 500;
    text-align: center;
    line-height: 1.5;
    border: 1px solid transparent;
}
.trk-phase-done    { background: #EAF3DE; border-color: #C0DD97; color: #2d5a0e; }
.trk-phase-active  {
    background: #E6F1FB; border-color: #378ADD; color: #0C447C;
    font-weight: 700;
    box-shadow: 0 0 0 2px #B5D4F4;
}
.trk-phase-pending { background: #f9fafb; border-color: #e5e7eb; color: #6b7280; }
.trk-phase-icon { display: block; font-size: 0.8rem; margin-bottom: 3px; }
.trk-spinner {
    display: inline-block;
    width: 12px; height: 12px;
    border: 2px solid #B5D4F4;
    border-top-color: #185FA5;
    border-radius: 50%;
    animation: trk-spin 0.75s linear infinite;
    margin-bottom: 3px;
    vertical-align: middle;
}
@keyframes trk-spin { to { transform: rotate(360deg); } }
</style>
"""


# ─────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────

def _get_phase_idx(step: int) -> int:
    for i, (smin, smax, _) in enumerate(PHASES):
        if smin <= step <= smax:
            return i
    return len(PHASES) - 1


def _get_phase_label(step: int) -> str:
    for smin, smax, label in PHASES:
        if smin <= step <= smax:
            return label
    return ""


def _parse_log(log_text: str) -> tuple[int, int]:
    pattern = re.compile(r"(?:Step|Etape)\s*(\d+)\s*/\s*(\d+)", re.IGNORECASE)
    matches = pattern.findall(log_text)
    if matches:
        last = matches[-1]
        return int(last[0]), int(last[1])
    return 0, TOTAL_STEPS


def _elapsed_str(seconds: float) -> str:
    s = int(seconds)
    if s < 60:
        return f"{s} s"
    return f"{s // 60} min {s % 60} s"


def _eta_str(elapsed: float, percent: int) -> str:
    if percent <= 0 or percent >= 100:
        return "—"
    total_est = elapsed / (percent / 100)
    remaining = total_est - elapsed
    if remaining < 0:
        return "—"
    return _elapsed_str(remaining)


def _build_tracker_html(
    percent: int,
    elapsed: float,
    eta: str,
    n_samp,
    current_step: int,
    total_steps: int,
    step_lbl: str,
    phase_lbl: str,
    phase_idx: int,
) -> str:
    """Construit tout le bloc HTML du tracker en une seule chaîne."""

    bar_width = max(percent, 1)
    fill_cls  = "trk-prog-fill done" if percent == 100 else "trk-prog-fill"
    pct_cls   = "trk-prog-pct done"  if percent == 100 else "trk-prog-pct"

    # ── Barre ──
    html = _CSS
    html += f"""
<div class="trk-section-label">Analyse bioinformatique en cours</div>
<div class="trk-prog-header">
  <span class="trk-prog-label">Progression</span>
  <span class="{pct_cls}">{percent}%</span>
</div>
<div class="trk-prog-outer">
  <div class="{fill_cls}" style="width:{bar_width}%;"></div>
</div>
"""

    # ── Métriques ──
    html += f"""
<div class="trk-metrics">
  <div class="trk-metric">
    <div class="trk-metric-label">Durée</div>
    <div class="trk-metric-value">{_elapsed_str(elapsed)}</div>
  </div>
  <div class="trk-metric">
    <div class="trk-metric-label">Temps restant</div>
    <div class="trk-metric-value">{eta}</div>
  </div>
  <div class="trk-metric">
    <div class="trk-metric-label">Échantillons</div>
    <div class="trk-metric-value">{n_samp}</div>
  </div>
  <div class="trk-metric">
    <div class="trk-metric-label">Étape</div>
    <div class="trk-metric-value">{current_step} / {total_steps}</div>
  </div>
</div>
"""

    # ── Étape courante ──
    if current_step > 0:
        html += f"""
<div class="trk-step-box">
  <div class="trk-step-badge">{current_step}</div>
  <div>
    <div class="trk-step-name">{step_lbl}</div>
    <div class="trk-step-phase">Phase : {phase_lbl}</div>
  </div>
</div>
"""

    # ── Phases ──
    html += '<div class="trk-phases-label">Phases du pipeline</div>'
    html += '<div class="trk-phases-grid">'
    for i, (smin, smax, label) in enumerate(PHASES):
        if percent == 100 or i < phase_idx:
            cls  = "trk-phase trk-phase-done"
            icon = "✓"
        elif i == phase_idx:
            cls  = "trk-phase trk-phase-active"
            html += (
                f'<div class="{cls}" title="Étapes {smin}–{smax}">'
                f'<div><span class="trk-spinner"></span></div>{label}'
                f'</div>'
            )
            continue
        else:
            cls  = "trk-phase trk-phase-pending"
            icon = "·"
        html += (
            f'<div class="{cls}" title="Étapes {smin}–{smax}">'
            f'<span class="trk-phase-icon">{icon}</span>{label}'
            f'</div>'
        )
    html += "</div>"

    return html


# ─────────────────────────────────────────────
# Rendu principal
# ─────────────────────────────────────────────

def render_progress_tracker(get_run_paths_fn, save_history_fn, notify_fn):
    """
    À appeler dans app.py à la place du bloc « if st.session_state.get("running"): ».

    Paramètres
    ----------
    get_run_paths_fn  : callable → dict avec clés "logs", "output", "zip"
    save_history_fn   : callable(run_id, names, status, duration, zip_path)
    notify_fn         : callable(run_id, n_samples, duration_sec, success)
    """
    names  = st.session_state.get("sample_names", [])
    run_id = st.session_state.get("run_id")

    # ── Bouton Arrêter ────────────────────────
    col_warn, col_stop = st.columns([3, 1])
    with col_warn:
        st.warning("Analyse en cours — upload et lancement temporairement désactivés.")
    with col_stop:
        if st.button("Arrêter", use_container_width=True):
            id_res = subprocess.run(
                ["docker", "ps", "--filter", "ancestor=bioinfo_pipeline", "--format", "{{.ID}}"],
                capture_output=True, text=True,
            )
            cid = id_res.stdout.strip()
            if cid:
                subprocess.run(["docker", "stop", cid], capture_output=True)
            st.session_state["running"] = False
            elapsed = int(datetime.now().timestamp() - st.session_state.get("launch_time", 0))
            _nfh = st.session_state.get("names_for_history", names)
            if run_id:
                save_history_fn(run_id, _nfh, "failed", elapsed, None)
                notify_fn(run_id, len(_nfh), elapsed, success=False)
            st.warning("Analyse arrêtée.")
            st.rerun()

    # ── Lecture du log ────────────────────────
    log_file    = Path(st.session_state.get("log_file", ""))
    launch_time = st.session_state.get("launch_time", 0)
    elapsed     = datetime.now().timestamp() - launch_time


    logs = ""
    if log_file.exists():
        with open(log_file, "r", errors="replace") as f:
            logs = f.read()

    current_step, total_steps = _parse_log(logs)
    st.session_state["current_step"] = current_step  # ← ajouter cette ligne
    percent   = int(current_step / total_steps * 100) if total_steps > 0 and current_step > 0 else 0
    phase_idx = _get_phase_idx(current_step) if current_step > 0 else -1
    step_lbl  = STEP_LABELS[current_step - 1] if 1 <= current_step <= TOTAL_STEPS else "Démarrage…"
    phase_lbl = _get_phase_label(current_step) if current_step > 0 else "—"
    eta       = _eta_str(elapsed, percent)
    n_samp    = len(names) if isinstance(names, list) else "?"

    # ── Rendu HTML unique ─────────────────────
    st.markdown(
        _build_tracker_html(
            percent=percent,
            elapsed=elapsed,
            eta=eta,
            n_samp=n_samp,
            current_step=current_step,
            total_steps=total_steps,
            step_lbl=step_lbl,
            phase_lbl=phase_lbl,
            phase_idx=phase_idx,
        ),
        unsafe_allow_html=True,
    )

    # ── Logs (50 dernières lignes) ────────────
    if logs.strip():
        with st.expander("Logs pipeline (50 dernières lignes)", expanded=False):
            st.code("\n".join(logs.splitlines()[-50:]), language="bash")
    else:
        st.info("Démarrage du container Docker… En attente des premiers logs.")

    # ── Détection de fin ──────────────────────
    if elapsed > 10:
        result = subprocess.run(
            ["docker", "ps", "--filter", "ancestor=bioinfo_pipeline", "--format", "{{.ID}}"],
            capture_output=True, text=True,
        )
        docker_running = result.stdout.strip() != ""

        if not docker_running:
            st.session_state["running"] = False
            fatal_patterns = [
                "command not found", "permission denied", "no such file or directory",
                "killed", "oom", "segmentation fault", "traceback",
                "exception:", "pipeline failed", "pipeline aborted",
            ]
            has_fatal = any(pat in logs.lower() for pat in fatal_patterns)
            duration  = int(elapsed)
            _nfh      = st.session_state.get("names_for_history", names)

            if has_fatal:
                st.session_state["pipeline_done"] = False
                if run_id:
                    save_history_fn(run_id, _nfh, "failed", duration, None)
                    notify_fn(run_id, len(_nfh), duration, success=False)
            else:
                st.session_state["pipeline_done"] = True
                paths_zip = get_run_paths_fn(run_id)
                zip_path  = None
                if paths_zip["output"].exists() and not st.session_state.get("zip_created"):
                    zip_path = paths_zip["zip"]
                    with zipfile.ZipFile(zip_path, "w") as zipf:
                        for root, _, files in os.walk(paths_zip["output"]):
                            for file in files:
                                full_path = os.path.join(root, file)
                                zipf.write(full_path, os.path.relpath(full_path, paths_zip["output"]))
                    st.session_state["zip_created"]   = True
                    st.session_state["show_download"] = True
                if run_id:
                    save_history_fn(run_id, _nfh, "success", duration, zip_path)
                    notify_fn(run_id, len(_nfh), duration, success=True)
            st.rerun()
