import subprocess
from datetime import datetime, timedelta
from pathlib import Path

import streamlit as st


# ─────────────────────────────────────────────
# CSS injecté une seule fois
# ─────────────────────────────────────────────
_DASHBOARD_CSS = """
<style>
/* ── KPI cards ── */
.dash-kpi-grid {
    display: grid;
    grid-template-columns: repeat(4, 1fr);
    gap: 14px;
    margin-bottom: 24px;
}
.dash-kpi {
    background: #ffffff;
    border: 1px solid #e2e8f0;
    border-radius: 10px;
    padding: 18px 20px 14px;
    position: relative;
}
.dash-kpi-label {
    font-size: 0.75rem;
    font-weight: 600;
    text-transform: uppercase;
    letter-spacing: 0.06em;
    color: #6b7280;
    margin-bottom: 8px;
}
.dash-kpi-value {
    font-size: 2rem;
    font-weight: 700;
    color: #1a2a4a;
    line-height: 1;
}
.dash-kpi-delta {
    font-size: 0.75rem;
    margin-top: 6px;
    color: #6b7280;
}
.dash-kpi-icon {
    position: absolute;
    top: 16px;
    right: 18px;
    font-size: 1.4rem;
    opacity: 0.18;
}

/* ── Docker status badge ── */
.docker-badge {
    display: inline-flex;
    align-items: center;
    gap: 7px;
    padding: 5px 14px;
    border-radius: 20px;
    font-size: 0.8rem;
    font-weight: 600;
    margin-bottom: 20px;
}
.docker-ok   { background: #edfaf0; color: #155724; border: 1px solid #a8d5b5; }
.docker-warn { background: #fff8e6; color: #7a4f00; border: 1px solid #f0d080; }
.docker-err  { background: #fdf0f0; color: #721c24; border: 1px solid #f0b8b8; }
.docker-dot  { width: 8px; height: 8px; border-radius: 50%; }
.docker-dot-ok   { background: #28a745; }
.docker-dot-warn { background: #f0ad00; }
.docker-dot-err  { background: #dc3545; }

/* ── Section label (réutilise le style global) ── */
.dash-section {
    display: flex;
    align-items: center;
    gap: 10px;
    font-size: 0.8rem;
    font-weight: 700;
    text-transform: uppercase;
    letter-spacing: 0.08em;
    color: #888888;
    margin: 22px 0 14px;
}
.dash-section::after {
    content: '';
    flex: 1;
    height: 1px;
    background: #d8d8d8;
}

/* ── Last run card ── */
.last-run-card {
    background: #f8faff;
    border: 1px solid #cddcf5;
    border-left: 4px solid #1f70b8;
    border-radius: 8px;
    padding: 16px 20px;
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: 16px;
    margin-bottom: 20px;
}
.last-run-id {
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.8rem;
    color: #6b7280;
    margin-bottom: 4px;
}
.last-run-info {
    font-size: 0.88rem;
    color: #374151;
    font-weight: 500;
}
.last-run-meta {
    font-size: 0.78rem;
    color: #888888;
    margin-top: 4px;
}

/* ── Activity mini-bar chart ── */
.activity-bar-wrap {
    display: flex;
    align-items: flex-end;
    gap: 5px;
    height: 72px;
    margin-bottom: 6px;
}
.act-bar-col {
    flex: 1;
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: flex-end;
    gap: 2px;
}
.act-bar {
    width: 100%;
    border-radius: 3px 3px 0 0;
    min-height: 4px;
    transition: height 0.3s;
}
.act-bar-ok   { background: #28a745; }
.act-bar-fail { background: #dc3545; }
.act-bar-zero { background: #e2e8f0; min-height: 4px; }
.act-label {
    font-size: 0.62rem;
    color: #aaaaaa;
    text-align: center;
    white-space: nowrap;
}

/* ── Recent runs mini-list ── */
.mini-run {
    display: flex;
    align-items: center;
    justify-content: space-between;
    padding: 9px 0;
    border-bottom: 1px solid #f0f0f0;
    font-size: 0.83rem;
}
.mini-run:last-child { border-bottom: none; }
.mini-run-id {
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.75rem;
    color: #1a2a4a;
    font-weight: 600;
}
.mini-run-meta { color: #6b7280; font-size: 0.75rem; }
.mini-badge {
    font-size: 0.68rem;
    font-weight: 600;
    padding: 2px 8px;
    border-radius: 10px;
}
.mini-ok   { background: #edfaf0; color: #155724; }
.mini-fail { background: #fdf0f0; color: #721c24; }
.mini-run_ { background: #e8f2fb; color: #1a3a6a; }
</style>
"""


# ─────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────

def _docker_status() -> dict:
    """Retourne l'état Docker : disponible, image présente, container actif."""
    result = {"daemon": False, "image": False, "running": False, "container_id": ""}
    try:
        r = subprocess.run(["docker", "info"], capture_output=True, text=True, timeout=4)
        result["daemon"] = (r.returncode == 0)
    except Exception:
        return result
    try:
        r = subprocess.run(
            ["docker", "images", "bioinfo_pipeline", "--format", "{{.Repository}}"],
            capture_output=True, text=True, timeout=4,
        )
        result["image"] = bool(r.stdout.strip())
    except Exception:
        pass
    try:
        r = subprocess.run(
            ["docker", "ps", "--filter", "ancestor=bioinfo_pipeline", "--format", "{{.ID}}"],
            capture_output=True, text=True, timeout=4,
        )
        cid = r.stdout.strip()
        result["running"] = bool(cid)
        result["container_id"] = cid
    except Exception:
        pass
    return result


def _activity_last_n_days(history: list, n: int = 14) -> list:
    """
    Retourne une liste de n entrées (les n derniers jours) :
    [{"label": "J-13", "ok": 2, "fail": 1}, ...]
    """
    today = datetime.now().date()
    buckets: dict = {}
    for i in range(n):
        d = today - timedelta(days=n - 1 - i)
        buckets[d] = {"ok": 0, "fail": 0}

    for h in history:
        try:
            d = datetime.strptime(h.get("date", ""), "%Y-%m-%d %H:%M").date()
        except Exception:
            continue
        if d in buckets:
            if h.get("status") == "success":
                buckets[d]["ok"] += 1
            elif h.get("status") == "failed":
                buckets[d]["fail"] += 1

    result = []
    for i, (d, v) in enumerate(sorted(buckets.items())):
        if i == n - 1:
            label = "Auj."
        elif i == n - 2:
            label = "J-1"
        elif (n - 1 - i) % 7 == 0 or i == 0:
            label = f"J-{n - 1 - i}"
        else:
            label = ""
        result.append({"label": label, "ok": v["ok"], "fail": v["fail"], "date": d})
    return result


def _kpi_trend(history: list, field: str, status_filter=None) -> str:
    """Compare le nombre de runs cette semaine vs semaine précédente."""
    today = datetime.now().date()
    this_week  = today - timedelta(days=7)
    last_week  = today - timedelta(days=14)
    tw = lw = 0
    for h in history:
        try:
            d = datetime.strptime(h.get("date", ""), "%Y-%m-%d %H:%M").date()
        except Exception:
            continue
        match = (status_filter is None) or (h.get("status") == status_filter)
        if match:
            if d >= this_week:
                tw += 1
            elif d >= last_week:
                lw += 1
    if lw == 0:
        return ""
    delta = tw - lw
    if delta > 0:
        return f"▲ +{delta} vs sem. préc."
    if delta < 0:
        return f"▼ {delta} vs sem. préc."
    return "= stable vs sem. préc."


# ─────────────────────────────────────────────
# Rendu principal
# ─────────────────────────────────────────────

def render_dashboard(load_history_fn, get_run_paths_fn):
    """
    Point d'entrée à appeler dans app.py :

        if active_page == "home":
            render_dashboard(load_history, get_run_paths)
    """
    # Auto-refresh toutes les 30 secondes (uniquement sur cette page)
    from streamlit_autorefresh import st_autorefresh
    st_autorefresh(interval=30_000, key="dashboard_refresh")

    # Injection CSS
    st.markdown(_DASHBOARD_CSS, unsafe_allow_html=True)

    # ── Titre ──────────────────────────────────
    st.markdown(
        '<div class="section-label" style="margin-bottom:20px;">Tableau de bord</div>',
        unsafe_allow_html=True,
    )

    # ── Statut Docker ──────────────────────────
    ds = _docker_status()
    if not ds["daemon"]:
        badge_cls, dot_cls = "docker-err", "docker-dot-err"
        badge_msg = "Docker — démon non disponible"
    elif not ds["image"]:
        badge_cls, dot_cls = "docker-warn", "docker-dot-warn"
        badge_msg = "Docker actif — image bioinfo_pipeline introuvable"
    elif ds["running"]:
        badge_cls, dot_cls = "docker-warn", "docker-dot-warn"
        badge_msg = f"Container actif ({ds['container_id'][:12]})"
    else:
        badge_cls, dot_cls = "docker-ok", "docker-dot-ok"
        badge_msg = "Docker prêt — image bioinfo_pipeline disponible"

    st.markdown(
        f'<span class="docker-badge {badge_cls}">'
        f'<span class="docker-dot {dot_cls}"></span>{badge_msg}</span>',
        unsafe_allow_html=True,
    )

    # ── Chargement de l'historique ────────────
    history = load_history_fn()

    # ── KPI cards ─────────────────────────────
    total       = len(history)
    successes   = sum(1 for h in history if h.get("status") == "success")
    failures    = sum(1 for h in history if h.get("status") == "failed")
    total_samp  = sum(h.get("samples", 0) for h in history)
    success_pct = f"{round(successes / total * 100)}%" if total else "—"

    durations = [
        h["duration_sec"]
        for h in history
        if h.get("duration_sec") and h["duration_sec"] > 0
    ]
    avg_dur_str = (
        f"{int(sum(durations)/len(durations))//60} min "
        f"{int(sum(durations)/len(durations))%60} s"
        if durations else "—"
    )

    trend_runs    = _kpi_trend(history, "run_id")
    trend_success = _kpi_trend(history, "run_id", status_filter="success")

    kpis = [
        ("Runs totaux",          str(total),       "🔬", trend_runs),
        ("Taux de succès",       success_pct,      "✅", trend_success),
        ("Échantillons traités", str(total_samp),  "🧬", ""),
        ("Durée moyenne",        avg_dur_str,       "⏱", ""),
    ]

    st.markdown('<div class="dash-kpi-grid">', unsafe_allow_html=True)
    for label, value, icon, delta in kpis:
        delta_html = f'<div class="dash-kpi-delta">{delta}</div>' if delta else ""
        st.markdown(
            f"""
            <div class="dash-kpi">
                <span class="dash-kpi-icon">{icon}</span>
                <div class="dash-kpi-label">{label}</div>
                <div class="dash-kpi-value">{value}</div>
                {delta_html}
            </div>
            """,
            unsafe_allow_html=True,
        )
    st.markdown("</div>", unsafe_allow_html=True)

    # ── Layout 2 colonnes ─────────────────────
    col_left, col_right = st.columns([3, 2])

    # ── Dernier run ───────────────────────────
    with col_left:
        st.markdown('<div class="dash-section">Dernier run</div>', unsafe_allow_html=True)

        last = next(
            (h for h in history if h.get("status") in ("success", "failed")), None
        )
        if last:
            run_id    = last["run_id"]
            status    = last.get("status", "—")
            n_samp    = last.get("samples", "?")
            date_str  = last.get("date", "—")
            dur_s     = last.get("duration_sec")
            dur_str   = f"{dur_s//60} min {dur_s%60} s" if dur_s else "—"
            badge_ok  = status == "success"
            badge_lbl = "✓ Succès" if badge_ok else "✗ Échec"
            badge_col = "#edfaf0" if badge_ok else "#fdf0f0"
            badge_tc  = "#155724" if badge_ok else "#721c24"
            badge_bc  = "#a8d5b5" if badge_ok else "#f0b8b8"

            zip_path = Path(last["zip_path"]) if last.get("zip_path") else None

            st.markdown(
                f"""
                <div class="last-run-card">
                    <div>
                        <div class="last-run-id">{run_id}</div>
                        <div class="last-run-info">
                            <span style="background:{badge_col};color:{badge_tc};
                                  border:1px solid {badge_bc};border-radius:10px;
                                  padding:2px 10px;font-size:0.72rem;font-weight:600;">
                                {badge_lbl}
                            </span>
                        </div>
                        <div class="last-run-meta">
                            📅 {date_str} &nbsp;·&nbsp;
                            🧬 {n_samp} échantillon(s) &nbsp;·&nbsp;
                            ⏱ {dur_str}
                        </div>
                    </div>
                </div>
                """,
                unsafe_allow_html=True,
            )

            btn_c1, btn_c2, btn_c3 = st.columns(3)
            with btn_c1:
                if st.button("📋 Voir l'historique", key="dash_to_hist", width="stretch"):
                    st.session_state["active_page"] = "history"
                    st.rerun()
            with btn_c2:
                if st.button("📊 Voir le QC", key="dash_to_qc", width="stretch"):
                    st.session_state["prev_page"]   = "home"
                    st.session_state["qc_run_id"]   = run_id
                    st.session_state["active_page"] = "qc_detail"
                    st.rerun()
            with btn_c3:
                if zip_path and zip_path.exists():
                    with open(zip_path, "rb") as f:
                        st.download_button(
                            "💾 Télécharger",
                            data=f,
                            file_name=f"resultats_{run_id}.zip",
                            mime="application/zip",
                            key="dash_dl_last",
                            use_container_width=True,
                        )
                else:
                    st.button("💾 Télécharger", key="dash_dl_last_dis",
                              disabled=True, use_container_width=True)
        else:
            st.info("Aucun run enregistré. Lancez votre premier pipeline.")

        # ── Runs récents (mini-liste) ──────────
        st.markdown('<div class="dash-section">Runs récents</div>', unsafe_allow_html=True)
        recent = [h for h in history if h.get("status") in ("success", "failed")][:6]
        if recent:
            rows_html = ""
            for h in recent:
                s      = h.get("status", "")
                lbl    = "✓ Succès" if s == "success" else ("✗ Échec" if s == "failed" else "⟳")
                mcls   = "mini-ok" if s == "success" else ("mini-fail" if s == "failed" else "mini-run_")
                dur_s  = h.get("duration_sec")
                dur    = f"{dur_s//60}m {dur_s%60}s" if dur_s else "—"
                rows_html += f"""
                <div class="mini-run">
                    <div>
                        <div class="mini-run-id">{h['run_id']}</div>
                        <div class="mini-run-meta">
                            {h.get('date','—')} · {h.get('samples','?')} éch. · {dur}
                        </div>
                    </div>
                    <span class="mini-badge {mcls}">{lbl}</span>
                </div>"""
            st.markdown(rows_html, unsafe_allow_html=True)
        else:
            st.markdown(
                '<p style="font-size:0.83rem;color:#aaaaaa;">Aucun run disponible.</p>',
                unsafe_allow_html=True,
            )

    # ── Graphique d'activité ──────────────────
    with col_right:
        st.markdown('<div class="dash-section">Activité — 14 jours</div>', unsafe_allow_html=True)

        days  = _activity_last_n_days(history, 14)
        max_v = max((d["ok"] + d["fail"] for d in days), default=1) or 1
        MAX_H = 64  # hauteur max des barres en px

        bars_html = '<div class="activity-bar-wrap">'
        for d in days:
            total_d = d["ok"] + d["fail"]
            if total_d == 0:
                bars_html += (
                    '<div class="act-bar-col">'
                    f'<div class="act-bar act-bar-zero" style="height:4px;"></div>'
                    f'<div class="act-label">{d["label"]}</div>'
                    "</div>"
                )
            else:
                h_ok   = max(4, int(d["ok"]   / max_v * MAX_H))
                h_fail = max(4, int(d["fail"] / max_v * MAX_H)) if d["fail"] else 0
                bar_inner = ""
                if h_fail:
                    bar_inner += f'<div class="act-bar act-bar-fail" style="height:{h_fail}px;"></div>'
                if h_ok:
                    bar_inner += f'<div class="act-bar act-bar-ok" style="height:{h_ok}px;"></div>'
                bars_html += (
                    f'<div class="act-bar-col" title="{d["date"]}: {d["ok"]} succès, {d["fail"]} échec(s)">'
                    f"{bar_inner}"
                    f'<div class="act-label">{d["label"]}</div>'
                    "</div>"
                )
        bars_html += "</div>"

        # Légende
        bars_html += (
            '<div style="display:flex;gap:14px;margin-top:4px;">'
            '<span style="font-size:0.7rem;color:#28a745;">■ Succès</span>'
            '<span style="font-size:0.7rem;color:#dc3545;">■ Échec</span>'
            "</div>"
        )
        st.markdown(bars_html, unsafe_allow_html=True)

        # ── Stats rapides ──────────────────────
        st.markdown(
            '<div class="dash-section" style="margin-top:24px;">Répartition</div>',
            unsafe_allow_html=True,
        )
        if total > 0:
            pct_ok   = round(successes / total * 100)
            pct_fail = round(failures  / total * 100)
            pct_run  = 100 - pct_ok - pct_fail

            bar_html = (
                '<div style="height:10px;border-radius:6px;overflow:hidden;'
                'display:flex;background:#f0f0f0;margin-bottom:10px;">'
                f'<div style="width:{pct_ok}%;background:#28a745;"></div>'
                f'<div style="width:{pct_fail}%;background:#dc3545;"></div>'
                f'<div style="width:{pct_run}%;background:#f0ad00;"></div>'
                "</div>"
                '<div style="display:flex;gap:14px;font-size:0.75rem;color:#6b7280;">'
                f'<span><b style="color:#28a745">{pct_ok}%</b> succès</span>'
                f'<span><b style="color:#dc3545">{pct_fail}%</b> échec</span>'
                f'<span><b style="color:#f0ad00">{pct_run}%</b> en cours</span>'
                "</div>"
            )
            st.markdown(bar_html, unsafe_allow_html=True)
        else:
            st.markdown(
                '<p style="font-size:0.83rem;color:#aaaaaa;">Aucune donnée.</p>',
                unsafe_allow_html=True,
            )

    # ── Bouton Nouveau pipeline ────────────────
    st.markdown("<br>", unsafe_allow_html=True)
    if st.button("🔬 Lancer un nouveau pipeline", key="dash_to_pipeline", type="primary",
                 use_container_width=True):
        st.session_state["active_page"] = "pipeline"
        st.rerun()
