"""
qc_dashboard.py — Dashboard QC enrichi MaRS-py-upgrade
=======================================================
Remplace le bloc « elif active_page == "qc_detail": » dans app.py.

Intégration dans app.py :
  1. Importer en haut du fichier :
       from qc_dashboard import render_qc_dashboard
  2. Remplacer le bloc complet (de la ligne
         elif active_page == "qc_detail":
       jusqu'à la dernière ligne du bloc incluse) par :
         elif active_page == "qc_detail":
             render_qc_dashboard(
                 get_qc_data_fn   = get_qc_data_for_run,
                 get_run_paths_fn = get_run_paths,
                 get_badge_fn     = get_module_badge,
             )

Dépendances (déjà présentes dans ton environment.yml) :
  - pandas
  - plotly
"""

import base64
import io
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
import streamlit as st

# ─────────────────────────────────────────────
# CSS
# ─────────────────────────────────────────────
_QC_CSS = """
<style>
/* ── KPI cards QC ── */
.qc-dash-kpi-grid {
    display: grid;
    grid-template-columns: repeat(4, 1fr);
    gap: 12px;
    margin-bottom: 22px;
}
.qc-dash-kpi {
    background: #ffffff;
    border: 1px solid #e2e8f0;
    border-radius: 10px;
    padding: 14px 16px;
}
.qc-dash-kpi-label {
    font-size: 0.72rem;
    font-weight: 600;
    text-transform: uppercase;
    letter-spacing: 0.06em;
    color: #6b7280;
    margin-bottom: 6px;
}
.qc-dash-kpi-value {
    font-size: 1.6rem;
    font-weight: 700;
    color: #1a2a4a;
    line-height: 1.1;
}
.qc-dash-kpi-sub {
    font-size: 0.7rem;
    color: #aaaaaa;
    margin-top: 3px;
}

/* ── Statut badge global ── */
.qc-global-pass { color: #155724; background: #edfaf0; border: 1px solid #a8d5b5; }
.qc-global-warn { color: #7a4f00; background: #fff8e6; border: 1px solid #f0d080; }
.qc-global-fail { color: #721c24; background: #fdf0f0; border: 1px solid #f0b8b8; }
.qc-global-badge {
    display: inline-block;
    padding: 3px 10px;
    border-radius: 12px;
    font-size: 0.72rem;
    font-weight: 700;
}

/* ── Heatmap cell ── */
.hm-cell {
    display: inline-block;
    width: 18px; height: 18px;
    border-radius: 3px;
    margin: 1px;
    cursor: default;
}
.hm-pass { background: #28a745; }
.hm-warn { background: #f0ad00; }
.hm-fail { background: #dc3545; }
.hm-none { background: #e2e8f0; }
</style>
"""

# Modules FastQC dans l'ordre standard
_FASTQC_MODULE_ORDER = [
    "Basic Statistics",
    "Per base sequence quality",
    "Per sequence quality scores",
    "Per base sequence content",
    "Per sequence GC content",
    "Per base N content",
    "Sequence Length Distribution",
    "Sequence Duplication Levels",
    "Overrepresented sequences",
    "Adapter Content",
]


# ─────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────

def _global_status(mods: dict) -> str:
    vals = [v.lower() for v in mods.values()]
    if "fail" in vals:
        return "fail"
    if "warn" in vals:
        return "warn"
    if "pass" in vals:
        return "pass"
    return "—"


def _parse_gc(gc_str) -> float | None:
    try:
        return float(str(gc_str).replace("%", "").strip())
    except Exception:
        return None


def _parse_reads(reads_str) -> int | None:
    try:
        return int(str(reads_str).replace(",", "").replace(" ", "").strip())
    except Exception:
        return None


def _build_summary_df(qc_data: dict) -> pd.DataFrame:
    rows = []
    for name, d in qc_data.items():
        m    = d.get("txt_metrics", {})
        mods = m.get("modules", {})
        gs   = _global_status(mods)
        rows.append({
            "Échantillon":    name,
            "Lectures":       m.get("total_reads", "—"),
            "Long. lecture":  m.get("read_length", "—"),
            "%GC":            m.get("gc_pct", "—"),
            "Rapports HTML":  len(d.get("html_files", [])),
        })
    return pd.DataFrame(rows)


def _build_heatmap_data(qc_data: dict) -> tuple[list, list, list]:
    """Retourne (samples, modules, matrix) pour la heatmap."""
    samples = list(qc_data.keys())
    # collecter tous les modules présents
    all_mods: set = set()
    for d in qc_data.values():
        all_mods.update(d.get("txt_metrics", {}).get("modules", {}).keys())
    # ordre standard puis alphabétique pour les autres
    ordered = [m for m in _FASTQC_MODULE_ORDER if m in all_mods]
    ordered += sorted(m for m in all_mods if m not in ordered)

    status_to_num = {"pass": 1, "warn": 0.5, "fail": 0, "—": -1}
    matrix = []
    for samp in samples:
        mods = qc_data[samp].get("txt_metrics", {}).get("modules", {})
        row  = [status_to_num.get(mods.get(mod, "—").lower(), -1) for mod in ordered]
        matrix.append(row)
    return samples, ordered, matrix


def _coverage_df(paths: dict) -> pd.DataFrame | None:
    """
    Charge le fichier Summary/Reads_Metrics_Samples.csv si disponible,
    sinon tente de lire les fichiers _coverage.txt par échantillon.
    """
    summary_csv = paths["output"] / "Summary" / "Reads_Metrics_Samples.csv"
    if summary_csv.exists():
        try:
            return pd.read_csv(summary_csv)
        except Exception:
            pass

    # Fallback : lire samtools coverage par échantillon
    cov_dir = paths["output"] / "samtools_coverage"
    if not cov_dir.exists():
        return None
    rows = []
    for txt in sorted(cov_dir.glob("**/*_coverage.txt")):
        sample_id = txt.stem.replace("_coverage", "")
        try:
            df = pd.read_csv(txt, sep="\t")
            # colonnes samtools coverage : #rname, startpos, endpos, numreads,
            #   covbases, coverage, meandepth, meanbaseq, meanmapq
            if "coverage" in df.columns and "#rname" in df.columns:
                for _, r in df.iterrows():
                    rows.append({
                        "Sample":    sample_id,
                        "Gène":      r.get("#rname", "?"),
                        "Couverture (%)": round(float(r.get("coverage", 0)), 2),
                        "Profondeur moy.": round(float(r.get("meandepth", 0)), 1),
                    })
        except Exception:
            continue
    return pd.DataFrame(rows) if rows else None


# ─────────────────────────────────────────────
# Onglet 1 — Synthèse
# ─────────────────────────────────────────────

def _tab_synthese(qc_data: dict, run_id: str, get_badge_fn):
    df = _build_summary_df(qc_data)

    # KPI cards
    total    = len(df)

    gc_vals  = [_parse_gc(v) for v in df["%GC"] if _parse_gc(v) is not None]
    gc_mean  = f"{sum(gc_vals)/len(gc_vals):.1f}%" if gc_vals else "—"

    rd_vals  = [_parse_reads(v) for v in df["Lectures"] if _parse_reads(v) is not None]
    rd_total = f"{sum(rd_vals):,}" if rd_vals else "—"

    kpis = [
        ("Échantillons",  str(total), "total analysés"),
        ("%GC moyen",     gc_mean,    "tous échantillons"),
        ("Lectures tot.", rd_total,   "reads cumulés"),
    ]
    st.markdown('<div class="qc-dash-kpi-grid">', unsafe_allow_html=True)
    for label, val, sub in kpis:
        st.markdown(
            f'<div class="qc-dash-kpi">'
            f'<div class="qc-dash-kpi-label">{label}</div>'
            f'<div class="qc-dash-kpi-value">{val}</div>'
            f'<div class="qc-dash-kpi-sub">{sub}</div>'
            f'</div>',
            unsafe_allow_html=True,
        )
    st.markdown("</div>", unsafe_allow_html=True)

    # Tableau avec statut coloré
    st.markdown('<div class="section-label">Tableau récapitulatif</div>', unsafe_allow_html=True)

    st.dataframe(df, width="stretch", hide_index=True)

    # Export CSV
    csv_bytes = df.to_csv(index=False).encode("utf-8")
    st.download_button(
        "⬇️ Exporter en CSV",
        data=csv_bytes,
        file_name=f"qc_summary_{run_id}.csv",
        mime="text/csv",
        key="dl_qc_csv",
    )




# ─────────────────────────────────────────────
# Onglet 2 — Heatmap modules
# ─────────────────────────────────────────────

def _tab_heatmap(qc_data: dict):
    samples, modules, matrix = _build_heatmap_data(qc_data)
    if not samples or not modules:
        st.info("Données insuffisantes pour la heatmap.")
        return

    st.markdown('<div class="section-label">Heatmap FastQC — tous échantillons</div>', unsafe_allow_html=True)
    st.caption("🟢 PASS  🟡 WARN  🔴 FAIL  ⬜ Absent")

    colorscale = [
        [0.0,  "#dc3545"],   # fail
        [0.25, "#dc3545"],
        [0.26, "#f0ad00"],   # warn
        [0.74, "#f0ad00"],
        [0.75, "#28a745"],   # pass
        [1.0,  "#28a745"],
    ]

    # Texte hover
    status_map = {1: "PASS", 0.5: "WARN", 0: "FAIL", -1: "—"}
    hover = [
        [f"{samples[i]}<br>{modules[j]}<br>{status_map.get(matrix[i][j], '—')}"
         for j in range(len(modules))]
        for i in range(len(samples))
    ]

    fig = go.Figure(go.Heatmap(
        z=matrix,
        x=modules,
        y=samples,
        colorscale=colorscale,
        zmin=-1, zmax=1,
        showscale=False,
        hovertext=hover,
        hovertemplate="%{hovertext}<extra></extra>",
        xgap=2, ygap=2,
    ))
    h = max(300, len(samples) * 32 + 100)
    fig.update_layout(
        height=h,
        margin=dict(t=20, b=120, l=160, r=20),
        plot_bgcolor="white",
        xaxis=dict(tickangle=-40, tickfont=dict(size=11)),
        yaxis=dict(tickfont=dict(size=11)),
    )
    st.plotly_chart(fig, width="stretch")

    # Détail par module
    st.markdown('<div class="section-label" style="margin-top:10px;">Détail par module</div>', unsafe_allow_html=True)
    sel_mod = st.selectbox("Module FastQC", modules, key="hm_mod_sel")
    if sel_mod:
        mod_rows = []
        for samp in samples:
            mods = qc_data[samp].get("txt_metrics", {}).get("modules", {})
            mod_rows.append({"Échantillon": samp, "Statut": mods.get(sel_mod, "—").upper()})
        mod_df = pd.DataFrame(mod_rows)

        def color_s(val):
            return {"PASS": "background-color:#edfaf0", "WARN": "background-color:#fff8e6",
                    "FAIL": "background-color:#fdf0f0"}.get(val, "")

        try:
            st.dataframe(mod_df.style.map(color_s, subset=["Statut"]),
                         width="stretch", hide_index=True)
        except AttributeError:
            st.dataframe(mod_df.style.applymap(color_s, subset=["Statut"]),
                         width="stretch", hide_index=True)


# ─────────────────────────────────────────────
# Onglet 3 — Distribution %GC
# ─────────────────────────────────────────────

def _tab_gc(qc_data: dict):
    st.markdown('<div class="section-label">Distribution %GC par échantillon</div>', unsafe_allow_html=True)

    names, gc_vals = [], []
    for samp, d in qc_data.items():
        gc = _parse_gc(d.get("txt_metrics", {}).get("gc_pct", ""))
        if gc is not None:
            names.append(samp)
            gc_vals.append(gc)

    if not gc_vals:
        st.info("Aucune valeur %GC disponible.")
        return

    mean_gc = sum(gc_vals) / len(gc_vals)

    # Bar chart par échantillon
    colors = ["#dc3545" if v < 35 or v > 65 else ("#f0ad00" if v < 40 or v > 60 else "#1f70b8")
              for v in gc_vals]
    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=names, y=gc_vals,
        marker_color=colors,
        text=[f"{v}%" for v in gc_vals],
        textposition="outside",
        name="%GC",
    ))
    fig.add_hline(y=mean_gc, line_dash="dash", line_color="#888888",
                  annotation_text=f"Moy. {mean_gc:.1f}%",
                  annotation_position="bottom right")
    fig.add_hrect(y0=40, y1=60, fillcolor="#edfaf0", opacity=0.3, line_width=0,
                  annotation_text="Zone normale (40–60%)", annotation_position="top left")
    fig.update_layout(
        height=320,
        margin=dict(t=20, b=100, l=40, r=20),
        plot_bgcolor="white",
        yaxis=dict(range=[0, 100], title="%GC", gridcolor="#f0f0f0"),
        xaxis=dict(tickangle=-40),
        showlegend=False,
    )
    st.plotly_chart(fig, width="stretch")
    st.caption("🔵 Normal (40–60%)  🟡 Limite (35–40% / 60–65%)  🔴 Anormal")

    # Histogramme de distribution
    st.markdown('<div class="section-label" style="margin-top:10px;">Histogramme de distribution</div>', unsafe_allow_html=True)
    fig2 = go.Figure()
    fig2.add_trace(go.Histogram(
        x=gc_vals, nbinsx=20,
        marker_color="#1f70b8", opacity=0.75,
        name="%GC",
    ))
    fig2.add_vline(x=mean_gc, line_dash="dash", line_color="#dc3545",
                   annotation_text=f"Moy. {mean_gc:.1f}%")
    fig2.update_layout(
        height=240,
        margin=dict(t=10, b=40, l=40, r=20),
        plot_bgcolor="white",
        xaxis=dict(title="%GC", range=[0, 100]),
        yaxis=dict(title="Effectif", gridcolor="#f0f0f0"),
        showlegend=False,
    )
    st.plotly_chart(fig2, width="stretch")


# ─────────────────────────────────────────────
# Onglet 4 — Couverture
# ─────────────────────────────────────────────

def _tab_coverage(paths: dict):
    st.markdown('<div class="section-label">Couverture par gène MaRS</div>', unsafe_allow_html=True)

    df_cov = _coverage_df(paths)
    if df_cov is None or df_cov.empty:
        st.info("Aucun fichier de couverture trouvé (output/Summary/Reads_Metrics_Samples.csv ou output/samtools_coverage/).")
        return

    # Cas fichier Summary enrichi
    if "Sample" in df_cov.columns and "Gène" in df_cov.columns and "Couverture (%)" in df_cov.columns:
        samples_list = sorted(df_cov["Sample"].unique())
        sel = st.selectbox("Échantillon", ["Tous"] + samples_list, key="cov_samp_sel")
        df_plot = df_cov if sel == "Tous" else df_cov[df_cov["Sample"] == sel]

        fig = go.Figure()
        if sel == "Tous":
            for samp in samples_list:
                sub = df_plot[df_plot["Sample"] == samp]
                fig.add_trace(go.Bar(name=samp, x=sub["Gène"], y=sub["Couverture (%)"]))
            fig.update_layout(barmode="group")
        else:
            fig.add_trace(go.Bar(
                x=df_plot["Gène"], y=df_plot["Couverture (%)"],
                marker_color="#1f70b8",
                text=[f"{v}%" for v in df_plot["Couverture (%)"]],
                textposition="outside",
            ))
        fig.add_hline(y=80, line_dash="dash", line_color="#28a745",
                      annotation_text="Seuil 80%")
        fig.update_layout(
            height=340,
            margin=dict(t=20, b=80, l=50, r=20),
            plot_bgcolor="white",
            yaxis=dict(range=[0, 110], title="Couverture (%)", gridcolor="#f0f0f0"),
            xaxis=dict(title="Gène MaRS"),
        )
        st.plotly_chart(fig, width="stretch")

        # Profondeur moyenne
        if "Profondeur moy." in df_cov.columns:
            st.markdown('<div class="section-label" style="margin-top:10px;">Profondeur moyenne</div>', unsafe_allow_html=True)
            df_depth = df_plot.groupby("Gène")["Profondeur moy."].mean().reset_index()
            fig2 = go.Figure(go.Bar(
                x=df_depth["Gène"], y=df_depth["Profondeur moy."],
                marker_color="#6f42c1",
                text=[f"{v:.0f}x" for v in df_depth["Profondeur moy."]],
                textposition="outside",
            ))
            fig2.update_layout(
                height=280,
                margin=dict(t=10, b=60, l=50, r=20),
                plot_bgcolor="white",
                yaxis=dict(title="Profondeur moy. (x)", gridcolor="#f0f0f0"),
            )
            st.plotly_chart(fig2, width="stretch")

        # Export
        csv_bytes = df_cov.to_csv(index=False).encode("utf-8")
        st.download_button("⬇️ Exporter couverture CSV", data=csv_bytes,
                           file_name="coverage_summary.csv", mime="text/csv",
                           key="dl_cov_csv")
    else:
        # Affichage brut si colonnes différentes
        st.dataframe(df_cov, width="stretch")


# ─────────────────────────────────────────────
# Onglet 5 — Rapports HTML FastQC
# ─────────────────────────────────────────────

def _tab_html(qc_data: dict, run_id: str, get_badge_fn):
    samples = list(qc_data.keys())
    if not samples:
        st.info("Aucun échantillon disponible.")
        return

    selected = st.selectbox("Sélectionner un échantillon", samples, key="html_samp_sel")
    if not selected:
        return

    d          = qc_data[selected]
    m          = d.get("txt_metrics", {})
    mods       = m.get("modules", {})
    html_files = d.get("html_files", [])

    # KPIs de l'échantillon
    if m:
        c1, c2, c3 = st.columns(3)
        for col, label, value in [
            (c1, "Lectures totales", m.get("total_reads", "—")),
            (c2, "Longueur lecture", m.get("read_length", "—")),
            (c3, "% GC",            m.get("gc_pct", "—")),
        ]:
            col.markdown(
                f'<div class="qc-metric"><div class="qc-label">{label}</div>'
                f'<div class="qc-value">{value}</div></div>',
                unsafe_allow_html=True,
            )

    # Modules FastQC
    if mods:
        st.markdown('<div class="section-label" style="margin-top:16px;">Modules FastQC</div>', unsafe_allow_html=True)
        mod_cols = st.columns(2)
        for idx, (module, status) in enumerate(mods.items()):
            with mod_cols[idx % 2]:
                st.markdown(
                    f'<div style="display:flex;align-items:center;justify-content:space-between;'
                    f'padding:7px 12px;background:#ffffff;border:1px solid #d8d8d8;'
                    f'border-radius:4px;margin-bottom:4px;">'
                    f'<span style="font-size:0.82rem;color:#555;">{module}</span>'
                    f'{get_badge_fn(status)}</div>',
                    unsafe_allow_html=True,
                )

    # Rapports HTML
    st.markdown('<div class="section-label" style="margin-top:16px;">Rapports FastQC HTML</div>', unsafe_allow_html=True)
    if not html_files:
        st.info("Aucun fichier HTML FastQC trouvé pour cet échantillon.")
        return

    for html_path in sorted(html_files, key=lambda p: p.name):
        label  = html_path.name.replace("_fastqc.html", "")
        suffix = "R1" if "R1" in html_path.name.upper() else ("R2" if "R2" in html_path.name.upper() else "")
        with st.expander(f"📄 {label}{'  —  ' + suffix if suffix else ''}", expanded=(len(html_files) == 1)):
            try:
                html_content = html_path.read_text(errors="replace")
                st.download_button(
                    f"⬇️ Télécharger {html_path.name}",
                    data=html_content.encode("utf-8"),
                    file_name=html_path.name,
                    mime="text/html",
                    key=f"dl_html_{html_path.name}_{run_id}",
                )
                b64 = base64.b64encode(html_content.encode("utf-8")).decode()
                st.markdown(
                    f'<iframe src="data:text/html;base64,{b64}" width="100%" height="900px" '
                    f'style="border:1px solid #dee2e6;border-radius:8px;background:#fff;" '
                    f'sandbox="allow-scripts allow-same-origin"></iframe>',
                    unsafe_allow_html=True,
                )
            except Exception as e:
                st.error(f"Impossible de lire {html_path.name} : {e}")


# ─────────────────────────────────────────────
# Point d'entrée
# ─────────────────────────────────────────────

def render_qc_dashboard(get_qc_data_fn, get_run_paths_fn, get_badge_fn):
    """
    À appeler dans app.py :

        elif active_page == "qc_detail":
            render_qc_dashboard(
                get_qc_data_fn   = get_qc_data_for_run,
                get_run_paths_fn = get_run_paths,
                get_badge_fn     = get_module_badge,
            )
    """
    run_id = st.session_state.get("qc_run_id")
    if not run_id:
        st.session_state["active_page"] = "history"
        st.rerun()

    prev_page = st.session_state.get("prev_page", "home")
    prev_labels = {
        "home":     "🏠 Accueil",
        "pipeline": "🔬 Pipeline",
        "history":  "📋 Historique",
        "results":  "🧬 Résultats",
    }
    prev_label = prev_labels.get(prev_page, "← Retour")
    if st.button(f"← {prev_label}"):
        st.session_state["active_page"] = prev_page
        st.rerun()

    st.markdown(
        f'<div class="section-label">Dashboard QC — {run_id}</div>',
        unsafe_allow_html=True,
    )

    # Injecter le CSS à chaque rendu
    st.markdown(_QC_CSS, unsafe_allow_html=True)

    # Chargement des données
    @st.cache_data(ttl=120, show_spinner=False)
    def _load_qc(rid):
        return get_qc_data_fn(rid)

    qc_data = _load_qc(run_id)
    paths   = get_run_paths_fn(run_id)

    if not qc_data:
        st.info("Aucune donnée QC disponible pour ce run.")
        log_file = paths["logs"] / "pipeline.log"
        if log_file.exists():
            with st.expander("📄 Voir les logs"):
                st.code(log_file.read_text(errors="replace"), language="bash")
        return

    # Onglets
    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "📋 Synthèse",
        "🔥 Heatmap modules",
        "📊 Distribution %GC",
        "📈 Couverture",
        "📄 Rapports HTML",
    ])

    with tab1:
        _tab_synthese(qc_data, run_id, get_badge_fn)
    with tab2:
        _tab_heatmap(qc_data)
    with tab3:
        _tab_gc(qc_data)
    with tab4:
        _tab_coverage(paths)
    with tab5:
        _tab_html(qc_data, run_id, get_badge_fn)
