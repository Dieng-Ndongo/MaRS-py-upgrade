"""
results_page.py — Résultats moléculaires MaRS-py-upgrade
=========================================================
Corrections v3 :
  - Fix ValueError "truth value of DataFrame is ambiguous" (ligne df_summary = ...)
  - Haplotypes : Combined_Haplotypes.csv uniquement, individuel uniquement,
    haplotypes non valides (Null/nul*) exclus, fréquences + visualisation complète
  - Bouton export PDF intégré (via pdf_export.py)
  - GENE_DRUG avec vrais noms Pf du pipeline (Pfcrt, Pfdhfr, Pfdhps, Pfk13, Pfmdr1)
  - Gènes analysés cités explicitement avec médicament associé
"""

import re
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
import streamlit as st

# ─────────────────────────────────────────────
# Mapping gène → médicament(s)
# Noms exacts produits par le pipeline (colonne Gene de filtered_summary_clean.csv)
# ─────────────────────────────────────────────
GENE_DRUG = {
    "Pfk13":   ["Artémisinine"],
    "Pfcrt":   ["Chloroquine", "Amodiaquine"],
    "Pfmdr1":  ["Méfloquine", "Lumefantrine"],
    "Pfcytb":  ["Atovaquone"],
    "Pfdhps":  ["Sulfadoxine"],
    "Pfdhfr":  ["Pyriméthamine"],
    "Pfpfs47": ["Transmission"],
    "Pfmdr2":  ["Piperaquine"],
    # variantes minuscules / sans préfixe (robustesse)
    "pfk13":   ["Artémisinine"],
    "pfcrt":   ["Chloroquine", "Amodiaquine"],
    "pfmdr1":  ["Méfloquine", "Lumefantrine"],
    "pfcytb":  ["Atovaquone"],
    "pfdhps":  ["Sulfadoxine"],
    "pfdhfr":  ["Pyriméthamine"],
    "k13":     ["Artémisinine"],
    "crt":     ["Chloroquine", "Amodiaquine"],
    "mdr1":    ["Méfloquine", "Lumefantrine"],
    "cytb":    ["Atovaquone"],
    "dhps":    ["Sulfadoxine"],
    "dhfr":    ["Pyriméthamine"],
}

GENO_BG = {
    "WT":        ("#edfaf0", "#155724"),
    "WildType":  ("#edfaf0", "#155724"),
    "MT":        ("#fdf0f0", "#721c24"),
    "Mutant":    ("#fdf0f0", "#721c24"),
    "MIX":       ("#fff8e6", "#7a4f00"),
    "mutant":    ("#fdf0f0", "#721c24"),
    "wild-type": ("#edfaf0", "#155724"),
    "Null":      ("#f8f9fa", "#aaaaaa"),
}

# ─────────────────────────────────────────────
# CSS
# ─────────────────────────────────────────────
_CSS = """
<style>
.res-kpi-grid {
    display: grid;
    grid-template-columns: repeat(4, 1fr);
    gap: 12px;
    margin-bottom: 22px;
}
.res-kpi {
    background: #ffffff;
    border: 1px solid #e2e8f0;
    border-radius: 10px;
    padding: 14px 16px;
}
.res-kpi-label {
    font-size: 0.72rem; font-weight: 600;
    text-transform: uppercase; letter-spacing: 0.06em;
    color: #6b7280; margin-bottom: 6px;
}
.res-kpi-value { font-size: 1.6rem; font-weight: 700; color: #1a2a4a; line-height: 1.1; }
.res-kpi-sub   { font-size: 0.7rem; color: #aaaaaa; margin-top: 3px; }
.res-section {
    display: flex; align-items: center; gap: 10px;
    font-size: 0.78rem; font-weight: 700;
    text-transform: uppercase; letter-spacing: 0.08em;
    color: #888888; margin: 20px 0 12px;
}
.res-section::after { content:''; flex:1; height:1px; background:#d8d8d8; }
.haplo-freq-bar {
    height: 20px; border-radius: 4px; display: inline-block;
    margin-right: 6px; vertical-align: middle;
}
</style>
"""

# ─────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────

def _gene_clean(raw: str) -> str:
    s = str(raw)
    return s.split("|")[-1].strip().lower() if "|" in s else s.strip().lower()


def _drugs(gene: str) -> str:
    g = str(gene).strip()
    # correspondance exacte d'abord
    if g in GENE_DRUG:
        return ", ".join(GENE_DRUG[g])
    # fallback insensible à la casse
    g_lower = g.lower()
    for k, v in GENE_DRUG.items():
        if k.lower() in g_lower:
            return ", ".join(v)
    return "—"


def _vaf(val) -> float:
    try:
        return float(str(val).replace("%", "").replace(",", ".").strip())
    except Exception:
        return 0.0


def _geno_style(val: str):
    bg, tc = GENO_BG.get(str(val).strip(), ("#ffffff", "#333333"))
    return f"background-color:{bg};color:{tc};font-weight:700"


def _detect_geno(row: pd.Series, type_col: str, vartype_col: str) -> str:
    if type_col and type_col in row.index:
        v = str(row[type_col]).strip()
        if v in ("Mutant", "mutant", "missense_variant"):    return "MT"
        if v in ("WildType", "wild-type", "synonymous_variant"): return "WT"
    if vartype_col and vartype_col in row.index:
        v = str(row[vartype_col]).strip().upper()
        if v == "MT":    return "MT"
        if v == "WT":    return "WT"
        if "MIX" in v:   return "MIX"
    return "—"


def _is_individual(sample_name: str) -> bool:
    """Individuel si pas de P suivi de chiffres > 1 dans le nom AMD."""
    m = re.search(r'P(\d+)', str(sample_name))
    if m:
        return int(m.group(1)) <= 1
    return True


def _is_null_haplo(val: str) -> bool:
    """Haplotype non valide si commence par 'nul' (insensible à la casse)."""
    return str(val).strip().lower().startswith("nul") or str(val).strip() == ""


# ─────────────────────────────────────────────
# Loaders
# ─────────────────────────────────────────────

@st.cache_data(ttl=300, show_spinner=False)
def _load_csv(path: str) -> pd.DataFrame | None:
    p = Path(path)
    if not p.exists():
        return None
    try:
        return pd.read_csv(p, dtype=str).fillna("")
    except Exception:
        return None


def _load_all_svaf(sample_vaf_dir: Path) -> pd.DataFrame | None:
    merged = sample_vaf_dir.parent / "Sample_VAF_merge" / "Sample_VAF_merge.csv"
    if merged.exists():
        return _load_csv(str(merged))
    files = list(sample_vaf_dir.glob("*_SVAF.csv"))
    if not files:
        return None
    try:
        df = pd.concat([pd.read_csv(f, dtype=str) for f in files], ignore_index=True)
        return df.fillna("")
    except Exception:
        return None


# ─────────────────────────────────────────────
# Onglet 1 — Vue d'ensemble
# ─────────────────────────────────────────────

def _tab_overview(df: pd.DataFrame, run_id: str):
    id_col      = "ID"             if "ID"             in df.columns else None
    gene_col    = "Gene"           if "Gene"           in df.columns else \
                  "CHROM"          if "CHROM"          in df.columns else None
    vaf_col     = "Average_VAF(%)" if "Average_VAF(%)" in df.columns else \
                  "AVG_VAF"        if "AVG_VAF"        in df.columns else None
    type_col    = "Type"           if "Type"           in df.columns else \
                  "Annotation"     if "Annotation"     in df.columns else None
    vartype_col = "VARTYPE"        if "VARTYPE"        in df.columns else None

    n_samples  = df[id_col].nunique()         if id_col   else "—"
    genes_list = sorted(df[gene_col].unique()) if gene_col else []
    n_genes    = len(genes_list)
    n_snps     = len(df)

    mt = mix = wt = 0
    for _, row in df.iterrows():
        g = _detect_geno(row, type_col, vartype_col)
        if g == "MT":    mt  += 1
        elif g == "MIX": mix += 1
        elif g == "WT":  wt  += 1

    # KPI cards
    st.markdown('<div class="res-kpi-grid">', unsafe_allow_html=True)
    for label, val, sub in [
        ("Échantillons",   str(n_samples), "dans ce run"),
        ("Gènes analysés", str(n_genes),   "marqueurs moleculaires"),
        ("SNPs détectés",  str(n_snps),    "variants totaux"),
        ("Mutations MT",   str(mt),        f"dont MIX : {mix}"),
    ]:
        st.markdown(
            f'<div class="res-kpi">'
            f'<div class="res-kpi-label">{label}</div>'
            f'<div class="res-kpi-value">{val}</div>'
            f'<div class="res-kpi-sub">{sub}</div>'
            f'</div>',
            unsafe_allow_html=True,
        )
    st.markdown("</div>", unsafe_allow_html=True)

    # Gènes analysés — pills avec médicament
    if genes_list:
        st.markdown('<div class="res-section">Gènes analysés dans ce run</div>',
                    unsafe_allow_html=True)
        pills = ""
        for g in genes_list:
            drug = _drugs(g)
            pills += (
                f'<span style="display:inline-flex;flex-direction:column;align-items:center;'
                f'background:#f0f6ff;border:1px solid #cddcf5;border-radius:8px;'
                f'padding:8px 14px;margin:4px;min-width:100px;">'
                f'<span style="font-weight:700;color:#1a2a4a;font-size:0.85rem;">{g}</span>'
                f'<span style="font-size:0.7rem;color:#6b7280;margin-top:3px;">{drug}</span>'
                f'</span>'
            )
        st.markdown(
            f'<div style="display:flex;flex-wrap:wrap;gap:4px;margin-bottom:16px;">'
            f'{pills}</div>',
            unsafe_allow_html=True,
        )

    # Bar chart répartition
    st.markdown('<div class="res-section">Répartition des génotypes</div>',
                unsafe_allow_html=True)
    fig = go.Figure()
    for label, cnt, color in [("WT", wt, "#28a745"), ("MT", mt, "#dc3545"), ("MIX", mix, "#f0ad00")]:
        fig.add_trace(go.Bar(
            name=label, x=[label], y=[cnt],
            marker_color=color,
            text=[cnt], textposition="outside",
        ))
    fig.update_layout(
        height=240, margin=dict(t=10, b=10, l=30, r=10),
        plot_bgcolor="white", showlegend=False, bargap=0.45,
        yaxis=dict(gridcolor="#f0f0f0", title="Nombre de SNPs"),
    )
    st.plotly_chart(fig, width="stretch")

    # Tableau résumé par gène
    if gene_col:
        st.markdown('<div class="res-section">Résumé par gène</div>', unsafe_allow_html=True)
        rows = []
        for gene, grp in df.groupby(gene_col):
            mt_g = mix_g = wt_g = 0
            for _, r in grp.iterrows():
                g = _detect_geno(r, type_col, vartype_col)
                if g == "MT":    mt_g  += 1
                elif g == "MIX": mix_g += 1
                elif g == "WT":  wt_g  += 1
            vafs = [_vaf(v) for v in grp[vaf_col] if vaf_col and v != ""] if vaf_col else []
            rows.append({
                "Gène":       gene,
                "Médicament": _drugs(gene),
                "SNPs":       len(grp),
                "WT":         wt_g,
                "MT":         mt_g,
                "MIX":        mix_g,
                "VAF moy.":   f"{sum(vafs)/len(vafs):.1f}%" if vafs else "—",
            })
        df_g = pd.DataFrame(rows)

        def _row_color(row):
            if row["MT"] > 0:
                return ["background-color:#fff5f5"] * len(row)
            return [""] * len(row)

        st.dataframe(
            df_g.style.apply(_row_color, axis=1),
            width="stretch", hide_index=True,
        )


# ─────────────────────────────────────────────
# Onglet 2 — Tableau génotypes
# ─────────────────────────────────────────────

def _tab_genotypes(df: pd.DataFrame):
    id_col      = "ID"              if "ID"              in df.columns else \
                  "Sample_name"     if "Sample_name"     in df.columns else None
    gene_col    = "Gene"            if "Gene"            in df.columns else \
                  "CHROM"           if "CHROM"           in df.columns else None
    pos_col     = "POS"             if "POS"             in df.columns else None
    aa_col      = "Gene_Annotation" if "Gene_Annotation" in df.columns else \
                  "AA_change"       if "AA_change"       in df.columns else None
    vaf_col     = "Average_VAF(%)"  if "Average_VAF(%)"  in df.columns else \
                  "AVG_VAF"         if "AVG_VAF"         in df.columns else None
    type_col    = "Type"            if "Type"            in df.columns else \
                  "Annotation"      if "Annotation"      in df.columns else None
    vartype_col = "VARTYPE"         if "VARTYPE"         in df.columns else None
    cov_col     = "Average_Coverage" if "Average_Coverage" in df.columns else \
                  "AVG_COV"         if "AVG_COV"         in df.columns else None

    df = df.copy()
    df["Génotype"] = df.apply(lambda r: _detect_geno(r, type_col, vartype_col), axis=1)

    c1, c2, c3 = st.columns(3)
    genes = sorted(df[gene_col].unique()) if gene_col else []
    with c1:
        sel_gene = st.selectbox("Gène", ["Tous"] + genes, key="gt_gene")
    with c2:
        sel_geno = st.selectbox("Génotype", ["Tous", "MT", "MIX", "WT"], key="gt_geno")
    with c3:
        all_drugs = sorted({d for g in genes for d in _drugs(g).split(", ") if d != "—"})
        sel_drug  = st.selectbox("Médicament", ["Tous"] + all_drugs, key="gt_drug")

    df_f = df.copy()
    if sel_gene != "Tous" and gene_col:
        df_f = df_f[df_f[gene_col] == sel_gene]
    if sel_drug != "Tous" and gene_col:
        df_f = df_f[df_f[gene_col].apply(lambda g: sel_drug in _drugs(g))]
    if sel_geno != "Tous":
        df_f = df_f[df_f["Génotype"] == sel_geno]

    st.caption(f"{len(df_f)} variant(s) affiché(s)")

    if df_f.empty:
        st.info("Aucun variant ne correspond aux filtres sélectionnés.")
        return

    disp = {}
    if id_col:   disp["Échantillon"] = df_f[id_col]
    if gene_col: disp["Gène"]        = df_f[gene_col]
    if pos_col:  disp["Position"]    = df_f[pos_col]
    if aa_col:   disp["Mutation"]    = df_f[aa_col]
    disp["Génotype"] = df_f["Génotype"]
    if vaf_col:  disp["VAF (%)"]     = df_f[vaf_col].apply(
        lambda v: f"{_vaf(v):.2f}%" if v != "" else "—")
    if cov_col:  disp["Couv."]       = df_f[cov_col]
    if type_col: disp["Type"]        = df_f[type_col]

    df_disp = pd.DataFrame(disp)
    try:
        st.dataframe(
            df_disp.style.map(_geno_style, subset=["Génotype"]),
            width="stretch", hide_index=True,
        )
    except Exception:
        st.dataframe(df_disp, width="stretch", hide_index=True)

    csv = df_disp.to_csv(index=False).encode("utf-8")
    st.download_button("⬇️ Exporter CSV", data=csv,
                       file_name="genotypes.csv", mime="text/csv",
                       key="dl_geno")


# ─────────────────────────────────────────────
# Onglet 3 — Graphique VAF
# ─────────────────────────────────────────────

def _tab_vaf(df_svaf: pd.DataFrame | None, df_summary: pd.DataFrame | None):
    df = df_svaf if df_svaf is not None else df_summary
    if df is None:
        st.info("Aucun fichier VAF disponible.")
        return

    gene_col    = "CHROM"           if "CHROM"           in df.columns else \
                  "Gene"            if "Gene"            in df.columns else None
    aa_col      = "AA_change"       if "AA_change"       in df.columns else \
                  "Gene_Annotation" if "Gene_Annotation" in df.columns else None
    vaf_col     = "AVG_VAF"         if "AVG_VAF"         in df.columns else \
                  "Average_VAF(%)"  if "Average_VAF(%)"  in df.columns else None
    sample_col  = "Sample_name"     if "Sample_name"     in df.columns else \
                  "ID"              if "ID"              in df.columns else None
    type_col    = "Annotation"      if "Annotation"      in df.columns else \
                  "Type"            if "Type"            in df.columns else None
    vartype_col = "VARTYPE"         if "VARTYPE"         in df.columns else None

    if not vaf_col:
        st.info("Colonne VAF introuvable.")
        return

    genes = sorted(df[gene_col].unique()) if gene_col else []

    c1, c2 = st.columns(2)
    with c1:
        sel_gene = st.selectbox("Gène", genes or ["—"], key="vaf_gene")
    with c2:
        view = st.radio("Vue", ["Par mutation", "Par échantillon"],
                        horizontal=True, key="vaf_view")

    df_g = df.copy()
    if gene_col and sel_gene != "—":
        df_g = df_g[df_g[gene_col] == sel_gene]
    df_g["_vaf"]  = df_g[vaf_col].apply(_vaf)
    df_g["_geno"] = df_g.apply(lambda r: _detect_geno(r, type_col, vartype_col), axis=1)

    if df_g.empty:
        st.info("Aucune donnée pour ce gène.")
        return

    COLOR = {"MT": "#dc3545", "MIX": "#f0ad00", "WT": "#28a745", "—": "#1f70b8"}

    if view == "Par mutation":
        x_col = aa_col if aa_col else ("POS" if "POS" in df_g.columns else None)
        if not x_col:
            st.info("Colonne AA_change/POS introuvable.")
            return
        df_agg = df_g.groupby([x_col, "_geno"])["_vaf"].mean().reset_index()
        fig = go.Figure()
        for geno in ["WT", "MT", "MIX", "—"]:
            sub = df_agg[df_agg["_geno"] == geno]
            if sub.empty: continue
            fig.add_trace(go.Bar(
                name=geno, x=sub[x_col], y=sub["_vaf"],
                marker_color=COLOR[geno],
                text=[f"{v:.1f}%" for v in sub["_vaf"]],
                textposition="outside",
            ))
        fig.add_hline(y=50, line_dash="dash", line_color="#888",
                      annotation_text="50%", annotation_position="right")
        fig.update_layout(
            height=360, barmode="group",
            margin=dict(t=20, b=120, l=50, r=10),
            plot_bgcolor="white",
            yaxis=dict(title="AVG VAF (%)", range=[0, 115], gridcolor="#f0f0f0"),
            xaxis=dict(tickangle=-40),
            legend=dict(orientation="h", yanchor="bottom", y=1.02),
        )
        st.plotly_chart(fig, width="stretch")
    else:
        if not sample_col or not aa_col:
            st.info("Colonnes Sample_name ou AA_change manquantes.")
            return
        samples = sorted(df_g[sample_col].unique())
        sel_s = st.multiselect("Échantillons", samples, default=samples[:6], key="vaf_samp")
        df_s  = df_g[df_g[sample_col].isin(sel_s)] if sel_s else df_g
        fig = go.Figure()
        for samp in (sel_s or samples[:6]):
            sub = df_s[df_s[sample_col] == samp]
            fig.add_trace(go.Scatter(
                name=samp, x=sub[aa_col], y=sub["_vaf"],
                mode="markers+lines", marker=dict(size=8),
            ))
        fig.add_hline(y=50, line_dash="dash", line_color="#888", annotation_text="50%")
        fig.update_layout(
            height=360,
            margin=dict(t=20, b=120, l=50, r=10),
            plot_bgcolor="white",
            yaxis=dict(title="AVG VAF (%)", range=[0, 115], gridcolor="#f0f0f0"),
            xaxis=dict(title="Mutation AA", tickangle=-40),
        )
        st.plotly_chart(fig, width="stretch")

    st.caption("🟢 WT  🔴 MT  🟡 MIX  — Ligne pointillée : seuil 50%")


# ─────────────────────────────────────────────
# Onglet 4 — Haplotypes
# Individuel uniquement + haplotypes non valides exclus
# ─────────────────────────────────────────────

def _tab_haplotypes(paths: dict, run_id: str):
    # Source unique : Combined_Haplotypes.csv
    combined_path = paths["output"] / "haplotypes" / "Combined_Haplotypes.csv"
    df_hap = _load_csv(str(combined_path))

    if df_hap is None or df_hap.empty:
        st.info("Fichier Combined_Haplotypes.csv introuvable pour ce run.")
        return

    st.caption(f"Source : `Combined_Haplotypes.csv`")

    sample_col = next((c for c in ["Sample", "Sample_name", "SampleID"]
                       if c in df_hap.columns), None)
    if sample_col is None:
        st.error("Colonne Sample introuvable.")
        return

    haplo_cols = [c for c in df_hap.columns
                  if c not in {sample_col, "SITE", "YEAR", "COUNTRY", "SampleID"}]

    if not haplo_cols:
        st.dataframe(df_hap, width="stretch", hide_index=True)
        return

    # ── Filtre individuel ──────────────────────
    df_indiv = df_hap[df_hap[sample_col].apply(_is_individual)].copy()
    n_total  = len(df_hap)
    n_indiv  = len(df_indiv)
    n_pool   = n_total - n_indiv

    st.markdown(
        f'<div style="background:#f0f6ff;border:1px solid #cddcf5;border-radius:8px;'
        f'padding:10px 16px;margin-bottom:16px;font-size:0.83rem;">'
        f'📊 <b>{n_indiv}</b> échantillons individuels retenus sur {n_total} total '
        f'(<b>{n_pool}</b> poolés exclus du calcul)'
        f'</div>',
        unsafe_allow_html=True,
    )

    if df_indiv.empty:
        st.warning("Aucun échantillon individuel détecté.")
        return

    # ── Sélecteur gène ─────────────────────────
    sel_gene = st.selectbox("Gène / Combinaison", haplo_cols, key="hap_gene")

    # ── Calcul fréquences (individuel + valide) ─
    serie       = df_indiv[sel_gene].copy()
    mask_null   = serie.apply(_is_null_haplo)
    serie_valid = serie[~mask_null]
    n_valid     = len(serie_valid)
    n_null      = int(mask_null.sum())

    if n_null > 0:
        st.caption(f"⚠️ {n_null} haplotype(s) non valide(s) exclus (Null/nul…) — N valide = {n_valid}")

    if serie_valid.empty:
        st.info(f"Aucun haplotype valide pour {sel_gene}.")
        return

    vc = serie_valid.value_counts().reset_index()
    vc.columns = ["Haplotype", "N"]
    vc["% (individuel)"] = (vc["N"] / n_valid * 100).round(1).astype(str) + "%"

    # ── Layout : tableau + graphiques ──────────
    col_tbl, col_pie = st.columns([2, 1])

    with col_tbl:
        st.markdown('<div class="res-section">Fréquences</div>', unsafe_allow_html=True)

        def _style_haplo_row(row):
            h = str(row["Haplotype"]).strip().lower()
            if h == "wt" or h == "wildtype":
                return ["background-color:#edfaf0;color:#155724;font-weight:600"] * len(row)
            return ["background-color:#fff8e6;color:#7a4f00;font-weight:600"] * len(row)

        try:
            st.dataframe(
                vc.style.apply(_style_haplo_row, axis=1),
                width="stretch", hide_index=True,
            )
        except Exception:
            st.dataframe(vc, width="stretch", hide_index=True)

        # Export fréquences CSV
        csv_freq = vc.to_csv(index=False).encode("utf-8")
        st.download_button(
            "⬇️ Exporter fréquences CSV",
            data=csv_freq,
            file_name=f"haplotype_freq_{sel_gene}.csv",
            mime="text/csv",
            key="dl_haplo_freq",
        )

    with col_pie:
        st.markdown('<div class="res-section">Distribution</div>', unsafe_allow_html=True)
        PALETTE = ["#28a745", "#f0ad00", "#dc3545", "#1f70b8",
                   "#6f42c1", "#17a2b8", "#e83e8c", "#aaaaaa"]
        fig_pie = go.Figure(go.Pie(
            labels=vc["Haplotype"],
            values=vc["N"],
            hole=0.42,
            marker=dict(colors=PALETTE[:len(vc)]),
            textinfo="label+percent",
            textfont_size=11,
            hovertemplate="<b>%{label}</b><br>N = %{value}<br>%{percent}<extra></extra>",
        ))
        fig_pie.update_layout(
            height=300,
            margin=dict(t=10, b=10, l=10, r=10),
            showlegend=False,
        )
        st.plotly_chart(fig_pie, width="stretch")

    # ── Bar chart horizontal ────────────────────
    st.markdown('<div class="res-section">Effectifs par haplotype</div>',
                unsafe_allow_html=True)
    fig_bar = go.Figure()
    for i, row in vc.iterrows():
        color = "#28a745" if str(row["Haplotype"]).lower() in ("wt","wildtype") else \
                PALETTE[i % len(PALETTE)]
        fig_bar.add_trace(go.Bar(
            name=row["Haplotype"],
            x=[row["N"]],
            y=[row["Haplotype"]],
            orientation="h",
            marker_color=color,
            text=[f'N={row["N"]}  ({row["% (individuel)"]})'  ],
            textposition="outside",
            showlegend=False,
        ))
    fig_bar.update_layout(
        height=max(200, len(vc) * 45 + 60),
        margin=dict(t=10, b=20, l=160, r=80),
        plot_bgcolor="white",
        xaxis=dict(title=f"N individuel (total valide = {n_valid})",
                   gridcolor="#f0f0f0"),
        yaxis=dict(autorange="reversed"),
    )
    st.plotly_chart(fig_bar, width="stretch")

    # ── Tableau synthèse TOUS GÈNES ────────────
    st.markdown('<div class="res-section">Synthèse tous gènes — individuel uniquement</div>',
                unsafe_allow_html=True)

    all_gene_rows = []
    for gene in haplo_cols:
        serie_g       = df_indiv[gene].copy()
        mask_null_g   = serie_g.apply(_is_null_haplo)
        serie_valid_g = serie_g[~mask_null_g]
        n_valid_g     = len(serie_valid_g)
        n_null_g      = int(mask_null_g.sum())

        if serie_valid_g.empty:
            all_gene_rows.append({
                "Gène":          gene,
                "Haplotype":     "—",
                "N":             0,
                "% (individuel)":"—",
                "N valide":      0,
                "N exclus (Null)": n_null_g,
            })
            continue

        for haplo, count in serie_valid_g.value_counts().items():
            pct = f"{count / n_valid_g * 100:.1f}%"
            all_gene_rows.append({
                "Gène":            gene,
                "Haplotype":       haplo,
                "N":               int(count),
                "% (individuel)":  pct,
                "N valide":        n_valid_g,
                "N exclus (Null)": n_null_g,
            })

    df_all_genes = pd.DataFrame(all_gene_rows)

    def _style_all_genes(row):
        h = str(row["Haplotype"]).strip().lower()
        if h in ("wt", "wildtype"):
            return ["background-color:#edfaf0;color:#155724;font-weight:600"] * len(row)
        if h == "—":
            return ["color:#aaaaaa"] * len(row)
        return ["background-color:#fff8e6;color:#7a4f00;font-weight:600"] * len(row)

    try:
        st.dataframe(
            df_all_genes.style.apply(_style_all_genes, axis=1),
            width="stretch", hide_index=True,
        )
    except Exception:
        st.dataframe(df_all_genes, width="stretch", hide_index=True)

    # Export synthèse tous gènes
    csv_all_genes = df_all_genes.to_csv(index=False).encode("utf-8")
    st.download_button(
        "⬇️ Exporter synthèse tous gènes CSV",
        data=csv_all_genes,
        file_name=f"haplotype_freq_tous_genes_{run_id}.csv",
        mime="text/csv",
        key="dl_haplo_all_genes",
    )

    # ── Tableau complet ─────────────────────────
    with st.expander("📋 Tableau complet (tous échantillons)", expanded=False):
        df_show = df_hap[[sample_col] + haplo_cols].copy()
        df_show.insert(1, "Individuel", df_show[sample_col].apply(
            lambda s: "✓" if _is_individual(s) else "Pool"))
        st.dataframe(df_show, width="stretch", hide_index=True)

        csv_all = df_hap.to_csv(index=False).encode("utf-8")
        st.download_button(
            "⬇️ Exporter Combined_Haplotypes CSV complet",
            data=csv_all,
            file_name="Combined_Haplotypes_export.csv",
            mime="text/csv",
            key="dl_hap_all",
        )


# ─────────────────────────────────────────────
# Point d'entrée
# ─────────────────────────────────────────────

def render_results_page(get_run_paths_fn):
    # Injecter le CSS à chaque rendu (pas de cache sur le CSS)
    st.markdown(_CSS, unsafe_allow_html=True)

    # ── Sélecteur de run ──────────────────────
    runs_base = Path.home() / "MaRS-py-upgrade" / "runs"
    available = sorted(
        [d.name.replace("run_", "") for d in runs_base.glob("run_*") if d.is_dir()],
        reverse=True,
    ) if runs_base.exists() else []

    if not available:
        st.info("Aucun run disponible. Lancez d'abord un pipeline.")
        if st.button("🔬 Lancer un pipeline"):
            st.session_state["active_page"] = "pipeline"
            st.rerun()
        return

    default = st.session_state.get("results_run_id", available[0])
    idx     = available.index(default) if default in available else 0

    c_back, c_sel = st.columns([1, 3])
    with c_back:
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
    with c_sel:
        run_id = st.selectbox("Run", available, index=idx, key="res_run_sel")

    st.markdown(
        f'<div class="section-label">Résultats moléculaires — {run_id}</div>',
        unsafe_allow_html=True,
    )

    paths = get_run_paths_fn(run_id)
    out   = paths["output"]

    # ── Chargement des fichiers ───────────────
    summary_clean = out / "Summary_merge_filtered" / "filtered_summary_clean.csv"
    summary_merge = out / "Summary_merge_filtered" / "filtered_summary_merge.csv"
    svaf_dir      = out / "Sample_VAF"

    # FIX : pas de "or" sur DataFrame
    df_summary = _load_csv(str(summary_clean))
    if df_summary is None or df_summary.empty:
        df_summary = _load_csv(str(summary_merge))

    df_svaf = _load_all_svaf(svaf_dir)
    df_main = df_summary

    if df_main is None:
        is_running = (
            st.session_state.get("running", False)
            and st.session_state.get("run_id") == run_id
        )
        if is_running:
            current_step = st.session_state.get("current_step", 0)
            st.info(
                f"⏳ Pipeline en cours (étape {current_step}/37) — "
                "les résultats seront disponibles à partir de l'étape 22."
            )
            if st.button("🔬 Voir la progression", key="res_to_running"):
                st.session_state["active_page"] = "pipeline"
                st.rerun()
        else:
            st.warning(
                "Aucun fichier de résultats SNP trouvé pour ce run. "
                "Le pipeline s'est peut-être interrompu avant l'étape `filtered_summary_merge`."
            )
            st.code(f"Fichiers attendus :\n  {summary_clean}\n  {summary_merge}")
        return

    # ── Bouton export PDF ─────────────────────
    st.markdown("---")
    col_pdf1, col_pdf2 = st.columns([3, 1])
    with col_pdf1:
        st.markdown(
            '<div style="font-size:0.85rem;color:#6b7280;padding-top:6px;">'
            '📄 Exporter les résultats en PDF (SNPs + haplotypes individuels)'
            '</div>',
            unsafe_allow_html=True,
        )
    with col_pdf2:
        if st.button("📥 Générer le PDF", key="gen_pdf", width="stretch"):
            with st.spinner("Génération du PDF en cours..."):
                try:
                    from pdf_export import generate_pdf_report
                    combined_hap = out / "haplotypes" / "Combined_Haplotypes.csv"
                    df_hap_pdf   = None
                    if combined_hap.exists():
                        df_hap_pdf = pd.read_csv(combined_hap, dtype=str).fillna("")
                    else:
                        st.warning("Combined_Haplotypes.csv introuvable — PDF sans haplotypes.")
                    logo = Path.home() / "MaRS-py-upgrade" / "images" / "logoCIGASS.png"
                    pdf_bytes = generate_pdf_report(
                        run_id     = run_id,
                        output_dir = out,
                        df_summary = df_main,
                        df_hap     = df_hap_pdf,
                        logo_path  = logo if logo.exists() else None,
                    )
                    st.session_state["_pdf_bytes"]  = pdf_bytes
                    st.session_state["_pdf_run_id"] = run_id
                    st.toast("✅ PDF généré !", icon="📄")
                except Exception as e:
                    st.error(f"Erreur lors de la génération du PDF : {e}")

    if st.session_state.get("_pdf_bytes") and \
       st.session_state.get("_pdf_run_id") == run_id:
        st.download_button(
            label     = "⬇️ Télécharger le PDF",
            data      = st.session_state["_pdf_bytes"],
            file_name = f"resultats_moleculaires_{run_id}.pdf",
            mime      = "application/pdf",
            key       = "dl_pdf",
            width     = "stretch",
        )
    st.markdown("---")

    # ── Onglets ───────────────────────────────
    tab1, tab2, tab3, tab4 = st.tabs([
        "📊 Vue d'ensemble",
        "🧬 Tableau génotypes",
        "📈 Graphique VAF",
        "🔗 Haplotypes",
    ])

    with tab1:
        _tab_overview(df_main, run_id)
    with tab2:
        _tab_genotypes(df_main)
    with tab3:
        _tab_vaf(df_svaf, df_summary)
    with tab4:
        _tab_haplotypes(paths, run_id)
