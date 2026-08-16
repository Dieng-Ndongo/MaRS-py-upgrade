"""
pdf_export.py — Export PDF des résultats moléculaires MaRS-py-upgrade
======================================================================
Génère un PDF structuré avec :
  - Résumé global (N échantillons individuels, poolés, sites)
  - Tableau SNPs par gène (N, WT%, MT%, MIX%, VAF)
  - Tableau haplotypes sur ÉCHANTILLONS INDIVIDUELS UNIQUEMENT
    (individuel = get_pool_size(sample_name) == 1, i.e. pas de P+chiffres)
  - Graphiques en barres WT/MT/MIX par gène

Dépendances :
  pip install reportlab matplotlib
  (ou conda install reportlab matplotlib)

Intégration dans results_page.py :
  from pdf_export import generate_pdf_report
  Appeler dans _tab_overview() ou comme bouton séparé.
"""

import io
import re
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.platypus import (
    BaseDocTemplate, Frame, Image, PageTemplate,
    Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle
)

# ─────────────────────────────────────────────
# Constantes
# ─────────────────────────────────────────────
GENE_DRUG = {
    # Noms exacts produits par le pipeline
    "Pfk13":   "Artémisinine",
    "Pfcrt":   "Chloroquine / Amodiaquine",
    "Pfmdr1":  "Méfloquine / Lumefantrine",
    "Pfcytb":  "Atovaquone",
    "Pfdhps":  "Sulfadoxine",
    "Pfdhfr":  "Pyriméthamine",
    "Pfpfs47": "Transmission",
    "Pfmdr2":  "Piperaquine",
    # Variantes minuscules / sans préfixe (robustesse)
    "pfk13":   "Artémisinine",
    "pfcrt":   "Chloroquine / Amodiaquine",
    "pfmdr1":  "Méfloquine / Lumefantrine",
    "pfcytb":  "Atovaquone",
    "pfdhps":  "Sulfadoxine",
    "pfdhfr":  "Pyriméthamine",
    "k13":     "Artémisinine",
    "crt":     "Chloroquine / Amodiaquine",
    "mdr1":    "Méfloquine / Lumefantrine",
    "cytb":    "Atovaquone",
    "dhps":    "Sulfadoxine",
    "dhfr":    "Pyriméthamine",
}

HAPLO_GENES = ["DHFR", "DHPS", "CRT", "MDR"]  # gènes haplotypes standard MaRS

# ─────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────

def _is_individual(sample_name: str) -> bool:
    """
    Individuel si get_pool_size == 1.
    Logique AMD : présence de P suivi de chiffres → poolé.
    Ex: 18USGA00A1000  → individuel (pas de P+chiffres)
        18USGA00A000P05 → poolé (P05 = pool de 5)
    """
    match = re.search(r'P(\d+)', str(sample_name))
    if match:
        return int(match.group(1)) <= 1
    return True


def _gene_clean(raw: str) -> str:
    s = str(raw)
    return s.split("|")[-1].strip().lower() if "|" in s else s.strip().lower()


def _drug_for_gene(gene: str) -> str:
    g = str(gene).strip()
    # 1. Correspondance exacte
    if g in GENE_DRUG:
        return GENE_DRUG[g]
    # 2. Fallback insensible à la casse
    g_lower = g.lower()
    for k, v in GENE_DRUG.items():
        if k.lower() in g_lower:
            return v
    return "—"


def _pct(n, total):
    return f"{n/total*100:.1f}%" if total > 0 else "—"


def _bar_chart(labels, wt_vals, mt_vals, mix_vals, title: str) -> str:
    """Génère un bar chart WT/MT/MIX et retourne le chemin temporaire."""
    fig, ax = plt.subplots(figsize=(max(6, len(labels) * 0.9), 3.5))
    x = range(len(labels))
    w = 0.25
    ax.bar([i - w for i in x], wt_vals,  width=w, color="#28a745", label="WT")
    ax.bar([i     for i in x], mt_vals,  width=w, color="#dc3545", label="MT")
    ax.bar([i + w for i in x], mix_vals, width=w, color="#f0ad00", label="MIX")
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, rotation=40, ha="right", fontsize=8)
    ax.set_ylabel("N échantillons", fontsize=9)
    ax.set_title(title, fontsize=10, fontweight="bold")
    ax.legend(fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    tmp = tempfile.NamedTemporaryFile(suffix=".png", delete=False)
    fig.savefig(tmp.name, dpi=130)
    plt.close(fig)
    return tmp.name


def _haplo_pie(haplo_counts: dict, title: str) -> str:
    """Génère un camembert pour les haplotypes."""
    labels = list(haplo_counts.keys())
    values = list(haplo_counts.values())
    palette = ["#28a745", "#f0ad00", "#dc3545", "#1f70b8",
               "#6f42c1", "#aaaaaa", "#e83e8c", "#17a2b8"]
    fig, ax = plt.subplots(figsize=(4, 3.5))
    wedges, texts, autotexts = ax.pie(
        values, labels=None,
        colors=palette[:len(values)],
        autopct="%1.1f%%", pctdistance=0.75,
        startangle=90,
    )
    for t in autotexts:
        t.set_fontsize(8)
    ax.legend(wedges, labels, loc="lower center", bbox_to_anchor=(0.5, -0.25),
              ncol=2, fontsize=7)
    ax.set_title(title, fontsize=9, fontweight="bold")
    plt.tight_layout()
    tmp = tempfile.NamedTemporaryFile(suffix=".png", delete=False)
    fig.savefig(tmp.name, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return tmp.name


# ─────────────────────────────────────────────
# Styles PDF
# ─────────────────────────────────────────────

def _get_styles():
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(
        name="CellCenter", alignment=1, fontSize=8, leading=10))
    styles.add(ParagraphStyle(
        name="CellLeft", alignment=0, fontSize=8, leading=10))
    styles.add(ParagraphStyle(
        name="HeaderCell", alignment=1, fontSize=8, leading=10,
        textColor=colors.white, fontName="Helvetica-Bold"))
    styles.add(ParagraphStyle(
        name="SectionTitle", fontSize=11, fontName="Helvetica-Bold",
        textColor=colors.HexColor("#1a2a4a"), spaceAfter=6))
    styles.add(ParagraphStyle(
        name="GeneTitle", fontSize=10, fontName="Helvetica-Bold",
        textColor=colors.HexColor("#1f70b8"), spaceAfter=4))
    return styles


def _table_style_base():
    return TableStyle([
        ("BACKGROUND",   (0, 0), (-1, 0), colors.HexColor("#1a2a4a")),
        ("GRID",         (0, 0), (-1, -1), 0.4, colors.HexColor("#cccccc")),
        ("ALIGN",        (0, 0), (-1, -1), "CENTER"),
        ("VALIGN",       (0, 0), (-1, -1), "MIDDLE"),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1),
         [colors.white, colors.HexColor("#f4f8ff")]),
    ])


# ─────────────────────────────────────────────
# Sections PDF
# ─────────────────────────────────────────────

def _section_summary(story, df_summary: pd.DataFrame, styles):
    """Tableau résumé global : N échantillons, sites, individuel/poolé."""
    story.append(Paragraph("Résumé global", styles["SectionTitle"]))

    id_col   = next((c for c in ["ID", "Sample_name"] if c in df_summary.columns), None)
    gene_col = next((c for c in ["Gene", "CHROM"] if c in df_summary.columns), None)

    if id_col is None:
        story.append(Paragraph("Colonnes ID/Sample_name introuvables.", styles["Normal"]))
        story.append(Spacer(1, 10))
        return

    samples = df_summary[id_col].unique()
    n_total      = len(samples)
    n_individual = sum(1 for s in samples if _is_individual(s))
    n_pooled     = n_total - n_individual
    n_genes      = df_summary[gene_col].nunique() if gene_col else "—"
    n_snps       = len(df_summary)

    data = [
        [Paragraph(h, styles["HeaderCell"]) for h in
         ["Métrique", "Valeur"]],
        [Paragraph("Échantillons totaux",     styles["CellLeft"]),
         Paragraph(str(n_total),              styles["CellCenter"])],
        [Paragraph("Échantillons individuels",styles["CellLeft"]),
         Paragraph(str(n_individual),         styles["CellCenter"])],
        [Paragraph("Échantillons poolés",     styles["CellLeft"]),
         Paragraph(str(n_pooled),             styles["CellCenter"])],
        [Paragraph("Gènes analysés",          styles["CellLeft"]),
         Paragraph(str(n_genes),              styles["CellCenter"])],
        [Paragraph("SNPs détectés",           styles["CellLeft"]),
         Paragraph(str(n_snps),              styles["CellCenter"])],
    ]
    t = Table(data, colWidths=[10*cm, 5*cm])
    t.setStyle(_table_style_base())
    story.append(t)
    story.append(Spacer(1, 16))


def _section_snps(story, df_summary: pd.DataFrame, styles, img_paths: list):
    """Tableau SNPs par gène : N, WT, MT, MIX, VAF."""
    story.append(Paragraph("Analyse des SNPs par gène", styles["SectionTitle"]))

    id_col      = next((c for c in ["ID", "Sample_name"]       if c in df_summary.columns), None)
    gene_col    = next((c for c in ["Gene", "CHROM"]           if c in df_summary.columns), None)
    aa_col      = next((c for c in ["Gene_Annotation","AA_change"] if c in df_summary.columns), None)
    type_col    = next((c for c in ["Type", "Annotation"]      if c in df_summary.columns), None)
    vartype_col = "VARTYPE" if "VARTYPE" in df_summary.columns else None
    vaf_col     = next((c for c in ["Average_VAF(%)","AVG_VAF"] if c in df_summary.columns), None)

    if not gene_col or not id_col:
        story.append(Paragraph("Données SNP insuffisantes.", styles["Normal"]))
        return

    for gene, grp in df_summary.groupby(gene_col):
        gene_clean = _gene_clean(gene)
        drug       = _drug_for_gene(gene)

        story.append(Paragraph(
            f"{gene_clean.upper()}  —  {drug}", styles["GeneTitle"]))

        # Détecter génotypes
        def geno(row):
            if type_col:
                v = str(row.get(type_col, "")).strip()
                if v in ("Mutant", "mutant", "missense_variant"):   return "MT"
                if v in ("WildType", "wild-type", "synonymous_variant"): return "WT"
            if vartype_col:
                v = str(row.get(vartype_col, "")).strip().upper()
                if v == "MT":   return "MT"
                if v == "WT":   return "WT"
                if "MIX" in v:  return "MIX"
            return "—"

        grp = grp.copy()
        grp["_geno"] = grp.apply(geno, axis=1)

        headers = ["SNP / Mutation", "N total", "WT (N, %)", "MT (N, %)", "MIX (N, %)",
                   "VAF moy. (%)"]
        data = [[Paragraph(h, styles["HeaderCell"]) for h in headers]]

        snp_plot = []
        iter_col = aa_col if aa_col else gene_col

        for snp, sgrp in grp.groupby(iter_col):
            n_total = len(sgrp)
            n_wt    = (sgrp["_geno"] == "WT").sum()
            n_mt    = (sgrp["_geno"] == "MT").sum()
            n_mix   = (sgrp["_geno"] == "MIX").sum()
            vafs    = [float(str(v).replace("%","")) for v in sgrp[vaf_col]
                       if vaf_col and v not in ("", "—")] if vaf_col else []
            vaf_str = f"{sum(vafs)/len(vafs):.1f}%" if vafs else "—"

            data.append([
                Paragraph(str(snp),                              styles["CellLeft"]),
                Paragraph(str(n_total),                          styles["CellCenter"]),
                Paragraph(f"{n_wt} ({_pct(n_wt,  n_total)})",   styles["CellCenter"]),
                Paragraph(f"{n_mt} ({_pct(n_mt,  n_total)})",   styles["CellCenter"]),
                Paragraph(f"{n_mix} ({_pct(n_mix, n_total)})",  styles["CellCenter"]),
                Paragraph(vaf_str,                               styles["CellCenter"]),
            ])
            snp_plot.append((str(snp), n_wt, n_mt, n_mix))

        t = Table(data, colWidths=[5*cm, 2*cm, 3.2*cm, 3.2*cm, 3.2*cm, 2.5*cm],
                  repeatRows=1)
        ts = _table_style_base()
        # Colorer MT en rouge clair
        for i, row_data in enumerate(data[1:], start=1):
            mt_val = snp_plot[i-1][2] if i-1 < len(snp_plot) else 0
            if mt_val > 0:
                ts.add("BACKGROUND", (3, i), (3, i), colors.HexColor("#fff0f0"))
        t.setStyle(ts)
        story.append(t)
        story.append(Spacer(1, 6))

        # Graphique
        if snp_plot:
            labels   = [x[0] for x in snp_plot]
            wt_vals  = [x[1] for x in snp_plot]
            mt_vals  = [x[2] for x in snp_plot]
            mix_vals = [x[3] for x in snp_plot]
            img_path = _bar_chart(labels, wt_vals, mt_vals, mix_vals,
                                  f"{gene_clean.upper()} — WT / MT / MIX")
            img_paths.append(img_path)
            story.append(Image(img_path, width=16*cm, height=5.5*cm))
            story.append(Spacer(1, 10))


def _is_null_haplo(val: str) -> bool:
    """
    Retourne True si la valeur est un haplotype non valide :
    Null, Nul, NULL, NUL, null, nul, NullXxx, etc.
    (toute chaîne commençant par 'nul', insensible à la casse)
    """
    return str(val).strip().lower().startswith("nul")


def _section_haplotypes(story, df_hap: pd.DataFrame, styles, img_paths: list):
    """
    Tableau haplotypes depuis Combined_Haplotypes.csv
    sur ÉCHANTILLONS INDIVIDUELS UNIQUEMENT.
    Les haplotypes non valides (Null, NUL, nul*, etc.) sont exclus
    du dénominateur et des calculs.

    Colonnes attendues de Combined_Haplotypes.csv :
      Sample, DHFR and DHPS, CRT, MDR
    """
    story.append(Paragraph("Haplotypes — Échantillons individuels",
                            styles["SectionTitle"]))
    story.append(Paragraph(
        "Source : <b>Combined_Haplotypes.csv</b>. "
        "Les calculs sont restreints aux échantillons individuels (non poolés). "
        "Les haplotypes non valides (Null, Nul, NULL…) sont exclus des effectifs "
        "et pourcentages.",
        styles["Normal"]))
    story.append(Spacer(1, 8))

    sample_col = next((c for c in ["Sample", "Sample_name", "SampleID"]
                       if c in df_hap.columns), None)
    if sample_col is None:
        story.append(Paragraph("Colonne Sample introuvable dans Combined_Haplotypes.csv.",
                               styles["Normal"]))
        return

    # ── Filtrer individuel uniquement ─────────
    df_indiv = df_hap[df_hap[sample_col].apply(_is_individual)].copy()
    n_indiv  = len(df_indiv)
    n_total  = len(df_hap)

    story.append(Paragraph(
        f"Échantillons individuels retenus : <b>{n_indiv}</b> sur {n_total} au total.",
        styles["Normal"]))
    story.append(Spacer(1, 8))

    if df_indiv.empty:
        story.append(Paragraph("Aucun échantillon individuel détecté.", styles["Normal"]))
        return

    # Colonnes haplotypes = tout sauf Sample et métadonnées
    haplo_cols = [c for c in df_hap.columns
                  if c not in {sample_col, "SITE", "YEAR", "COUNTRY",
                                "SampleID", "Sample_ID"}]

    for gene in haplo_cols:
        story.append(Paragraph(gene, styles["GeneTitle"]))

        # ── Exclure les haplotypes non valides ──
        serie = df_indiv[gene].copy()
        serie_valid = serie[~serie.apply(_is_null_haplo)]
        n_null  = serie.apply(_is_null_haplo).sum()
        n_valid = len(serie_valid)

        if serie_valid.empty:
            story.append(Paragraph(
                f"Aucun haplotype valide pour {gene} "
                f"({n_null} non valide(s) exclu(s)).",
                styles["Normal"]))
            story.append(Spacer(1, 8))
            continue

        vc = serie_valid.value_counts()

        # Note sur exclusions
        if n_null > 0:
            story.append(Paragraph(
                f"<i>{n_null} haplotype(s) non valide(s) exclu(s) du calcul "
                f"(Null/NUL/nul…). N valide = {n_valid}.</i>",
                styles["Normal"]))
            story.append(Spacer(1, 4))

        headers = ["Haplotype", "N", "%"]
        data = [[Paragraph(h, styles["HeaderCell"]) for h in headers]]

        for haplo, count in vc.items():
            pct = f"{count / n_valid * 100:.1f}%"
            data.append([
                Paragraph(str(haplo), styles["CellLeft"]),
                Paragraph(str(count),  styles["CellCenter"]),
                Paragraph(pct,         styles["CellCenter"]),
            ])

        t = Table(data, colWidths=[8*cm, 3*cm, 3*cm], repeatRows=1)
        ts = _table_style_base()
        # Coloration des lignes
        for i, (haplo, _) in enumerate(vc.items(), start=1):
            h = str(haplo).strip().lower()
            if h == "wt" or h == "wildtype":
                ts.add("BACKGROUND", (0, i), (-1, i), colors.HexColor("#edfaf0"))
            else:
                ts.add("BACKGROUND", (0, i), (-1, i), colors.HexColor("#fff8e6"))
        t.setStyle(ts)
        story.append(t)
        story.append(Spacer(1, 6))

        # ── Camembert (haplotypes valides uniquement) ──
        haplo_dict = {str(k): int(v) for k, v in vc.items()}
        img_path = _haplo_pie(haplo_dict, f"Distribution haplotypes — {gene} (n={n_valid})")
        img_paths.append(img_path)
        story.append(Image(img_path, width=9*cm, height=7*cm))
        story.append(Spacer(1, 12))


# ─────────────────────────────────────────────
# Point d'entrée principal
# ─────────────────────────────────────────────

def generate_pdf_report(
    run_id: str,
    output_dir: Path,
    df_summary: pd.DataFrame,
    df_hap: pd.DataFrame | None = None,
    logo_path: Path | None = None,
) -> bytes:
    """
    Génère le PDF et retourne les bytes pour st.download_button.

    Paramètres
    ----------
    run_id      : identifiant du run (ex: "20260624_123417")
    output_dir  : dossier de sortie du pipeline (pour chercher les fichiers)
    df_summary  : filtered_summary_clean.csv ou filtered_summary_merge.csv
    df_hap      : Combined_Haplotypes.csv uniquement (colonnes: Sample, DHFR and DHPS, CRT, MDR)
    logo_path   : chemin vers le logo CIGASS (optionnel)

    Retourne
    --------
    bytes du PDF généré en mémoire
    """
    buf      = io.BytesIO()
    styles   = _get_styles()
    story    = []
    img_paths = []  # chemins temporaires à nettoyer

    doc = SimpleDocTemplate(
        buf,
        pagesize=landscape(A4),
        leftMargin=1.5*cm, rightMargin=1.5*cm,
        topMargin=1.5*cm,  bottomMargin=1.5*cm,
        title=f"MaRS-py-upgrade — Résultats moléculaires — {run_id}",
    )

    # ── En-tête ──────────────────────────────
    if logo_path and Path(logo_path).exists():
        story.append(Image(str(logo_path), width=6*cm, height=2*cm))
        story.append(Spacer(1, 8))

    story.append(Paragraph(
        "MaRS-py-upgrade — Résultats moléculaires",
        styles["Title"]))
    story.append(Paragraph(
        f"Run : <b>{run_id}</b> &nbsp;·&nbsp; "
        f"Généré le : <b>{pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}</b>",
        styles["Normal"]))
    story.append(Spacer(1, 14))

    # ── Méthodologie ─────────────────────────
    story.append(Paragraph("Méthodologie", styles["SectionTitle"]))
    story.append(Paragraph(
        "Ce rapport présente une synthèse de l'analyse moléculaire des marqueurs de résistance "
        "de <i>Plasmodium falciparum</i> à partir de données de séquençage FASTQ. "
        "Les fréquences alléliques variants (VAF) ont été calculées à partir des résultats "
        "du variant calling. Pour chaque mutation, une fréquence allélique moyenne (AVG_VAF) "
        "a été calculée. Les haplotypes sont calculés exclusivement sur les échantillons "
        "individuels (non poolés).",
        styles["Normal"]))
    story.append(Spacer(1, 14))

    # ── Résumé ───────────────────────────────
    _section_summary(story, df_summary, styles)

    # ── SNPs ─────────────────────────────────
    _section_snps(story, df_summary, styles, img_paths)

    # ── Haplotypes ───────────────────────────
    if df_hap is not None and not df_hap.empty:
        story.append(Spacer(1, 10))
        _section_haplotypes(story, df_hap, styles, img_paths)

    # ── Build PDF ─────────────────────────────
    doc.build(story)
    pdf_bytes = buf.getvalue()

    # Nettoyage des images temporaires
    import os
    for p in img_paths:
        try:
            os.unlink(p)
        except Exception:
            pass

    return pdf_bytes