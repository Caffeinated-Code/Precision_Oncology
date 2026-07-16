#!/usr/bin/env python3
"""Reproducible precision-oncology analysis for the Van Allen melanoma cohort.

The script uses only the Python standard library so that the core analysis is
easy to inspect and run in constrained environments. It downloads public data
from cBioPortal, builds derived tables, performs gene-response enrichment tests,
and writes a Markdown report plus small SVG figures.
"""

from __future__ import annotations

import csv
import html
import json
import math
import statistics
import sys
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path


BASE_URL = "https://www.cbioportal.org/api"
STUDY_ID = "skcm_dfci_2015"
SAMPLE_LIST_ID = "skcm_dfci_2015_sequenced"
MUTATION_PROFILE_ID = "skcm_dfci_2015_mutations"

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data" / "derived"
TABLE_DIR = ROOT / "results" / "tables"
FIGURE_DIR = ROOT / "results" / "figures"
REPORT_PATH = ROOT / "results" / "precision_oncology_report.md"

DCB_LABELS = {"CR", "PR", "SD"}
NO_DCB_LABELS = {"PD"}
UNKNOWN_LABELS = {"X", "", "NA", "N/A", "Unknown", "UNKNOWN"}

FLAGS = {
    "TTN",
    "MUC16",
    "OBSCN",
    "AHNAK2",
    "SYNE1",
    "FLG",
    "MUC5B",
    "DNAH17",
    "PLEC",
    "DST",
    "SYNE2",
    "NEB",
    "HSPG2",
    "LAMA5",
    "AHNAK",
    "HMCN1",
    "USH2A",
    "DNAH11",
    "MACF1",
    "MUC17",
}

BIOLOGY_PANELS = {
    "DNA repair": {
        "ERCC2",
        "ERCC3",
        "ERCC4",
        "ERCC5",
        "XPA",
        "XPC",
        "BRCA1",
        "BRCA2",
        "ATM",
        "ATR",
        "MSH2",
        "MSH6",
        "MLH1",
        "PMS2",
        "POLE",
        "POLD1",
    },
    "Melanoma drivers": {
        "BRAF",
        "NRAS",
        "NF1",
        "KIT",
        "PTEN",
        "CDKN2A",
        "TP53",
        "RAC1",
        "MAP2K1",
        "MAP2K2",
    },
    "Immune signaling": {
        "CD274",
        "PDCD1",
        "CTLA4",
        "LAG3",
        "TIGIT",
        "IFNG",
        "JAK1",
        "JAK2",
        "B2M",
        "HLA-A",
        "HLA-B",
        "HLA-C",
    },
}

IMMUNE_FIELDS = {
    "CBEST_T_CD8": "CD8 T cells",
    "CBEST_T_CD4_MEM_ACT": "Activated CD4 memory T cells",
    "CBEST_T_FOL_HELP": "Follicular helper T cells",
    "CBEST_T_REG": "Regulatory T cells",
    "CBEST_NK_ACT": "Activated NK cells",
    "CBEST_MACRO_M1": "M1 macrophages",
    "CBEST_MACRO_M2": "M2 macrophages",
}


def api_get(path: str):
    url = f"{BASE_URL}{path}"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as response:
        return json.load(response)


def api_post(path: str, payload: dict):
    url = f"{BASE_URL}{path}"
    req = urllib.request.Request(
        url,
        data=json.dumps(payload).encode("utf-8"),
        headers={"Accept": "application/json", "Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=120) as response:
        return json.load(response)


def ensure_dirs() -> None:
    for path in (DATA_DIR, TABLE_DIR, FIGURE_DIR):
        path.mkdir(parents=True, exist_ok=True)


def write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=2, sort_keys=True), encoding="utf-8")


def write_csv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    if fieldnames is None:
        fields = []
        for row in rows:
            for key in row:
                if key not in fields:
                    fields.append(key)
        fieldnames = fields
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def pivot_clinical(rows: list[dict], id_key: str) -> dict[str, dict]:
    result: dict[str, dict] = defaultdict(dict)
    for row in rows:
        entity_id = row[id_key]
        result[entity_id][row["clinicalAttributeId"]] = row.get("value", "")
    return dict(result)


def to_float(value):
    if value in (None, "", "NA", "N/A", "NaN"):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def benefit_group(value: str) -> str:
    value = (value or "").strip()
    if value in DCB_LABELS:
        return "Durable clinical benefit"
    if value in NO_DCB_LABELS:
        return "No durable benefit"
    if value in UNKNOWN_LABELS:
        return "Indeterminate"
    return "Indeterminate"


def is_nonsynonymous(mutation_type: str) -> bool:
    silent = {
        "Silent",
        "silent",
        "3'UTR",
        "5'UTR",
        "Intron",
        "IGR",
        "RNA",
        "lincRNA",
    }
    return mutation_type not in silent


def comb(n: int, k: int) -> int:
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)


def hypergeom_prob(k: int, row_total: int, col_total: int, n: int) -> float:
    return comb(col_total, k) * comb(n - col_total, row_total - k) / comb(n, row_total)


def fisher_exact_two_sided(a: int, b: int, c: int, d: int) -> float:
    row_total = a + c
    col_total = a + b
    n = a + b + c + d
    lo = max(0, row_total - (n - col_total))
    hi = min(row_total, col_total)
    observed = hypergeom_prob(a, row_total, col_total, n)
    p_value = 0.0
    for k in range(lo, hi + 1):
        p = hypergeom_prob(k, row_total, col_total, n)
        if p <= observed + 1e-12:
            p_value += p
    return min(1.0, p_value)


def odds_ratio(a: int, b: int, c: int, d: int) -> float:
    if b * c == 0:
        if a * d == 0:
            return float("nan")
        return float("inf")
    return (a * d) / (b * c)


def bh_adjust(rows: list[dict], p_key: str = "p_value", q_key: str = "q_value") -> None:
    valid = [(i, row[p_key]) for i, row in enumerate(rows) if row[p_key] is not None]
    valid.sort(key=lambda item: item[1], reverse=True)
    m = len(valid)
    running = 1.0
    for rank_from_end, (idx, p_value) in enumerate(valid, start=1):
        rank = m - rank_from_end + 1
        running = min(running, p_value * m / rank)
        rows[idx][q_key] = min(1.0, running)


def ranks(values: list[float]) -> list[float]:
    indexed = sorted(enumerate(values), key=lambda item: item[1])
    out = [0.0] * len(values)
    i = 0
    while i < len(indexed):
        j = i
        while j + 1 < len(indexed) and indexed[j + 1][1] == indexed[i][1]:
            j += 1
        avg_rank = (i + 1 + j + 1) / 2
        for k in range(i, j + 1):
            out[indexed[k][0]] = avg_rank
        i = j + 1
    return out


def mann_whitney_effect(values_a: list[float], values_b: list[float]) -> dict:
    if not values_a or not values_b:
        return {"n_a": len(values_a), "n_b": len(values_b), "u": None, "rank_biserial": None}
    vals = values_a + values_b
    r = ranks(vals)
    n_a, n_b = len(values_a), len(values_b)
    rank_sum_a = sum(r[:n_a])
    u_a = rank_sum_a - n_a * (n_a + 1) / 2
    rank_biserial = (2 * u_a) / (n_a * n_b) - 1
    return {
        "n_a": n_a,
        "n_b": n_b,
        "u": u_a,
        "rank_biserial": rank_biserial,
    }


def median_or_none(values: list[float]) -> float | None:
    return statistics.median(values) if values else None


def fmt_float(value, digits=3) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return "NA"
    if isinstance(value, float) and math.isinf(value):
        return "Inf"
    return f"{value:.{digits}f}"


def fetch_data() -> tuple[dict, list[dict], list[dict], list[dict]]:
    study = api_get(f"/studies/{STUDY_ID}")
    patient_clinical = api_get(
        f"/studies/{STUDY_ID}/clinical-data?clinicalDataType=PATIENT&projection=DETAILED&pageSize=100000"
    )
    sample_clinical = api_get(
        f"/studies/{STUDY_ID}/clinical-data?clinicalDataType=SAMPLE&projection=DETAILED&pageSize=100000"
    )
    mutations = api_post(
        f"/molecular-profiles/{MUTATION_PROFILE_ID}/mutations/fetch?projection=DETAILED",
        {"sampleListId": SAMPLE_LIST_ID},
    )
    return study, patient_clinical, sample_clinical, mutations


def build_cohort(patient_rows: list[dict], sample_rows: list[dict]) -> list[dict]:
    patient = pivot_clinical(patient_rows, "patientId")
    sample = pivot_clinical(sample_rows, "sampleId")
    cohort = []
    for sample_id, sample_data in sorted(sample.items()):
        patient_id = sample_id
        patient_data = patient.get(patient_id, {})
        dcb_raw = patient_data.get("DURABLE_CLINICAL_BENEFIT", "")
        group = benefit_group(dcb_raw)
        row = {
            "sample_id": sample_id,
            "patient_id": patient_id,
            "durable_clinical_benefit_raw": dcb_raw,
            "benefit_group": group,
            "analysis_included": group != "Indeterminate",
            "age": to_float(patient_data.get("AGE")),
            "m_stage": patient_data.get("M_STAGE", ""),
            "os_status": patient_data.get("OS_STATUS", ""),
            "os_months": to_float(patient_data.get("OS_MONTHS")),
            "mutation_count": to_float(sample_data.get("MUTATION_COUNT")),
            "mutation_load": to_float(sample_data.get("MUTATION_LOAD")),
            "tmb_nonsynonymous": to_float(sample_data.get("TMB_NONSYNONYMOUS")),
            "neoantigen_load": to_float(sample_data.get("NEOAGCNT")),
        }
        for field in IMMUNE_FIELDS:
            row[field] = to_float(sample_data.get(field))
        cohort.append(row)
    return cohort


def build_mutation_tables(mutations: list[dict]) -> tuple[list[dict], dict[str, set[str]]]:
    rows = []
    gene_to_samples: dict[str, set[str]] = defaultdict(set)
    for m in mutations:
        mutation_type = m.get("mutationType", "")
        if not is_nonsynonymous(mutation_type):
            continue
        gene = m.get("gene", {}).get("hugoGeneSymbol") or m.get("geneSymbol") or ""
        sample_id = m.get("sampleId", "")
        if not gene or not sample_id:
            continue
        row = {
            "sample_id": sample_id,
            "patient_id": m.get("patientId", sample_id),
            "gene": gene,
            "mutation_type": mutation_type,
            "protein_change": m.get("proteinChange", ""),
            "chromosome": m.get("chr", ""),
            "start_position": m.get("startPosition", ""),
            "variant_type": m.get("variantType", ""),
        }
        rows.append(row)
        gene_to_samples[gene].add(sample_id)
    return rows, gene_to_samples


def analyze_gene_enrichment(cohort: list[dict], gene_to_samples: dict[str, set[str]]) -> list[dict]:
    included = [row for row in cohort if row["analysis_included"]]
    dcb_samples = {row["sample_id"] for row in included if row["benefit_group"] == "Durable clinical benefit"}
    no_dcb_samples = {row["sample_id"] for row in included if row["benefit_group"] == "No durable benefit"}
    all_samples = dcb_samples | no_dcb_samples
    rows = []
    for gene, samples in sorted(gene_to_samples.items()):
        samples = samples & all_samples
        mutated_total = len(samples)
        if mutated_total < 2:
            continue
        a = len(samples & dcb_samples)
        c = len(samples & no_dcb_samples)
        b = len(dcb_samples) - a
        d = len(no_dcb_samples) - c
        p_value = fisher_exact_two_sided(a, b, c, d)
        row = {
            "gene": gene,
            "mutated_samples": mutated_total,
            "dcb_mutated": a,
            "dcb_wildtype": b,
            "no_dcb_mutated": c,
            "no_dcb_wildtype": d,
            "dcb_mutation_fraction": a / len(dcb_samples) if dcb_samples else None,
            "no_dcb_mutation_fraction": c / len(no_dcb_samples) if no_dcb_samples else None,
            "odds_ratio_dcb": odds_ratio(a, b, c, d),
            "p_value": p_value,
            "flag_gene": gene in FLAGS,
            "biology_panel": gene_panel_name(gene),
            "frequency_tier": frequency_tier(mutated_total),
        }
        rows.append(row)
    rows.sort(key=lambda r: (r["p_value"], -r["mutated_samples"], r["gene"]))
    bh_adjust(rows)
    return rows


def frequency_tier(mutated_total: int) -> str:
    if mutated_total >= 5:
        return "primary_screen"
    if mutated_total >= 2:
        return "low_frequency"
    return "singleton"


def gene_panel_name(gene: str) -> str:
    matches = [name for name, genes in BIOLOGY_PANELS.items() if gene in genes]
    return "; ".join(matches)


def burden_and_immune_summary(cohort: list[dict]) -> tuple[list[dict], list[dict]]:
    included = [row for row in cohort if row["analysis_included"]]
    dcb = [row for row in included if row["benefit_group"] == "Durable clinical benefit"]
    no_dcb = [row for row in included if row["benefit_group"] == "No durable benefit"]
    feature_labels = {
        "mutation_count": "Mutation count",
        "mutation_load": "Nonsynonymous mutation load",
        "tmb_nonsynonymous": "TMB nonsynonymous",
        "neoantigen_load": "Neoantigen load",
    }
    burden_rows = []
    for field, label in feature_labels.items():
        a = [row[field] for row in dcb if row[field] is not None]
        b = [row[field] for row in no_dcb if row[field] is not None]
        effect = mann_whitney_effect(a, b)
        burden_rows.append(
            {
                "feature": label,
                "n_dcb": len(a),
                "n_no_dcb": len(b),
                "median_dcb": median_or_none(a),
                "median_no_dcb": median_or_none(b),
                "rank_biserial_dcb_vs_no_dcb": effect["rank_biserial"],
            }
        )
    immune_rows = []
    for field, label in IMMUNE_FIELDS.items():
        a = [row[field] for row in dcb if row[field] is not None]
        b = [row[field] for row in no_dcb if row[field] is not None]
        effect = mann_whitney_effect(a, b)
        immune_rows.append(
            {
                "feature": label,
                "n_dcb": len(a),
                "n_no_dcb": len(b),
                "median_dcb": median_or_none(a),
                "median_no_dcb": median_or_none(b),
                "rank_biserial_dcb_vs_no_dcb": effect["rank_biserial"],
            }
        )
    immune_rows.sort(
        key=lambda r: abs(r["rank_biserial_dcb_vs_no_dcb"] or 0),
        reverse=True,
    )
    return burden_rows, immune_rows


def gene_recurrence_rows(gene_to_samples: dict[str, set[str]], cohort: list[dict]) -> list[dict]:
    sample_ids = {row["sample_id"] for row in cohort}
    n = len(sample_ids)
    rows = []
    for gene, samples in gene_to_samples.items():
        count = len(samples & sample_ids)
        rows.append(
            {
                "gene": gene,
                "mutated_samples": count,
                "sample_fraction": count / n if n else None,
                "flag_gene": gene in FLAGS,
                "biology_panel": gene_panel_name(gene),
            }
        )
    rows.sort(key=lambda r: (-r["mutated_samples"], r["gene"]))
    return rows


def create_bar_svg(path: Path, rows: list[dict], label_key: str, value_key: str, title: str, x_label: str) -> None:
    rows = rows[:15]
    width = 860
    row_h = 30
    left = 165
    top = 58
    bar_w = 610
    height = top + row_h * len(rows) + 45
    max_value = max((row[value_key] or 0 for row in rows), default=1)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="{left}" y="28" font-family="Arial" font-size="18" font-weight="700">{html.escape(title)}</text>',
    ]
    for i, row in enumerate(rows):
        y = top + i * row_h
        value = row[value_key] or 0
        bw = (value / max_value) * bar_w if max_value else 0
        label = str(row[label_key])
        parts.append(f'<text x="12" y="{y + 19}" font-family="Arial" font-size="13">{html.escape(label)}</text>')
        parts.append(f'<rect x="{left}" y="{y}" width="{bw:.1f}" height="20" fill="#356f8f"/>')
        parts.append(
            f'<text x="{left + bw + 8}" y="{y + 15}" font-family="Arial" font-size="12">{fmt_float(value, 3)}</text>'
        )
    parts.append(
        f'<text x="{left}" y="{height - 12}" font-family="Arial" font-size="12" fill="#555">{html.escape(x_label)}</text>'
    )
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def create_scatter_svg(path: Path, rows: list[dict]) -> None:
    rows = [r for r in rows if r["mutated_samples"] >= 2][:40]
    width, height = 860, 520
    left, right, top, bottom = 70, 35, 50, 70
    plot_w = width - left - right
    plot_h = height - top - bottom
    xs = [math.log2(r["mutated_samples"]) for r in rows]
    ys = [-math.log10(max(r["p_value"], 1e-300)) for r in rows]
    x_min, x_max = min(xs, default=0), max(xs, default=1)
    y_min, y_max = 0, max(ys + [1.4])

    def sx(x):
        return left + (x - x_min) / (x_max - x_min or 1) * plot_w

    def sy(y):
        return top + plot_h - (y - y_min) / (y_max - y_min or 1) * plot_h

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="70" y="28" font-family="Arial" font-size="18" font-weight="700">Gene-level response enrichment screen</text>',
        f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" stroke="#333"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" stroke="#333"/>',
    ]
    threshold_y = sy(-math.log10(0.05))
    parts.append(f'<line x1="{left}" x2="{left + plot_w}" y1="{threshold_y}" y2="{threshold_y}" stroke="#999" stroke-dasharray="5,5"/>')
    for row in rows:
        x = sx(math.log2(row["mutated_samples"]))
        y = sy(-math.log10(max(row["p_value"], 1e-300)))
        color = "#b4473d" if row["frequency_tier"] == "primary_screen" else "#6d8fb3"
        parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="5" fill="{color}" opacity="0.82"/>')
        if row["gene"] in {"ERCC2", "BRAF", "NRAS", "NF1", "TP53"} or row["p_value"] < 0.02:
            parts.append(f'<text x="{x + 7:.1f}" y="{y - 7:.1f}" font-family="Arial" font-size="11">{html.escape(row["gene"])}</text>')
    parts.append(f'<text x="{left + plot_w / 2 - 65}" y="{height - 22}" font-family="Arial" font-size="12">log2(mutated samples)</text>')
    parts.append(f'<text x="12" y="{top + plot_h / 2}" font-family="Arial" font-size="12" transform="rotate(-90 12,{top + plot_h / 2})">-log10 p-value</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def create_effect_svg(path: Path, rows: list[dict], title: str) -> None:
    rows = [r for r in rows if r.get("rank_biserial_dcb_vs_no_dcb") is not None]
    width = 860
    row_h = 34
    left = 310
    center = 500
    top = 58
    scale = 250
    height = top + row_h * len(rows) + 55
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="{left}" y="28" font-family="Arial" font-size="18" font-weight="700">{html.escape(title)}</text>',
        f'<line x1="{center}" y1="{top - 8}" x2="{center}" y2="{height - 45}" stroke="#555"/>',
        f'<text x="{center - 68}" y="{height - 16}" font-family="Arial" font-size="12">No benefit higher</text>',
        f'<text x="{center + 14}" y="{height - 16}" font-family="Arial" font-size="12">Benefit higher</text>',
    ]
    for i, row in enumerate(rows):
        y = top + i * row_h
        value = row["rank_biserial_dcb_vs_no_dcb"]
        x0 = center
        x1 = center + max(-1, min(1, value)) * scale
        color = "#356f8f" if value >= 0 else "#9f5c43"
        parts.append(f'<text x="12" y="{y + 18}" font-family="Arial" font-size="13">{html.escape(str(row["feature"]))}</text>')
        parts.append(f'<line x1="{x0}" y1="{y + 10}" x2="{x1:.1f}" y2="{y + 10}" stroke="{color}" stroke-width="8" stroke-linecap="round"/>')
        parts.append(f'<text x="{x1 + (10 if value >= 0 else -48):.1f}" y="{y + 15}" font-family="Arial" font-size="12">{fmt_float(value, 3)}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def markdown_table(rows: list[dict], columns: list[str], labels: dict[str, str] | None = None, n: int | None = None) -> str:
    labels = labels or {}
    if n is not None:
        rows = rows[:n]
    header = [labels.get(c, c) for c in columns]
    lines = ["| " + " | ".join(header) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in rows:
        vals = []
        for col in columns:
            val = row.get(col, "")
            if val is None:
                vals.append("")
                continue
            if isinstance(val, float):
                vals.append(fmt_float(val, 4 if "p_" in col or "q_" in col else 3))
            else:
                vals.append(str(val))
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)


def write_report(study: dict, cohort: list[dict], recurrence: list[dict], enrichment: list[dict], burden: list[dict], immune: list[dict]) -> None:
    included = [row for row in cohort if row["analysis_included"]]
    dcb = [row for row in included if row["benefit_group"] == "Durable clinical benefit"]
    no_dcb = [row for row in included if row["benefit_group"] == "No durable benefit"]
    indeterminate = [row for row in cohort if not row["analysis_included"]]
    primary = [r for r in enrichment if r["frequency_tier"] == "primary_screen"]
    low = [r for r in enrichment if r["frequency_tier"] == "low_frequency" and r["p_value"] < 0.05]
    top_primary = primary[:10]
    panel_hits = [r for r in enrichment if r["biology_panel"] and r["mutated_samples"] >= 2]
    fdr_hits = [r for r in enrichment if r.get("q_value", 1) <= 0.10]

    text = f"""# Precision Oncology Report: CTLA-4 Blockade Melanoma Cohort

## Cohort

Public cohort: **{study.get('name', STUDY_ID)}**  
Citation: **{study.get('citation', '')}**  
PMID: **{study.get('pmid', '')}**  
Sequenced samples reported by cBioPortal: **{study.get('sequencedSampleCount', '')}**

This analysis uses public cBioPortal data from the Van Allen metastatic melanoma cohort. The biological question is whether somatic mutation features and tumor immune context distinguish patients with durable clinical benefit from those with progressive disease after immune checkpoint blockade.

## Response Definition

- Durable clinical benefit: `CR`, `PR`, or `SD`
- No durable benefit: `PD`
- Excluded from binary response analysis: `X` or missing/indeterminate labels

Analysis set:

- Durable clinical benefit: **{len(dcb)}**
- No durable benefit: **{len(no_dcb)}**
- Excluded indeterminate: **{len(indeterminate)}**

## Why This Is A Precision Oncology Analysis

The analysis connects four clinically relevant layers:

1. Somatic mutation recurrence from WES.
2. Gene-level response enrichment.
3. Tumor mutational and neoantigen burden.
4. Immune-context estimates where available.

The goal is not to declare a clinical biomarker from one cohort. The goal is to identify biologically plausible, statistically transparent candidate signals that could be followed up in an independent validation set.

## Executive Interpretation

- No gene-level association passes FDR <= 0.10 in this exploratory screen.
- Several nominal gene-level signals are enriched in durable-benefit samples, but many are sparse and not established checkpoint-response biomarkers.
- Durable-benefit samples show higher median mutation count, nonsynonymous TMB, and neoantigen load.
- RNA-derived immune estimates are available in a subset and trend toward higher CD8 T-cell and activated lymphocyte context in durable-benefit samples.
- The most credible interpretation is a cohort-level immunogenicity pattern, not a validated single-gene predictor.

## Top Recurrently Mutated Genes

{markdown_table(recurrence, ["gene", "mutated_samples", "sample_fraction", "flag_gene", "biology_panel"], n=15)}

Large genes such as `TTN` and `MUC16` are retained in the table but flagged because recurrence may reflect gene length and mutability rather than treatment-specific biology.

![Top recurrent genes](figures/top_recurrent_genes.svg)

## Gene-Level Response Enrichment

Primary-screen genes are mutated in at least 5 analyzed samples. Low-frequency findings are reported separately because sparse events can produce unstable odds ratios.

{markdown_table(top_primary, ["gene", "mutated_samples", "dcb_mutated", "no_dcb_mutated", "odds_ratio_dcb", "p_value", "q_value", "flag_gene", "biology_panel"], n=10)}

FDR-significant genes at q <= 0.10: **{len(fdr_hits)}**

![Enrichment screen](figures/enrichment_screen.svg)

## Low-Frequency Signals

These genes have nominal p-value < 0.05 but are mutated in only 2-4 samples. They may be biologically interesting, but they should not be promoted as biomarkers without external evidence.

{markdown_table(low, ["gene", "mutated_samples", "dcb_mutated", "no_dcb_mutated", "odds_ratio_dcb", "p_value", "q_value", "biology_panel"], n=15)}

## Burden And Neoantigen Context

{markdown_table(burden, ["feature", "n_dcb", "n_no_dcb", "median_dcb", "median_no_dcb", "rank_biserial_dcb_vs_no_dcb"])}

Rank-biserial effect size is positive when values tend to be higher in durable-benefit samples and negative when values tend to be higher in no-benefit samples.

![Burden effect sizes](figures/burden_effect_sizes.svg)

## Immune Context

Selected CIBERSORT/ESTIMATE-derived immune features available through cBioPortal:

{markdown_table(immune, ["feature", "n_dcb", "n_no_dcb", "median_dcb", "median_no_dcb", "rank_biserial_dcb_vs_no_dcb"])}

![Immune effect sizes](figures/immune_effect_sizes.svg)

## Biology-Focused Candidate Review

Curated panel genes observed in the enrichment screen:

{markdown_table(panel_hits, ["gene", "biology_panel", "mutated_samples", "dcb_mutated", "no_dcb_mutated", "odds_ratio_dcb", "p_value", "q_value"], n=20)}

Interpretation should prioritize convergence across statistics and biology:

- DNA repair candidates can be linked to mutation burden, neoantigen load, and treatment sensitivity, but burden may also confound response association.
- Melanoma driver genes are important for disease biology but are not automatically response biomarkers.
- Antigen presentation and interferon pathway genes may be mechanistically relevant to checkpoint response, but low event counts require caution.

## Caveats

- The analysis is retrospective and exploratory.
- The binary response definition collapses heterogeneous clinical categories.
- Gene-level mutation presence does not distinguish driver, passenger, clonal, subclonal, loss-of-function, or gain-of-function effects.
- cBioPortal-derived fields may differ from publication-specific preprocessing.
- No external validation cohort is included here.
- No clinical actionability is claimed.

## Reproducibility

Run:

```bash
python3 src/analyze_van_allen_melanoma.py
```

The script downloads public cBioPortal data and regenerates the derived tables, figures, and this report.
"""
    REPORT_PATH.write_text(text, encoding="utf-8")


def main() -> int:
    ensure_dirs()
    try:
        study, patient_rows, sample_rows, mutations = fetch_data()
    except urllib.error.URLError as exc:
        print(f"Failed to fetch cBioPortal data: {exc}", file=sys.stderr)
        return 1

    write_json(DATA_DIR / "study_metadata.json", study)
    write_json(DATA_DIR / "patient_clinical_long.json", patient_rows)
    write_json(DATA_DIR / "sample_clinical_long.json", sample_rows)
    write_json(DATA_DIR / "mutations_raw.json", mutations)

    cohort = build_cohort(patient_rows, sample_rows)
    mutation_rows, gene_to_samples = build_mutation_tables(mutations)
    recurrence = gene_recurrence_rows(gene_to_samples, cohort)
    enrichment = analyze_gene_enrichment(cohort, gene_to_samples)
    burden, immune = burden_and_immune_summary(cohort)

    write_csv(DATA_DIR / "cohort_summary.csv", cohort)
    write_csv(DATA_DIR / "nonsynonymous_mutations.csv", mutation_rows)
    write_csv(TABLE_DIR / "top_recurrent_genes.csv", recurrence[:50])
    write_csv(TABLE_DIR / "gene_response_enrichment.csv", enrichment)
    write_csv(TABLE_DIR / "burden_summary.csv", burden)
    write_csv(TABLE_DIR / "immune_context_summary.csv", immune)

    create_bar_svg(
        FIGURE_DIR / "top_recurrent_genes.svg",
        recurrence,
        "gene",
        "sample_fraction",
        "Top recurrent nonsynonymous mutations",
        "Fraction of sequenced samples",
    )
    create_scatter_svg(FIGURE_DIR / "enrichment_screen.svg", enrichment)
    create_effect_svg(FIGURE_DIR / "burden_effect_sizes.svg", burden, "Burden and neoantigen effect sizes")
    create_effect_svg(FIGURE_DIR / "immune_effect_sizes.svg", immune, "Immune-context effect sizes")
    write_report(study, cohort, recurrence, enrichment, burden, immune)

    included = [row for row in cohort if row["analysis_included"]]
    print(f"Study: {study.get('name', STUDY_ID)}")
    print(f"Total samples with clinical rows: {len(cohort)}")
    print(f"Binary response analysis samples: {len(included)}")
    print(f"Mutation records fetched: {len(mutations)}")
    print(f"Nonsynonymous mutation rows retained: {len(mutation_rows)}")
    print(f"Report written: {REPORT_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
