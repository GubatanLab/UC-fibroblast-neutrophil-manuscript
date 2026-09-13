import concurrent.futures
import json
import pathlib
import time
import urllib.parse
import urllib.request


ROOT = pathlib.Path(r"input_data/mouse\Figure 6 Spatial External Validation\data\SCP3818")
ROOT.mkdir(parents=True, exist_ok=True)

# Locked before analysis: fibroblast identity/adhesion/recruitment and neutrophil/feedback modules.
GENES = [
    "FAP", "PDPN", "THY1", "ITGA5", "ITGB1", "FN1",
    "CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL6", "CXCL8", "CSF3", "IL6", "ICAM1", "OSMR",
    "S100A8", "S100A9", "CSF3R", "CXCR1", "CXCR2", "FCGR3B", "CEACAM8", "MPO",
    "OSM", "NAMPT", "IL1B", "TNF",
]
SECTIONS = ["UC1 inflamed", "UC1 less inflamed"]
BASE = "https://singlecell.broadinstitute.org/single_cell/api/v1/studies/SCP3818/expression/violin"


def fetch(section, gene):
    query = urllib.parse.urlencode({
        "cluster": section,
        "genes": gene,
        "annotation_name": "cell_types",
        "annotation_type": "group",
        "annotation_scope": "study",
    })
    req = urllib.request.Request(BASE + "?" + query, headers={"Accept": "application/json", "User-Agent": "Mozilla/5.0"})
    last = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=300) as response:
                payload = json.load(response)
            if gene not in payload.get("gene_names", []):
                return section, gene, "absent", None
            rows = []
            for cell_type, block in payload["values"].items():
                ys = block.get("y", [])
                cells = block.get("cells", [])
                rows.extend([str(c), cell_type, float(y)] for c, y in zip(cells, ys))
            return section, gene, "ok", rows
        except Exception as exc:
            last = repr(exc)
            time.sleep(2 ** attempt)
    return section, gene, "error", last


results = {}
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
    futures = [pool.submit(fetch, section, gene) for section in SECTIONS for gene in GENES]
    for future in concurrent.futures.as_completed(futures):
        section, gene, status, value = future.result()
        results[(section, gene)] = (status, value)
        print(section, gene, status, flush=True)

manifest = []
for section in SECTIONS:
    safe_section = section.replace(" ", "_")
    for gene in GENES:
        status, value = results[(section, gene)]
        manifest.append({"section": section, "gene": gene, "status": status, "detail": value if status == "error" else None})
        if status == "ok":
            out = ROOT / f"{safe_section}_{gene}.tsv"
            with out.open("w", encoding="utf-8", newline="") as handle:
                handle.write("cell_id\tcell_type\texpression\n")
                for cell_id, cell_type, expression in value:
                    handle.write(f"{cell_id}\t{cell_type}\t{expression}\n")

(ROOT / "expression_download_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
print(json.dumps(manifest, indent=2))
