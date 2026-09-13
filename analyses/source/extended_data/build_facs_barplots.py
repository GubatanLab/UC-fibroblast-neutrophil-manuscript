from pathlib import Path
import csv
import hashlib
import json
import math
import shutil
import sys

ROOT = Path(__file__).resolve().parents[2]
PKG = ROOT / "tmp/canonical_nature_20260903/packages"
sys.path.insert(0, str(PKG))

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from PIL import Image
import pymupdf as fitz


MM = 1 / 25.4
ACTIVE = ROOT / "Extended Data"
OUT = ROOT / "output/Nature_Extended_Data_FACS_Redraw_2026-09-03"
FIG_OUT = OUT / "extended_data"
SOURCE_OUT = OUT / "source_data"
REPORT_OUT = OUT / "reports"
ARCHIVE = ROOT / "output/Archive_before_FACS_barplot_redraw_2026-09-03/Extended Data"
QA = ROOT / "tmp/facs_barplot_redraw_20260903/qa"

for directory in [FIG_OUT, SOURCE_OUT, REPORT_OUT, ARCHIVE, QA]:
    directory.mkdir(parents=True, exist_ok=True)

# Manuscript-wide visual language.
INK = "#2F3337"
MUTED = "#6B7280"
GRID = "#E5E7EB"
CONTROL = "#3B75B9"
INFLAMED = "#D44B50"
DSS = "#E76F00"
BLOCKADE = "#1FA987"
DUAL = "#6F4C9B"
NEUTROPHIL = "#454B50"

plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 6.5,
        "axes.titlesize": 7.0,
        "axes.labelsize": 6.3,
        "xtick.labelsize": 5.7,
        "ytick.labelsize": 5.7,
        "text.color": INK,
        "axes.labelcolor": INK,
        "axes.edgecolor": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.65,
        "lines.linewidth": 0.8,
        "xtick.major.width": 0.55,
        "ytick.major.width": 0.55,
        "xtick.major.size": 2.4,
        "ytick.major.size": 2.4,
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
    }
)


def values_from_pixels(ys, zero_px, px_per_unit):
    return [round((zero_px - y) / px_per_unit, 4) for y in ys]


def sem(values):
    values = np.asarray(values, dtype=float)
    return float(np.std(values, ddof=1) / math.sqrt(len(values))) if len(values) > 1 else 0.0


def tint(hex_color, fraction=0.90):
    rgb = np.array(matplotlib.colors.to_rgb(hex_color))
    return tuple(rgb * (1 - fraction) + np.ones(3) * fraction)


def add_panel_header(ax, letter, title, y=1.13):
    ax.text(-0.12, y, letter, transform=ax.transAxes, ha="left", va="bottom", fontsize=8.0, fontweight="bold")
    ax.text(-0.03, y, title, transform=ax.transAxes, ha="left", va="bottom", fontsize=6.8, fontweight="bold")


def add_bracket(ax, left, right, y, label, tick, lw=0.75):
    ax.plot([left, left, right, right], [y - tick, y, y, y - tick], color=INK, lw=lw, clip_on=False)
    ax.text((left + right) / 2, y + tick * 0.22, label, ha="center", va="bottom", fontsize=6.4, fontweight="bold")


def draw_barplot(ax, groups, colors, ylim, yticks, ylabel, brackets, zero_line=False):
    x = np.arange(len(groups), dtype=float)
    for index, ((group_name, values), color) in enumerate(zip(groups.items(), colors)):
        mean = float(np.mean(values))
        error = sem(values)
        ax.bar(
            index,
            mean,
            width=0.58,
            facecolor=tint(color, 0.93),
            edgecolor=color,
            linewidth=1.0,
            zorder=1,
        )
        ax.errorbar(
            index,
            mean,
            yerr=error,
            fmt="none",
            ecolor=color,
            elinewidth=0.8,
            capsize=2.6,
            capthick=0.8,
            zorder=3,
        )
        offsets = np.linspace(-0.19, 0.19, len(values)) if len(values) > 1 else np.array([0.0])
        ax.scatter(
            np.full(len(values), index) + offsets,
            values,
            s=13,
            color=color,
            edgecolor="white",
            linewidth=0.25,
            zorder=4,
        )
    ax.set_xlim(-0.52, len(groups) - 0.48)
    ax.set_ylim(*ylim)
    ax.set_yticks(yticks)
    ax.set_ylabel(ylabel, labelpad=2.2)
    ax.set_xticks(x)
    ax.set_xticklabels([])
    ax.tick_params(axis="x", length=0)
    ax.tick_params(axis="y", pad=1.5)
    ax.spines["bottom"].set_position(("data", 0)) if zero_line else None
    ax.axhline(0, color=INK, lw=0.65, zorder=0) if zero_line else None
    for left, right, y, label, tick in brackets:
        add_bracket(ax, left, right, y, label, tick)


def save_figure(fig, number):
    stem = FIG_OUT / f"Extended_Data_Figure_{number:02d}"
    metadata = {
        "Title": f"Extended Data Figure {number}",
        "Author": "UC fibroblast-neutrophil study",
        "Subject": "Nature-style native redraw of source-displayed FACS/NET barplots; no new inferential analysis",
    }
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches=None, metadata=metadata)
    fig.savefig(stem.with_suffix(".png"), dpi=450, bbox_inches=None)
    plt.close(fig)
    with Image.open(stem.with_suffix(".png")) as image:
        image.save(stem.with_suffix(".tiff"), compression="tiff_lzw", dpi=(450, 450))
    if stem.with_suffix(".tiff").stat().st_size > 10_000_000:
        with Image.open(stem.with_suffix(".png")) as image:
            target = (round(image.width * 300 / 450), round(image.height * 300 / 450))
            image.resize(target, Image.Resampling.LANCZOS).save(
                stem.with_suffix(".tiff"), compression="tiff_lzw", dpi=(300, 300)
            )
    doc = fitz.open(stem.with_suffix(".pdf"))
    assert len(doc) == 1
    doc[0].get_pixmap(matrix=fitz.Matrix(2.3, 2.3), alpha=False).save(QA / f"ED{number:02d}.png")
    return stem


# Display-coordinate transcription of the active source-only graphics. Values are
# recovered from the plotted dot centers and calibrated to the printed y-axis ticks.
# Coordinates obscured by coincident points or bar outlines were visually resolved
# at the source graphic's plotting precision; they are not raw FCS measurements.
ed9_pixels = {
    "a": {
        "zero": 714.0,
        "ppu": 12.25,
        "groups": {
            "Neutrophils only": [578.0, 580.0, 591.5, 602.0, 603.0, 604.0],
            "+ control fibroblasts": [581.3, 595.1, 600.9, 614.4, 623.0, 647.3],
            "+ UC fibroblasts": [298.8, 399.7, 408.5, 446.0, 473.0, 473.8],
        },
    },
    "b": {
        "zero": 706.0,
        "ppu": 18.9,
        "groups": {
            "Neutrophils only": [530.0, 611.5, 626.6, 641.7, 649.3, 668.2],
            "+ control fibroblasts": [558.6, 583.2, 639.9, 666.3, 683.3],
            "+ UC fibroblasts": [374.2, 419.5, 421.5, 424.3, 476.8, 489.5],
        },
    },
    "c": {
        "zero": 700.5,
        "ppu": 14.9,
        "groups": {
            "Neutrophils only": [605.7, 606.8, 627.0, 685.3, 685.9, 687.5],
            "+ control fibroblasts": [581.2, 597.0, 604.6, 623.5, 634.2, 667.8],
            "+ UC fibroblasts": [372.9, 432.0, 432.4, 462.2, 490.5, 495.5],
        },
    },
    "d": {
        "zero": 714.0,
        "ppu": 17.6,
        "groups": {
            "Neutrophils only": [568.2, 619.5, 632.7, 664.5, 671.4, 692.8],
            "+ control fibroblasts": [496.1, 503.3, 510.1, 554.0, 601.7, 618.1],
            "+ UC fibroblasts": [379.0, 389.3, 397.7, 496.1, 498.5, 504.0],
        },
    },
    "e": {
        "zero": 721.5,
        "ppu": 11.65,
        "groups": {
            "Neutrophils only": [379.8, 471.6, 483.0, 489.0, 524.5, 546.2],
            "+ control fibroblasts": [451.5, 473.5, 494.7, 547.0, 573.0, 574.0],
            "+ UC fibroblasts": [567.9, 576.7, 610.7, 611.0, 634.0, 675.8],
        },
    },
}

ed10_pixels = {
    "a": {
        "zero": 783.5,
        "ppu": 225.0,
        "groups": {
            "Control": [253.4, 316.7, 343.8, 522.0, 585.2, 628.5],
            "DSS + PBS": [325.8, 452.1, 472.4, 483.7, 551.4, 590.0, 652.9, 687.0, 725.2, 736.5],
            "DSS + α5β1 inhibitor": [650.7, 675.9, 720.4, 721.0, 725.0],
            "DSS + α5β1/NAMPT inhibitors": [709.5, 730.2, 730.0, 783.5, 783.5],
        },
    },
    "b": {
        "zero": 774.5,
        "ppu": 73.75,
        "groups": {
            "Control": [343.5, 362.5, 380.9, 497.1, 500.1, 566.4],
            "DSS + PBS": [362.5, 386.0, 401.1, 466.9, 490.9, 495.0, 516.4, 538.4, 570.2, 576.9, 602.5, 603.1, 698.6],
            "DSS + α5β1 inhibitor": [590.1, 605.2, 659.0, 669.0, 717.7, 733.5],
            "DSS + α5β1/NAMPT inhibitors": [623.0, 649.8, 706.5, 726.2, 739.5, 758.0],
        },
    },
    "c": {
        "zero": 2000.5,
        "ppu": 33.7,
        "groups": {
            "Control": [1481.8, 1589.7, 1748.0, 1776.0, 1901.0, 1916.0],
            "DSS + PBS": [1502.3, 1609.3, 1626.5, 1636.3, 1639.9, 1701.0, 1704.0, 1713.0, 1744.2, 1780.5, 1804.5, 1957.4],
            "DSS + α5β1 inhibitor": [1696.0, 1725.8, 1840.0, 1906.0, 1923.6, 1964.0],
            "DSS + α5β1/NAMPT inhibitors": [1708.4, 1747.9, 1784.7, 1809.5, 1949.7, 1954.0],
        },
    },
    "d": {
        "zero": 2003.0,
        "ppu": 45.3,
        "groups": {
            "Control": [1520.7, 1655.0, 1670.1, 1684.4, 1792.5, 1888.9],
            "DSS + PBS": [1626.4, 1650.7, 1693.3, 1717.0, 1720.6, 1749.2, 1813.3, 1847.3, 1854.7, 1945.5],
            "DSS + α5β1 inhibitor": [1622.0, 1924.5, 1942.0, 1949.0, 1979.2, 1998.7],
            "DSS + α5β1/NAMPT inhibitors": [1832.5, 1866.2, 1888.0, 1890.0, 1911.0, 1930.0],
        },
    },
    "e": {
        "zero": 2003.5,
        "ppu": 6.85,
        "groups": {
            "Control": [1603.9, 1636.0, 1660.1, 1702.2, 1723.5, 1771.6],
            "DSS + PBS": [1545.3, 1607.9, 1632.0, 1639.0, 1641.0, 1642.0, 1654.0, 1668.0, 1680.0, 1680.0],
            "DSS + α5β1 inhibitor": [1484.7, 1654.6, 1693.0, 1720.0, 1754.0, 1763.0],
            "DSS + α5β1/NAMPT inhibitors": [1615.0, 1615.0, 1624.0, 1679.0, 1750.5, 1757.5],
        },
    },
    "f": {
        "zero": 2559.5,
        "ppu": 0.034,
        "groups": {
            "Control": [2392.5, 2416.0, 2419.7, 2436.0, 2512.1, 2535.4],
            "DSS + PBS": [2299.6, 2386.4, 2399.0, 2403.0, 2407.5, 2416.0, 2426.0, 2426.0, 2433.0, 2451.7, 2471.0, 2473.0],
            "DSS + α5β1 inhibitor": [2477.4, 2509.2, 2531.5, 2541.0, 2565.7, 2576.2],
            "DSS + α5β1/NAMPT inhibitors": [2514.7, 2520.9, 2521.3, 2543.0, 2551.0, 2554.0],
        },
    },
    "g": {
        "zero": 2565.5,
        "ppu": 0.05575,
        "groups": {
            "Control": [2362.2, 2375.7, 2428.0, 2436.0, 2437.5, 2532.4],
            "DSS + PBS": [2271.6, 2279.9, 2364.5, 2367.6, 2372.6, 2382.5, 2398.0, 2431.0, 2433.0, 2441.0, 2442.0, 2450.0],
            "DSS + α5β1 inhibitor": [2506.4, 2542.7, 2543.0, 2567.0, 2633.2, 2655.8],
            "DSS + α5β1/NAMPT inhibitors": [2395.3, 2539.3, 2586.6, 2608.7, 2647.3],
        },
    },
    "h": {
        "zero": 2603.5,
        "ppu": 0.0073,
        "groups": {
            "Control": [2591.0, 2592.0, 2593.0, 2597.0, 2598.0, 2601.0],
            "DSS + PBS": [2313.9, 2358.3, 2393.3, 2398.0, 2403.0, 2563.1, 2582.0, 2585.0, 2590.0, 2591.0, 2591.4, 2592.0],
            "DSS + α5β1 inhibitor": [2395.7, 2443.7, 2450.5, 2539.8, 2563.8, 2587.7],
            "DSS + α5β1/NAMPT inhibitors": [2391.7, 2393.4, 2401.8, 2423.6, 2546.4, 2606.2],
        },
    },
}


def convert_panels(pixel_panels):
    result = {}
    for panel, specification in pixel_panels.items():
        result[panel] = {
            group: values_from_pixels(ys, specification["zero"], specification["ppu"])
            for group, ys in specification["groups"].items()
        }
    return result


ed9 = convert_panels(ed9_pixels)
ed10 = convert_panels(ed10_pixels)


def build_ed9():
    fig = plt.figure(figsize=(180 * MM, 75 * MM))
    gs = fig.add_gridspec(1, 5, left=0.055, right=0.992, bottom=0.185, top=0.765, wspace=0.48)
    axes = [fig.add_subplot(gs[0, i]) for i in range(5)]
    colors = [NEUTROPHIL, CONTROL, INFLAMED]
    titles = ["NETosis assay", "PADI4+ neutrophils", "OSM+ neutrophils", "MX1+ neutrophils", "CXCR4+ neutrophils"]
    settings = [
        ((0, 44), [0, 10, 20, 30, 40], "NET-derived elastase\n(mU ml$^{-1}$ per $10^6$ cells)", [(0, 1, 37.4, "ns", 1.0), (1, 2, 37.4, "****", 1.0)]),
        ((0, 25), [0, 5, 10, 15, 20], "Positive neutrophils (%)", [(0, 1, 20.5, "ns", 0.6), (1, 2, 20.5, "**", 0.6), (0, 2, 24.0, "**", 0.6)]),
        ((0, 30), [0, 5, 10, 15, 20, 25], "Positive neutrophils (%)", [(0, 1, 25.4, "ns", 0.7), (1, 2, 25.4, "****", 0.7), (0, 2, 28.7, "****", 0.7)]),
        ((0, 29), [0, 5, 10, 15, 20, 25], "Positive neutrophils (%)", [(0, 1, 22.3, "ns", 0.7), (1, 2, 22.3, "ns", 0.7), (0, 2, 27.5, "***", 0.7)]),
        ((0, 45), [0, 10, 20, 30, 40], "Positive neutrophils (%)", [(0, 1, 34.2, "ns", 1.0), (1, 2, 34.2, "*", 1.0), (0, 2, 43.0, "**", 1.0)]),
    ]
    for index, (ax, key, title, setting) in enumerate(zip(axes, "abcde", titles, settings)):
        ylim, yticks, ylabel, brackets = setting
        draw_barplot(ax, ed9[key], colors, ylim, yticks, ylabel, brackets)
        add_panel_header(ax, key, title)
    handles = [
        Line2D([0], [0], marker="o", color=c, markerfacecolor=c, markersize=4.2, lw=1.1, label=l)
        for c, l in zip(colors, ["Neutrophils only", "+ control fibroblasts", "+ UC fibroblasts"])
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.525, 0.955), ncol=3, frameon=False, handletextpad=0.5, columnspacing=1.8, fontsize=6.2)
    fig.text(0.995, 0.018, "Bars show mean ± s.e.m.; significance labels copied from the source display.", ha="right", va="bottom", fontsize=5.2, color=MUTED)
    return save_figure(fig, 9)


def build_ed10():
    fig = plt.figure(figsize=(180 * MM, 170 * MM))
    outer = fig.add_gridspec(3, 1, left=0.06, right=0.992, bottom=0.07, top=0.905, height_ratios=[1.05, 1, 1], hspace=0.55)
    top = outer[0].subgridspec(1, 2, wspace=0.32)
    middle = outer[1].subgridspec(1, 3, wspace=0.42)
    bottom = outer[2].subgridspec(1, 3, wspace=0.42)
    axes = [fig.add_subplot(top[0, i]) for i in range(2)]
    axes += [fig.add_subplot(middle[0, i]) for i in range(3)]
    axes += [fig.add_subplot(bottom[0, i]) for i in range(3)]
    colors = [CONTROL, DSS, BLOCKADE, DUAL]
    titles = [
        "α5β1+ fibroblasts",
        "α5β1+FAP+ fibroblasts",
        "OSM+ neutrophils",
        "PADI4+ neutrophils",
        "CXCR4+ neutrophils",
        "MPO in OSM+ cells",
        "MPO in PADI4+ cells",
        "MPO in CXCR4+ cells",
    ]
    settings = [
        ((0, 2.75), [0, 0.5, 1.0, 1.5, 2.0, 2.5], "PDPN+ α5β1+\nfibroblasts (%)", [(1, 2, 2.55, "*", 0.07)], False),
        ((0, 8.4), [0, 2, 4, 6, 8], "PDPN+ α5β1+FAP+\nfibroblasts (%)", [(1, 2, 6.65, "**", 0.18)], False),
        ((0, 20.5), [0, 5, 10, 15, 20], "Positive neutrophils (%)", [(1, 2, 17.8, "*", 0.45)], False),
        ((0, 15.5), [0, 5, 10, 15], "Positive neutrophils (%)", [(1, 2, 12.6, "*", 0.35)], False),
        ((0, 100), [0, 20, 40, 60, 80, 100], "Positive neutrophils (%)", [(1, 2, 88, "ns", 2.2)], False),
        ((-2000, 11000), [-2000, 0, 2000, 4000, 6000, 8000, 10000], "MPO MFI", [(1, 2, 9300, "***", 320)], True),
        ((-2000, 6500), [-2000, 0, 2000, 4000, 6000], "MPO MFI", [(1, 2, 5600, "***", 210)], True),
        ((-10000, 52000), [-10000, 0, 10000, 20000, 30000, 40000, 50000], "MPO MFI", [(1, 2, 45500, "ns", 1500)], True),
    ]
    for ax, key, title, setting in zip(axes, "abcdefgh", titles, settings):
        ylim, yticks, ylabel, brackets, zero_line = setting
        draw_barplot(ax, ed10[key], colors, ylim, yticks, ylabel, brackets, zero_line=zero_line)
        add_panel_header(ax, key, title, y=1.10)
        if key in "fgh":
            ax.ticklabel_format(style="plain", axis="y")
    handles = [
        Line2D([0], [0], marker="o", color=c, markerfacecolor=c, markersize=4.1, lw=1.1, label=l)
        for c, l in zip(colors, ["Control", "DSS + PBS", "DSS + α5β1 inhibitor", "DSS + α5β1/NAMPT inhibitors"])
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.525, 0.976), ncol=4, frameon=False, handletextpad=0.45, columnspacing=1.3, fontsize=5.85)
    fig.text(0.992, 0.012, "Bars show mean ± s.e.m.; significance labels copied from the source display.", ha="right", va="bottom", fontsize=5.2, color=MUTED)
    return save_figure(fig, 10)


stem9 = build_ed9()
stem10 = build_ed10()

# Machine-readable display transcription. These are plot-reconstruction values,
# explicitly not claimed as raw-event or instrument-export source data.
rows = []
panel_endpoints = {
    9: {"a": "NET-derived elastase", "b": "PADI4+ neutrophils", "c": "OSM+ neutrophils", "d": "MX1+ neutrophils", "e": "CXCR4+ neutrophils"},
    10: {"a": "α5β1+ fibroblasts", "b": "α5β1+FAP+ fibroblasts", "c": "OSM+ neutrophils", "d": "PADI4+ neutrophils", "e": "CXCR4+ neutrophils", "f": "MPO MFI in OSM+ cells", "g": "MPO MFI in PADI4+ cells", "h": "MPO MFI in CXCR4+ cells"},
}
for figure, pixel_data, value_data in [(9, ed9_pixels, ed9), (10, ed10_pixels, ed10)]:
    for panel, groups in value_data.items():
        for group, values in groups.items():
            ys = pixel_data[panel]["groups"][group]
            for replicate, (value, source_y) in enumerate(zip(values, ys), 1):
                rows.append(
                    {
                        "figure": figure,
                        "panel": panel,
                        "endpoint": panel_endpoints[figure][panel],
                        "group": group,
                        "display_point": replicate,
                        "display_transcribed_value": value,
                        "source_y_pixel": source_y,
                        "source_zero_y_pixel": pixel_data[panel]["zero"],
                        "source_pixels_per_unit": pixel_data[panel]["ppu"],
                        "provenance": "digitized_from_canonical_source_display; not raw FCS data",
                    }
                )
with (SOURCE_OUT / "Extended_Data_Figures_09_10_display_transcribed_values.csv").open("w", newline="", encoding="utf-8-sig") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)

# Updated legends accurately describe the native redraw and its limitations.
legend9 = """Extended Data Fig. 9 | NET-associated release and neutrophil protein-marker assays.\n\n(a) NET-derived elastase release after neutrophils alone, control-fibroblast coculture or UC-fibroblast coculture. (b–e) Percentages of PADI4-, OSM-, MX1- and CXCR4-positive neutrophils, respectively, across the same three conditions. Individual observations were transcribed from the canonical source display and replotted as native vectors; bars show the reconstructed mean ± s.e.m. Existing brackets and significance labels were copied verbatim, without recalculating P values. Single-marker gates and NET-associated release are complementary to multigene RNA states, not interchangeable measurements. No FACS reanalysis or new hypothesis testing was performed. [AUTHOR TO CONFIRM against primary assay records before submission: exact n, summary/error-bar definition, statistical test, correction family, exact adjusted P values, significance-key thresholds and parent gating denominator.]\n"""
legend10 = """Extended Data Fig. 10 | In-vivo stromal, neutrophil and MPO protein readouts.\n\n(a,b) Detectable α5β1-positive fractions among PDPN-positive and PDPN-positive/FAP-positive fibroblasts, respectively, in control, DSS+PBS, DSS+α5β1-inhibitor and DSS+combined α5β1/NAMPT-inhibitor groups. Reduced detection is compatible with target engagement but may reflect receptor occupancy or epitope interference. (c–e) OSM-, PADI4- and CXCR4-positive neutrophil frequencies. (f–h) MPO mean fluorescence intensity within the corresponding subsets. Individual observations were transcribed from the canonical source display and replotted as native vectors; bars show the reconstructed mean ± s.e.m. Existing brackets and significance labels were copied verbatim, without recalculating P values. CXCR4 frequency and MPO comparisons remain nonsignificant in the displayed source comparisons. Points represent biological mice. No FACS reanalysis, new hypothesis testing, or additive/synergistic inference was performed. [AUTHOR TO CONFIRM against primary assay records before submission: exact n, summary/error-bar definition, gating hierarchy, MFI transformation/background subtraction, statistical test, correction family, exact adjusted P values and significance-key thresholds.]\n"""
(FIG_OUT / "Extended_Data_Figure_09_legend.md").write_text(legend9, encoding="utf-8")
(FIG_OUT / "Extended_Data_Figure_10_legend.md").write_text(legend10, encoding="utf-8")

# Package the unchanged ED1–ED8 pages alongside the revised ED9–ED10 pages.
for number in range(1, 9):
    for suffix in [".pdf", ".png", ".tiff"]:
        shutil.copy2(ACTIVE / f"Extended_Data_Figure_{number:02d}{suffix}", FIG_OUT / f"Extended_Data_Figure_{number:02d}{suffix}")
    legend = ACTIVE / f"Extended_Data_Figure_{number:02d}_legend.md"
    if legend.exists():
        shutil.copy2(legend, FIG_OUT / legend.name)

combined = fitz.open()
for number in range(1, 11):
    source = FIG_OUT / f"Extended_Data_Figure_{number:02d}.pdf"
    with fitz.open(source) as document:
        combined.insert_pdf(document)
combined.set_metadata(
    {
        "title": "Extended Data Figures 1–10",
        "author": "UC fibroblast-neutrophil study",
        "subject": "Nature-style Extended Data set; Figures 9–10 are native display-transcribed redraws",
    }
)
combined_path = OUT / "Extended_Data_Figures_1_to_10.pdf"
combined.save(combined_path, garbage=4, deflate=True)

legends = []
for number in range(1, 11):
    legends.append((FIG_OUT / f"Extended_Data_Figure_{number:02d}_legend.md").read_text(encoding="utf-8").strip())
all_legends = "\n\n".join(legends) + "\n"
(OUT / "Extended_Data_Figure_Legends_1_to_10.md").write_text(all_legends, encoding="utf-8")
(FIG_OUT / "Extended_Data_Figure_Legends_1_to_10.md").write_text(all_legends, encoding="utf-8")

provenance = """# FACS barplot redraw provenance\n\nExtended Data Figures 9 and 10 were rebuilt as native vector charts to match the manuscript's Nature-style figure system. The canonical source graphics were raster-only. The available mouse Prism files describe a separate FAP/GCV experiment and were deliberately not substituted for the displayed α5β1-blockade groups.\n\nIndividual plotted observations were transcribed from the canonical source displays using their printed axis calibrations. Coincident points or points touching bar outlines were resolved at the plotting precision of the source graphic. Bars and error bars were reconstructed descriptively as mean ± s.e.m., which visually matches the source convention but remains subject to author confirmation from the primary assay records. Existing significance brackets and labels were copied verbatim. No raw FCS files were regated; no group labels were reassigned; no P values, multiple-testing adjustments, or interaction tests were calculated.\n\nThe CSV in `source_data` is a display-transcription table for figure reconstruction, not a replacement for primary instrument exports or raw-event source data. Before submission, the corresponding author should reconcile it against the original assay table and complete the bracketed reporting details in the legends.\n"""
(REPORT_OUT / "FACS_REDRAW_PROVENANCE.md").write_text(provenance, encoding="utf-8")

# Recoverable promotion into the active Extended Data set.
promotion_targets = [
    ACTIVE / "Extended_Data_Figure_09.pdf",
    ACTIVE / "Extended_Data_Figure_09.png",
    ACTIVE / "Extended_Data_Figure_09.tiff",
    ACTIVE / "Extended_Data_Figure_09_legend.md",
    ACTIVE / "Extended_Data_Figure_10.pdf",
    ACTIVE / "Extended_Data_Figure_10.png",
    ACTIVE / "Extended_Data_Figure_10.tiff",
    ACTIVE / "Extended_Data_Figure_10_legend.md",
    ACTIVE / "Extended_Data_Figure_Legends_1_to_10.md",
    ACTIVE / "Extended_Data_Figures_1_to_10.pdf",
]
for target in promotion_targets:
    if target.exists() and not (ARCHIVE / target.name).exists():
        shutil.copy2(target, ARCHIVE / target.name)

for number in [9, 10]:
    for suffix in [".pdf", ".png", ".tiff"]:
        shutil.copy2(FIG_OUT / f"Extended_Data_Figure_{number:02d}{suffix}", ACTIVE / f"Extended_Data_Figure_{number:02d}{suffix}")
    shutil.copy2(FIG_OUT / f"Extended_Data_Figure_{number:02d}_legend.md", ACTIVE / f"Extended_Data_Figure_{number:02d}_legend.md")
shutil.copy2(OUT / "Extended_Data_Figure_Legends_1_to_10.md", ACTIVE / "Extended_Data_Figure_Legends_1_to_10.md")
shutil.copy2(combined_path, ACTIVE / "Extended_Data_Figures_1_to_10.pdf")

inventory = {}
for number in [9, 10]:
    pdf = FIG_OUT / f"Extended_Data_Figure_{number:02d}.pdf"
    png = FIG_OUT / f"Extended_Data_Figure_{number:02d}.png"
    tiff = FIG_OUT / f"Extended_Data_Figure_{number:02d}.tiff"
    with fitz.open(pdf) as document:
        page = document[0]
        inventory[str(number)] = {
            "pages": len(document),
            "width_mm": round(page.rect.width * 25.4 / 72, 2),
            "height_mm": round(page.rect.height * 25.4 / 72, 2),
            "pdf_bytes": pdf.stat().st_size,
            "png_bytes": png.stat().st_size,
            "tiff_bytes": tiff.stat().st_size,
            "pdf_sha256": hashlib.sha256(pdf.read_bytes()).hexdigest(),
        }
with fitz.open(combined_path) as document:
    inventory["combined"] = {
        "pages": len(document),
        "bytes": combined_path.stat().st_size,
        "sha256": hashlib.sha256(combined_path.read_bytes()).hexdigest(),
    }
(REPORT_OUT / "build_manifest.json").write_text(
    json.dumps(
        {
            "canonical_promotion": True,
            "figures_redrawn": [9, 10],
            "source_mode": "display transcription from canonical raster source",
            "raw_fcs_reanalysis": False,
            "new_inferential_statistics": False,
            "archived_previous_active_files": str(ARCHIVE),
            "inventory": inventory,
        },
        indent=2,
    ),
    encoding="utf-8",
)

assert inventory["9"]["pages"] == 1 and inventory["10"]["pages"] == 1
assert inventory["combined"]["pages"] == 10
assert inventory["9"]["tiff_bytes"] <= 10_000_000 and inventory["10"]["tiff_bytes"] <= 10_000_000
print(json.dumps({"output": str(OUT), "active": str(ACTIVE), "inventory": inventory}, indent=2))
