import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


# This script generates a cohort-specific lollipop plot for variants in a
# selected target gene. Variants are positioned according to their genomic
# location and displayed according to the number of carriers in the three
# study cohorts.
#
# Gene-specific coordinates and exon/UTR positions can be adjusted in the
# configuration section below.


# ___ INPUT AND OUTPUT ___

EXCEL_FILE = "variant_data.xlsx"

# Columns containing carrier IDs for the three cohorts.
COHORT_COLUMNS = {
    "N0": "PATIENTS_HEALTHY",
    "N1": "PATIENTS_1xCANCER",
    "N2": "PATIENTS_2xCANCER",
}


# ___ GENE-SPECIFIC SETTINGS ___
# Adjust these parameters for the gene to be visualized.
# The configuration below corresponds to XRCC6 (GRCh38).

GENE = "XRCC6"

GENE_START = 41621295
GENE_END = 41664041

PLOT_MIN_X = 41616000
PLOT_MAX_X = 41670000

# 5' UTR regions
UTR_5_BLOCKS = [
    (41621295, 41621345),
    (41621990, 41622004),
]

# 3' UTR regions
UTR_3_BLOCKS = [
    (41663816, 41664041),
]

# Coding exon regions
CODING_EXONS = [
    (41622005, 41622086),
    (41628118, 41628230),
    (41636113, 41636251),
    (41636516, 41636770),
    (41637608, 41637791),
    (41646896, 41647082),
    (41650723, 41650891),
    (41653529, 41653690),
    (41656903, 41657032),
    (41658252, 41658352),
    (41661331, 41661444),
    (41663622, 41663815),
]

# Labels for the genomic regions shown on the left and right side.
# Adjust these if required for the orientation of the selected gene.
LEFT_FLANK_LABEL = "Upstream"
RIGHT_FLANK_LABEL = "Downstream"


# ___ DERIVED SETTINGS ___

SHEET_NAME = GENE

OUTPUT_PREFIX = GENE.lower()
VARIANT_KEY_FILE = f"{OUTPUT_PREFIX}_variant_key.xlsx"
PLOT_FILE = f"{OUTPUT_PREFIX}_lollipop_plot.png"


# ___ FONT SETTINGS ___
# Calibri is used if available on the system.

plt.rcParams["font.family"] = "Calibri"


# ___ LOAD AND PREPARE DATA ___

if not os.path.exists(EXCEL_FILE):
    raise FileNotFoundError(
        f"Excel file '{EXCEL_FILE}' was not found. "
        "Place the file in the same directory as the script "
        "or adjust EXCEL_FILE."
    )

df = pd.read_excel(
    EXCEL_FILE,
    sheet_name=SHEET_NAME
)

# Fill merged cells containing variant annotation information.
cols_to_fill = [
    "Location",
    "Consequence",
    "Existing_variation",
    "IMPACT",
    "VARIANT_CLASS",
    "gnomAD AF",
    "CLINVAR",
]

existing_cols_to_fill = [
    col
    for col in cols_to_fill
    if col in df.columns
]

df[existing_cols_to_fill] = df[existing_cols_to_fill].ffill()

# Extract GRCh38 genomic position from the Location column.
df["pos_38"] = df["Location"].apply(
    lambda x: int(str(x).split(":")[1].strip())
)


# ___ COUNT CARRIERS PER COHORT ___

def count_patients(cell_value):
    """Count comma-separated patient IDs in an Excel cell."""

    if pd.isna(cell_value):
        return 0

    value = str(cell_value).strip()

    if not value or value in {"-", "nan"}:
        return 0

    patients = [
        patient.replace("*", "").strip()
        for patient in value.split(",")
        if patient.strip()
    ]

    return len(patients)


required_patient_columns = list(COHORT_COLUMNS.values())

missing_columns = [
    col
    for col in required_patient_columns
    if col not in df.columns
]

if missing_columns:
    raise KeyError(
        "The following columns are missing from the Excel file: "
        + ", ".join(missing_columns)
    )

df["N0"] = df[COHORT_COLUMNS["N0"]].apply(count_patients)
df["N1"] = df[COHORT_COLUMNS["N1"]].apply(count_patients)
df["N2"] = df[COHORT_COLUMNS["N2"]].apply(count_patients)


# ___ GROUP VARIANTS ___
# Combine multiple rows corresponding to the same genomic position.

agg_rules = {
    "VARIANT_CLASS": "first",
    "Consequence": "first",
    "N0": "sum",
    "N1": "sum",
    "N2": "sum",
}

grouped = (
    df.groupby("pos_38")
    .agg(agg_rules)
    .reset_index()
)


# ___ PREPARE VARIANTS FOR PLOTTING ___

variants = []

for _, row in grouped.iterrows():

    pos = int(row["pos_38"])
    var_class = str(row["VARIANT_CLASS"]).strip()
    consequence = str(row["Consequence"]).strip().lower()

    class_lower = var_class.lower()

    if "sequence_alteration" in class_lower:
        class_clean = "seq. alter"
    elif "insertion" in class_lower:
        class_clean = "ins"
    elif "deletion" in class_lower:
        class_clean = "del"
    else:
        class_clean = var_class

    if "downstream" in consequence:
        consequence_key = "downstream"
    elif "3_prime" in consequence:
        consequence_key = "3_prime_UTR"
    elif "5_prime" in consequence:
        consequence_key = "5_prime_UTR"
    else:
        consequence_key = "intron"

    variants.append({
        "pos": pos,
        "class": class_clean,
        "consequence": consequence_key,
        "cohorts": {
            "N0": int(row["N0"]),
            "N1": int(row["N1"]),
            "N2": int(row["N2"]),
        },
    })

variants = sorted(
    variants,
    key=lambda variant: variant["pos"]
)


# ___ VARIANT SPACING ___
# Spread nearby variants horizontally while retaining their connection
# to the original genomic position.

n = len(variants)

if n:

    spread_left = PLOT_MIN_X + 2500
    spread_right = PLOT_MAX_X - 2500

    min_dist = min(
        4000,
        (spread_right - spread_left) / max(n - 1, 1)
    )

    x_spaced = np.array(
        [variant["pos"] for variant in variants],
        dtype=float
    )

    for _ in range(50):

        x_spaced[0] = max(
            x_spaced[0],
            spread_left
        )

        for i in range(1, n):
            x_spaced[i] = max(
                x_spaced[i],
                x_spaced[i - 1] + min_dist
            )

        x_spaced[-1] = min(
            x_spaced[-1],
            spread_right
        )

        for i in range(n - 2, -1, -1):
            x_spaced[i] = min(
                x_spaced[i],
                x_spaced[i + 1] - min_dist
            )

    for variant, x in zip(variants, x_spaced):
        variant["pos_spaced"] = float(x)


# ___ PLOT SETTINGS ___

cohort_colors = {
    "N0": "#81C784",
    "N1": "#FFB74D",
    "N2": "#E57373",
}

fig, ax = plt.subplots(
    figsize=(9.1, 5.7),
    dpi=300
)


# ___ GENE STRUCTURE ___

# Upstream and downstream regions
ax.plot(
    [PLOT_MIN_X, GENE_START],
    [0, 0],
    color="#BDBDBD",
    linestyle="--",
    linewidth=1.8,
    zorder=1
)

ax.plot(
    [GENE_END, PLOT_MAX_X],
    [0, 0],
    color="#BDBDBD",
    linestyle="--",
    linewidth=1.8,
    zorder=1
)

# Gene body
ax.plot(
    [GENE_START, GENE_END],
    [0, 0],
    color="#616161",
    linewidth=2.5,
    zorder=2
)

# UTR regions
for start, end in UTR_5_BLOCKS + UTR_3_BLOCKS:

    ax.add_patch(
        patches.Rectangle(
            (start, -0.075),
            end - start,
            0.15,
            facecolor="#E0E0E0",
            edgecolor="#616161",
            linewidth=1,
            zorder=3
        )
    )

# Coding exons
for start, end in CODING_EXONS:

    ax.add_patch(
        patches.Rectangle(
            (start, -0.10),
            end - start,
            0.20,
            facecolor="#424242",
            edgecolor="#212121",
            linewidth=1.2,
            zorder=4
        )
    )


# ___ LOLLIPOP SETTINGS ___

spacing_y = 0.11
circle_size = 30
first_dot_y = 0.32
lollipop_base = first_dot_y - spacing_y / 2

plotted_variants = [
    variant
    for variant in variants
    if sum(variant["cohorts"].values()) > 0
]

class_names = {
    "ins": "Ins",
    "del": "Del",
    "seq. alter": "Seq",
    "snv": "SNV",
}


# ___ EXPORT VARIANT KEY ___

variant_key = []

for number, variant in enumerate(
    plotted_variants,
    start=1
):

    variant_key.append({
        "Number": number,
        "Genomic position": variant["pos"],
        "Variant class": variant["class"],
        "Healthy cohort (N0)": variant["cohorts"]["N0"],
        "Single-cancer cohort (N1)": variant["cohorts"]["N1"],
        "Multiple-cancer cohort (N2+)": variant["cohorts"]["N2"],
    })

variant_key_df = pd.DataFrame(variant_key)

variant_key_df.to_excel(
    VARIANT_KEY_FILE,
    index=False
)


# ___ DRAW LOLLIPOPS ___

for number, variant in enumerate(
    plotted_variants,
    start=1
):

    x_genomic = variant["pos"]
    x_spaced = variant["pos_spaced"]

    n0 = variant["cohorts"]["N0"]
    n1 = variant["cohorts"]["N1"]
    n2 = variant["cohorts"]["N2"]

    total = n0 + n1 + n2

    top_y = first_dot_y + (total - 1) * spacing_y

    # Connect original genomic position to displaced lollipop position.
    ax.plot(
        [x_genomic, x_genomic, x_spaced, x_spaced],
        [0, 0.18, lollipop_base, top_y],
        color="#9E9E9E",
        linewidth=0.9,
        zorder=1
    )

    # Stack one point per carrier.
    level = 0

    for cohort, count in (
        ("N0", n0),
        ("N1", n1),
        ("N2", n2)
    ):

        for _ in range(count):

            cy = first_dot_y + level * spacing_y

            ax.scatter(
                x_spaced,
                cy,
                s=circle_size,
                facecolor=cohort_colors[cohort],
                edgecolor="black",
                linewidths=0.6,
                zorder=5
            )

            level += 1

    # Add variant number and variant class above each lollipop.
    label_y = top_y + 0.18

    ax.plot(
        [x_spaced, x_spaced],
        [top_y + 0.03, label_y],
        color="#BDBDBD",
        linestyle=":",
        linewidth=0.75,
        zorder=1
    )

    variant_class_raw = str(
        variant["class"]
    ).strip().lower()

    variant_class = class_names.get(
        variant_class_raw,
        variant_class_raw
    )

    number_y = label_y + 0.02

    ax.text(
        x_spaced,
        number_y,
        str(number),
        rotation=0,
        ha="center",
        va="bottom",
        fontsize=10,
        fontweight="light",
        color="#BDBDBD",
        clip_on=True
    )

    ax.text(
        x_spaced,
        number_y + 0.30,
        variant_class,
        ha="center",
        va="bottom",
        rotation=-45,
        fontsize=9,
        fontweight="semibold",
        color="#BDBDBD",
        clip_on=True
    )


# ___ AXIS LIMITS ___
# Adjust plot height according to the maximum number of carriers.

max_total = max(
    (
        sum(variant["cohorts"].values())
        for variant in variants
    ),
    default=1
)

max_stack_y = (
    first_dot_y
    + (max_total - 1) * spacing_y
)

upper_y = max(
    5.5,
    max_stack_y + 1.0
)

lower_y = -0.75

ax.set_ylim(
    lower_y,
    upper_y
)

ax.set_xlim(
    PLOT_MIN_X,
    PLOT_MAX_X
)

ax.set_xticks(
    np.arange(
        PLOT_MIN_X,
        PLOT_MAX_X + 1,
        10000
    )
)


# ___ REGION LABELS ___

region_flank_y = -0.20
region_utr_y = -0.50

ax.text(
    PLOT_MIN_X + 1000,
    region_flank_y,
    LEFT_FLANK_LABEL,
    ha="left",
    va="center",
    fontsize=9,
    fontweight="semibold",
    color="#616161"
)

if UTR_5_BLOCKS:

    utr_5_positions = [
        coordinate
        for block in UTR_5_BLOCKS
        for coordinate in block
    ]

    ax.text(
        np.mean(utr_5_positions),
        region_utr_y,
        "5' UTR",
        ha="center",
        va="center",
        fontsize=9,
        fontweight="semibold",
        color="#616161"
    )

if UTR_3_BLOCKS:

    utr_3_positions = [
        coordinate
        for block in UTR_3_BLOCKS
        for coordinate in block
    ]

    ax.text(
        np.mean(utr_3_positions),
        region_utr_y,
        "3' UTR",
        ha="center",
        va="center",
        fontsize=9,
        fontweight="semibold",
        color="#616161"
    )

ax.text(
    PLOT_MAX_X - 1000,
    region_flank_y,
    RIGHT_FLANK_LABEL,
    ha="right",
    va="center",
    fontsize=9,
    fontweight="semibold",
    color="#616161"
)


# ___ AXIS LABELS ___

ax.xaxis.set_major_formatter(
    plt.FuncFormatter(
        lambda x, _: f"{int(x):,}"
    )
)

ax.set_xlabel(
    "Genomic location (bp)",
    fontsize=11,
    labelpad=15,
    fontweight="bold"
)

ax.set_ylabel(
    "Number of patients ($n$)",
    fontsize=11,
    labelpad=25,
    fontweight="bold",
    rotation=270
)

ax.set_yticks([])

ax.tick_params(
    axis="x",
    labelsize=10
)

for side in ("top", "right"):
    ax.spines[side].set_visible(False)

for side in ("left", "bottom"):
    ax.spines[side].set_color("#616161")


# ___ LEGEND ___

legend_elements = [
    Line2D(
        [0], [0],
        marker="o",
        color="w",
        label="= 1 Patient",
        markerfacecolor="none",
        markeredgecolor="black",
        markeredgewidth=0.6,
        markersize=6
    ),
    Line2D(
        [0], [0],
        marker="o",
        color="w",
        label="Healthy cohort (N0)",
        markerfacecolor=cohort_colors["N0"],
        markeredgecolor="black",
        markeredgewidth=0.6,
        markersize=6
    ),
    Line2D(
        [0], [0],
        marker="o",
        color="w",
        label="Single-cancer cohort (N1)",
        markerfacecolor=cohort_colors["N1"],
        markeredgecolor="black",
        markeredgewidth=0.6,
        markersize=6
    ),
    Line2D(
        [0], [0],
        marker="o",
        color="w",
        label="Multiple-cancer cohort (N2+)",
        markerfacecolor=cohort_colors["N2"],
        markeredgecolor="black",
        markeredgewidth=0.6,
        markersize=6
    ),
]

ax.legend(
    handles=legend_elements,
    loc="upper right",
    frameon=False,
    fontsize=9
)


# ___ SAVE PLOT ___

ax.grid(
    axis="x",
    linestyle=":",
    alpha=0.30,
    color="#CCCCCC"
)

ax.set_axisbelow(True)

fig.subplots_adjust(
    left=0.08,
    right=0.98,
    bottom=0.18,
    top=0.97
)

plt.savefig(
    PLOT_FILE,
    dpi=300,
    facecolor="white"
)

plt.show()

print(f"Done: {PLOT_FILE}")
print(f"Variant key: {VARIANT_KEY_FILE}")