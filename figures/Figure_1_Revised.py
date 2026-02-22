"""
Figure 1 Revised — ASV Analysis Workflow
Addresses Reviewer 2, Comment #8: Improved formatting, matched M&M section numbers.

Changes from original:
- Box numbers now match Methods sections 2.1–2.6
- Increased text size (minimum 9pt)
- Compact legend as horizontal bar at bottom
- Fixed text alignment
- Clearer arrow flow
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig, ax = plt.subplots(1, 1, figsize=(12, 14))
ax.set_xlim(0, 12)
ax.set_ylim(0, 16)
ax.axis('off')

# Color scheme
colors = {
    'specimen': '#4472C4',      # Blue - Specimen & Molecular
    'bioinformatics': '#70AD47', # Green - Bioinformatics
    'phylogenetic': '#FFC000',   # Gold - Phylogenetic Framework
    'classification': '#ED7D31',  # Orange - Classification
    'statistics': '#A5A5A5',     # Gray - Statistics
    'output': '#5B9BD5',         # Light blue - Final Output
}
text_colors = {
    'specimen': 'white',
    'bioinformatics': 'white',
    'phylogenetic': 'black',
    'classification': 'white',
    'statistics': 'white',
    'output': 'white',
}

def draw_box(ax, x, y, w, h, title, subtitle, color_key, section_num=None):
    """Draw a workflow box with section number."""
    box = FancyBboxPatch((x, y), w, h,
                         boxstyle="round,pad=0.1",
                         facecolor=colors[color_key],
                         edgecolor='black', linewidth=1.2)
    ax.add_patch(box)

    tc = text_colors[color_key]
    if section_num:
        ax.text(x + w/2, y + h - 0.25, f"Section {section_num}",
                ha='center', va='top', fontsize=8, fontstyle='italic',
                color=tc, alpha=0.85)
        ax.text(x + w/2, y + h/2 - 0.05, title,
                ha='center', va='center', fontsize=10, fontweight='bold',
                color=tc)
    else:
        ax.text(x + w/2, y + h/2 + 0.15, title,
                ha='center', va='center', fontsize=10, fontweight='bold',
                color=tc)

    ax.text(x + w/2, y + 0.25, subtitle,
            ha='center', va='bottom', fontsize=7.5,
            color=tc, alpha=0.9)

def draw_arrow(ax, x1, y1, x2, y2, label=None):
    """Draw an arrow with optional label."""
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color='#333333',
                               lw=1.5, connectionstyle='arc3,rad=0'))
    if label:
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        ax.text(mx + 0.15, my, label, fontsize=7, fontstyle='italic',
                color='#555555', ha='left', va='center')

# === Draw boxes (top to bottom) ===

# Box 2.1: Specimen Collection & Processing
draw_box(ax, 1, 13.5, 4.5, 1.5,
         "Specimen Collection\n& Processing",
         "Sampling, Sorting, Imaging, ID, DNA Extraction",
         'specimen', section_num="2.1")

# Box 2.2: Molecular Data Generation
draw_box(ax, 1, 11.2, 4.5, 1.5,
         "Molecular Data\nGeneration",
         "COI PCR Amplification, Illumina Sequencing",
         'specimen', section_num="2.2")

# Box 2.3: Bioinformatics Pipeline
draw_box(ax, 1, 8.9, 4.5, 1.5,
         "Bioinformatics\nPipeline",
         "Cutadapt, PEAR, VSEARCH, filtertranslate",
         'bioinformatics', section_num="2.3")

# Box 2.4: Phylogenetic Framework (right side)
draw_box(ax, 6.5, 8.9, 4.5, 1.5,
         "Phylogenetic Framework\nConstruction",
         "Mitogenome Assembly, MAFFT, FastTree",
         'phylogenetic', section_num="2.4")

# Box 2.5: ASV Authentication & Classification
draw_box(ax, 3.5, 6.3, 5, 1.7,
         "ASV Authentication\n& Classification",
         "MRCT Filtering, Phylo Placement, Taxonomic Congruence",
         'classification', section_num="2.5")

# Box 2.6: Statistical Analysis
draw_box(ax, 3.5, 4.0, 5, 1.5,
         "Statistical Analysis",
         "Chi-square, Cramér's V, Feature Importance (Python/scipy)",
         'statistics', section_num="2.6")

# Final Output box
draw_box(ax, 3.5, 1.8, 5, 1.5,
         "Final Outputs",
         "Authenticated Barcodes, Classifications, Statistics, Images",
         'output')

# === Draw arrows ===
# 2.1 → 2.2
draw_arrow(ax, 3.25, 13.5, 3.25, 12.7, "DNA Extract")
# 2.2 → 2.3
draw_arrow(ax, 3.25, 11.2, 3.25, 10.4, "Raw Paired-End Reads")
# 2.3 → 2.5
draw_arrow(ax, 3.25, 8.9, 4.8, 8.0, "ASV FASTA")
# 2.4 → 2.5
draw_arrow(ax, 8.75, 8.9, 7.2, 8.0, "Reference Tree")
# Metadata input to 2.5
ax.annotate('', xy=(6.0, 8.0), xytext=(6.0, 8.6),
            arrowprops=dict(arrowstyle='->', color='#333333',
                           lw=1.0, linestyle='dashed'))
ax.text(6.15, 8.35, "Metadata", fontsize=7, fontstyle='italic',
        color='#555555')
# 2.5 → 2.6
draw_arrow(ax, 6.0, 6.3, 6.0, 5.5, "Classification Table")
# 2.6 → Output
draw_arrow(ax, 6.0, 4.0, 6.0, 3.3, "Results")

# === Legend (compact horizontal bar at bottom) ===
legend_y = 0.6
legend_items = [
    ('Specimen & Molecular', 'specimen'),
    ('Bioinformatics', 'bioinformatics'),
    ('Phylogenetic Framework', 'phylogenetic'),
    ('Classification', 'classification'),
    ('Statistics', 'statistics'),
    ('Final Output', 'output'),
]

ax.text(6, 0.05, "Methodological Components", fontsize=8,
        ha='center', va='bottom', fontstyle='italic', color='#666666')

x_start = 0.5
for label, key in legend_items:
    box = FancyBboxPatch((x_start, legend_y), 0.3, 0.3,
                         boxstyle="round,pad=0.02",
                         facecolor=colors[key], edgecolor='black', linewidth=0.5)
    ax.add_patch(box)
    ax.text(x_start + 0.45, legend_y + 0.15, label, fontsize=7,
            va='center', ha='left')
    x_start += 2.0

# Title
ax.text(6, 15.7, "ASV Analysis Workflow", fontsize=14,
        ha='center', va='center', fontweight='bold')

plt.tight_layout()
plt.savefig('/Users/sarawut/Desktop/Manuscript_ASV_selection/figures/Figure_1_Revised.png',
            dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig('/Users/sarawut/Desktop/Manuscript_ASV_selection/figures/Figure_1_Revised.pdf',
            bbox_inches='tight', facecolor='white')
print("Figure 1 (Revised) saved successfully.")
plt.show()
