from pycirclize import Circos
from pycirclize.parser import Genbank
from matplotlib.lines import Line2D
import warnings
# Suppress pycirclize warnings about plot range
warnings.filterwarnings("ignore", message="r_lim=.* is unexpected plot range.*")


def draw_plasmid_map_from_genbank_file(genbank_file_path, map_file_path, plasmid):

    # Read Genbank file
    gbk = Genbank(genbank_file_path)

    # Initialize Circos instance with genome size
    genome_size = gbk.range_size
    circos = Circos(sectors={gbk.name: genome_size})
    circos.text(f"{plasmid} plasmid", size=15, r=40)
    circos.rect(r_lim=(80, 85), fc="grey", ec="none", alpha=0.5)
    sector = circos.sectors[0]

    # Define category colors
    category_colors = {
        'Conjugation': '#4AA532',
        'Toxin-Antitoxin System': '#2067BF',
        'Origin of Replication': "#FF0066",
        'Origin of Transfer': '#E97451',
        'Plasmid Maintenance, Replication and Regulation': '#ED7E7E',
        'Metabolism': '#DBA602',
        'Stress Response': '#7859D3',
        'Other': '#6d4058',
        'Non-conjugative DNA mobility': 'black',
        'Antibiotic Resistance': 'green',
        'Metal and Biocide Resistance': 'red',
        'Open Reading Frame': 'skyblue',
        'Virulence and Defense Mechanism':'#85b90b'
    }

    # Function to add features to a track
    def add_features_to_track(track, features, default_color='blue', lw=0.5):
        for feat in features:
            category = feat.qualifiers.get('category', [None])[0]  # Adjust the 'category' key as necessary
            color = category_colors.get(category, default_color)  # Use a default color if category not in mapping
            track.genomic_features(feat, plotstyle="arrow", fc=color, lw=lw)
    
    # Safely extract features, handling empty strands
    def safe_extract_features(feature_type, strand):
        try:
            if gbk.get_seqid2features(feature_type, target_strand=strand):
                return gbk.extract_features(feature_type, target_strand=strand)
            return []
        except IndexError:
            return []  # Return empty list if there are no features

    # Extract features with checks for existence
    f_cds_feats = safe_extract_features("CDS", 1)
    f_oriC_feats = safe_extract_features("oriC", 1)
    f_oriT_feats = safe_extract_features("oriT", 1)
    f_all_feats = f_cds_feats + f_oriC_feats + f_oriT_feats

    r_cds_feats = safe_extract_features("CDS", -1)
    r_oriC_feats = safe_extract_features("oriC", -1)
    r_oriT_feats = safe_extract_features("oriT", -1)
    r_all_feats = r_cds_feats + r_oriC_feats + r_oriT_feats

    # Plot forward strand CDS with color based on category
    f_cds_track = sector.add_track((80, 85))
    if f_all_feats:
        add_features_to_track(f_cds_track, f_all_feats)

    # Plot reverse strand CDS only if features exist
    r_cds_track = sector.add_track((75, 80))
    if r_all_feats:
        add_features_to_track(r_cds_track, r_all_feats)

    # Plot 'gene' qualifier label if exists
    labels, label_pos_list = [], []
    all_cds_feats = safe_extract_features("CDS", None)  # Extract all CDS features regardless of strand
    for feat in all_cds_feats:
        start = int(str(feat.location.start))
        end = int(str(feat.location.end))
        label_pos = (start + end) / 2
        gene_name = feat.qualifiers.get("gene", [None])[0]
        if gene_name is not None:
            labels.append(gene_name)
            label_pos_list.append(label_pos)

    # Mobile element track
    if gbk.get_seqid2features("MGE", target_strand=1):
        t_mobile_track = sector.add_track((100, 105))
        t_mobile_feats = safe_extract_features("MGE", 1)
        t_mobile_track.genomic_features(t_mobile_feats, plotstyle="arrow", fc="yellow", lw=2)

        tlabels, tlabel_pos_list = [], []
        for feat in safe_extract_features("MGE", None):
            start = int(str(feat.location.start))
            end = int(str(feat.location.end))
            tlabel_pos = (start + end) / 2
            gene_name = feat.qualifiers.get("gene", [None])[0]
            product_name = feat.qualifiers.get("product", [None])[0]
            final_tn_name = (gene_name or "")
            if gene_name or product_name:
                tlabels.append(final_tn_name)
                tlabel_pos_list.append(tlabel_pos)

        t_mobile_track.xticks(tlabel_pos_list, tlabels, label_size=7.5, label_margin=1.5, label_orientation="horizontal")

    # Add labels if any exist
    if labels and label_pos_list:
        # Adjust the label margins and tick lengths
        labels_group1 = [label for i, label in enumerate(labels) if i % 2 == 0]
        label_pos_group1 = [pos for i, pos in enumerate(label_pos_list) if i % 2 == 0]

        labels_group2 = [label for i, label in enumerate(labels) if i % 2 != 0]
        label_pos_group2 = [pos for i, pos in enumerate(label_pos_list) if i % 2 != 0]

        # Plot the first group with one margin
        if labels_group1 and label_pos_group1:
            f_cds_track.xticks(
                label_pos_group1,
                labels_group1,
                label_size=6.5,
                label_margin=2.0,  # Margin for the first group
                label_orientation="vertical"
            )

        # Plot the second group with a different margin
        if labels_group2 and label_pos_group2:
            f_cds_track.xticks(
                label_pos_group2,
                labels_group2,
                label_size=6.5,
                label_margin=2.0,  # Margin for the second group
                label_orientation="vertical"
            )

    # Set tick intervals based on genome size
    if genome_size <= 100000:  # <= 100 kbp
        major_ticks_interval = 10000  # 10 kbp
    elif 100000 < genome_size <= 150000:  # 100 kbp - 150 kbp
        major_ticks_interval = 15000  # 15 kbp
    elif 150000 < genome_size <= 200000:  # 150 kbp - 200 kbp
        major_ticks_interval = 20000  # 20 kbp
    elif 200000 < genome_size <= 250000:   # > 200 kbp
        major_ticks_interval = 25000 
    else:
        major_ticks_interval = 30000  # 30 kbp

    minor_ticks_interval = major_ticks_interval / 5

    outer_track = sector.add_track((75, 85))
    outer_track.axis(fc="lightgrey")

    def skip_zero_label(value):
        if value == 0:
            return ""
        return f"{value / 1000:.1f} kb"

    outer_track.xticks_by_interval(
        major_ticks_interval,
        outer=False,
        label_formatter=skip_zero_label
    )
    outer_track.xticks_by_interval(
        minor_ticks_interval,
        outer=False,
        tick_length=1,
        show_label=False
    )

    # Plot replicon features if they exist
    rep_feats = safe_extract_features("Replicon", 1)
    if rep_feats:  # Check if replicon features are present before proceeding
        rep_track = sector.add_track((62, 67))
        rep_track.genomic_features(rep_feats, plotstyle="arrow", fc="yellow", lw=2)
        
        rlabels, rlabel_pos_list = [], []
        for feat in rep_feats:
            start = int(str(feat.location.start))
            end = int(str(feat.location.end))
            rlabel_pos = (start + end) / 2
            gene_name = feat.qualifiers.get("gene", [None])[0]
            product_name = feat.qualifiers.get("product", [None])[0]
            final_rp_name = (gene_name or "") 
            if gene_name or product_name:
                rlabels.append(final_rp_name)
                rlabel_pos_list.append(rlabel_pos)

        if rlabels and rlabel_pos_list:
            rep_track.xticks(rlabel_pos_list, rlabels, label_size=7.5, label_margin=-12, label_orientation="horizontal")

    fig = circos.plotfig()

    # Create legend
    line_handles = [
        Line2D([], [], color=color, label=category, lw=4)
        for category, color in category_colors.items()
    ]
    line_legend = circos.ax.legend(
        handles=line_handles,
        bbox_to_anchor=(0.275, 0.65),
        fontsize=8.5,
        handlelength=3,
    )

    fig.savefig(map_file_path, dpi=300)