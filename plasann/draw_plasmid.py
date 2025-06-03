from pycirclize import Circos
from pycirclize.parser import Genbank
from matplotlib.lines import Line2D
import warnings


def draw_plasmid_map_from_genbank_file(genbank_file_path, map_file_path, plasmid):
    try:
        # Suppress pycirclize warnings
        warnings.filterwarnings("ignore", category=UserWarning, module="pycirclize")
        
        # Read Genbank file
        gbk = Genbank(genbank_file_path)

        # Initialize Circos instance with genome size
        genome_size = gbk.range_size  # Make sure this is correctly retrieving the total genome size
        
        # Safety check for genome size
        if genome_size <= 0:
            print(f"⚠️ Invalid genome size: {genome_size}")
            warnings.resetwarnings()
            return False
            
        circos = Circos(sectors={gbk.name: genome_size})
        
        # Safety check for sectors
        if len(circos.sectors) == 0:
            print("⚠️ No sectors created in Circos")
            warnings.resetwarnings()
            return False
            
        circos.text(f"{plasmid}", size=15, r=0)
        circos.rect(r_lim=(70, 75), fc="grey", ec="none", alpha=0.5)
        sector = circos.sectors[0]

        # Define category colors
        category_colors = {
            'Conjugation': '#4AA532',
            'Toxin-Antitoxin System': '#2067BF',
            'Origin of Replication': "#FF0066",
            'Origin of Transfer': 'skyblue',
            'Plasmid Maintenance, Replication and Regulation': '#ED7E7E',
            'Metabolism': '#DBA602',
            'Stress Response': '#7859D3',
            'Other': '#6d4058',
            'Non-conjugative DNA mobility': 'black',
            'Antibiotic Resistance': 'green',
            'Metal and Biocide Resistance': 'red',
            'Open reading frame': 'blue',
            'Virulence and Defense Mechanism':'#85b90b',
            'Non coding RNA/Regulatory elements': 'purple',
            'Replicon': 'orange'
        }

        # Function to safely get qualifier value
        def safe_get_qualifier(feat, key, default=None):
            """Safely get qualifier value, handling missing keys and empty lists"""
            try:
                values = feat.qualifiers.get(key, [])
                if values and len(values) > 0:
                    return values[0]
                return default
            except (IndexError, AttributeError):
                return default

        # Function to add features to a track
        def add_features_to_track(track, features, default_color='blue', lw=0.5):
            for feat in features:
                category = safe_get_qualifier(feat, 'category')
                color = category_colors.get(category, default_color)  # Use a default color if category not in mapping
                track.genomic_features(feat, plotstyle="arrow", fc=color, lw=lw)

        # Safe feature extraction function
        def safe_extract_features(gbk, feature_type, target_strand=None):
            """Safely extract features, return empty list if none found"""
            try:
                if target_strand is not None:
                    seqid2features = gbk.get_seqid2features(feature_type, target_strand)
                    if seqid2features:  # Check if any features found
                        return gbk.extract_features(feature_type, target_strand=target_strand)
                else:
                    # No target strand specified
                    seqid2features = gbk.get_seqid2features(feature_type)
                    if seqid2features:
                        return gbk.extract_features(feature_type)
                return []
            except (IndexError, KeyError, AttributeError) as e:
                print(f"⚠️ Warning: Could not extract {feature_type} features: {e}")
                return []

        # Extract features safely - keep origins separate from CDS
        f_cds_feats = safe_extract_features(gbk, "CDS", target_strand=1)
        f_all_feats = f_cds_feats  # Only CDS for forward track

        r_cds_feats = safe_extract_features(gbk, "CDS", target_strand=-1)
        r_all_feats = r_cds_feats  # Only CDS for reverse track

        # Only create tracks if we have features
        if f_all_feats:
            # Plot forward strand CDS with color based on category
            f_cds_track = sector.add_track((70, 75))
            add_features_to_track(f_cds_track, f_all_feats)
        else:
            f_cds_track = None

        if r_all_feats:
            # Repeat for reverse strand CDS
            r_cds_track = sector.add_track((65, 70))
            add_features_to_track(r_cds_track, r_all_feats)
        else:
            r_cds_track = None

        # Plot 'gene' qualifier label if exists
        labels, label_pos_list = [], []
        all_cds_feats = safe_extract_features(gbk, "CDS")
        
        for feat in all_cds_feats:
            try:
                start = int(str(feat.location.start))
                end = int(str(feat.location.end))
                label_pos = (start + end) / 2
                gene_name = safe_get_qualifier(feat, "gene")
                if gene_name and gene_name not in ["ORF"]:
                    labels.append(gene_name)
                    label_pos_list.append(label_pos)
            except (ValueError, AttributeError) as e:
                print(f"⚠️ Warning: Could not process CDS feature: {e}")
                continue

        # Mobile element track
        mge_feats = safe_extract_features(gbk, "MGE", target_strand=1)
        if mge_feats:
            t_mobile_track = sector.add_track((100, 105))
            t_mobile_track.genomic_features(mge_feats, plotstyle="arrow", fc="yellow", lw=2)

            tlabels, tlabel_pos_list = [], []
            all_mge_feats = safe_extract_features(gbk, "MGE")
            for feat in all_mge_feats:
                try:
                    start = int(str(feat.location.start))
                    end = int(str(feat.location.end))
                    tlabel_pos = (start + end) / 2
                    gene_name = safe_get_qualifier(feat, "gene")
                    if gene_name:
                        tlabels.append(gene_name)
                        tlabel_pos_list.append(tlabel_pos)
                except (ValueError, AttributeError):
                    continue

            if tlabels:
                t_mobile_track.xticks(tlabel_pos_list, tlabels, label_size=10, label_margin=1, label_orientation="vertical")

        # Origin track - with category-based colors like CDS
        oriv_feats = safe_extract_features(gbk, "ORIV")
        orit_feats = safe_extract_features(gbk, "ORIT")
        ori_feats = oriv_feats + orit_feats
        
        if ori_feats:
            ori_track = sector.add_track((95, 100))
            
            # Color each feature based on its category (like CDS)
            for feat in ori_feats:
                category = safe_get_qualifier(feat, 'category')
                color = category_colors.get(category, 'red')
                ori_track.genomic_features(feat, plotstyle="arrow", fc=color, lw=2)

            ori_labels, ori_label_pos_list = [], []
            for feat in ori_feats:
                try:
                    start = int(str(feat.location.start))
                    end = int(str(feat.location.end))
                    ori_label_pos = (start + end) / 2
                    gene_name = safe_get_qualifier(feat, "gene")
                    if gene_name:
                        ori_labels.append(gene_name)
                        ori_label_pos_list.append(ori_label_pos)
                except (ValueError, AttributeError):
                    continue

            if ori_labels:
                ori_track.xticks(ori_label_pos_list, ori_labels, label_size=10, label_margin=1, label_orientation="vertical")

        # Only add CDS labels if we have tracks and labels
        if labels and label_pos_list:
            # Adjust the label margins and tick lengths for CDS
            labels_group1 = [label for i, label in enumerate(labels) if i % 2 == 0]
            label_pos_group1 = [pos for i, pos in enumerate(label_pos_list) if i % 2 == 0]

            labels_group2 = [label for i, label in enumerate(labels) if i % 2 != 0]
            label_pos_group2 = [pos for i, pos in enumerate(label_pos_list) if i % 2 != 0]

            # Plot individual labels in group 1 (forward track) with alternating margins/tick lengths
            if f_cds_track and labels_group1:
                for i, (label, pos) in enumerate(zip(labels_group1, label_pos_group1)):
                    margin = 0.5 if i % 2 == 0 else 1
                    tick_len = 2 if i % 2 == 0 else 5
                    f_cds_track.xticks(
                        [pos], [label],
                        tick_length=tick_len,
                        label_size=10,
                        label_margin=margin,
                        label_orientation="vertical"
                    )
            
            # Plot individual labels in group 2 (reverse track) with alternating margins/tick lengths
            if r_cds_track and labels_group2:
                for i, (label, pos) in enumerate(zip(labels_group2, label_pos_group2)):
                    margin = 0.5 if i % 2 == 0 else 1
                    tick_len = 2 if i % 2 == 0 else 12
                    r_cds_track.xticks(
                        [pos], [label],
                        outer=False,
                        tick_length=tick_len,
                        label_size=10,
                        label_margin=margin,
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
        else :
            major_ticks_interval= 30000 # 30 kbp

        minor_ticks_interval = major_ticks_interval / 5

        outer_track = sector.add_track((65, 75))
        outer_track.axis(fc="lightgrey")

        tick_track = sector.add_track((30, 28))
        tick_track.axis(fc="black")

        def skip_zero_label(value):
            if value == 0:
                return ""
            return f"{value / 1000:.1f} kb"

        tick_track.xticks_by_interval(
            major_ticks_interval,
            outer=False,
            label_formatter=skip_zero_label
        )
        tick_track.xticks_by_interval(
            minor_ticks_interval,
            outer=False,
            tick_length=1,
            show_label=False
        )

        # Add inner track for replicon feature
        rep_feats = safe_extract_features(gbk, "Replicon", target_strand=1)
        if rep_feats:  # Check if replicon features are present before proceeding
            rep_track = sector.add_track((97, 102))
            rep_track.genomic_features(rep_feats, plotstyle="box", fc="yellow", lw=2)
            
            rlabels, rlabel_pos_list = [], []
            for feat in rep_feats:
                try:
                    start = int(str(feat.location.start))
                    end = int(str(feat.location.end))
                    rlabel_pos = (start + end) / 2
                    gene_name = safe_get_qualifier(feat, "gene")
                    if gene_name:
                        rlabels.append(gene_name)
                        rlabel_pos_list.append(rlabel_pos)
                except (ValueError, AttributeError):
                    continue

            if rlabels:
                rep_track.xticks(rlabel_pos_list, rlabels, tick_length=0, label_size=8.5, label_margin=-10, label_orientation="horizontal")

        fig = circos.plotfig()

        # Legend
        legend_lines = [
            Line2D([], [], color=color, lw=4, label=label)
            for label, color in category_colors.items()
        ]
        fig.legend(
            handles=legend_lines,
            bbox_to_anchor=(1.05, 0.5),
            loc="center left",
            borderaxespad=0,
            fontsize=14,
            title="Feature Categories",
            title_fontsize=16
        )

        fig.savefig(map_file_path, dpi=300)
        print(f"✅ Plasmid map saved: {map_file_path}")
        
        # Reset warnings
        warnings.resetwarnings()
        
        return True
        
    except Exception as e:
        print(f"❌ Error creating plasmid map: {e}")
        print(f"   GenBank file: {genbank_file_path}")
        print(f"   Plasmid: {plasmid}")
        
        # Reset warnings even on error
        warnings.resetwarnings()
        
        return False