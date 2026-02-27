#!/usr/bin/env python

import sys
import os
from collections import defaultdict
from datetime import datetime # For dynamic timestamp

# A list of predefined colors to cycle through for the rules in colors.conf
RULE_COLORS = [
    '220,20,60', '218,112,214', '148,0,211', '65,105,225', '0,191,255', '72,209,204', '127,255,170', '173,255,47',
    '255,255,0', '255,215,0', 'dred', 'dblue', 'dgreen', 'dpurple', 'dorange',
    'lred', 'lblue', 'lgreen', 'lpurple', 'lorange'
]

def find_sensible_prefix(chromosome_names_list):
    """
    Finds a sensible common prefix for a list of chromosome names.
    Prioritizes prefixes ending with '_' or '-'.
    """
    if not chromosome_names_list:
        return ""

    # Use os.path.commonprefix to find the raw longest common prefix
    # This function requires a list of strings.
    # If chromosome_names_list is a set, convert it to a list first.
    str_list = list(chromosome_names_list)
    if not str_list: # Handle empty list after conversion (e.g. if input was empty set)
        return ""

    lcp = os.path.commonprefix(str_list)
    if not lcp:
        return ""

    # Attempt 1: Check if LCP can be trimmed to end with '_' and still be a common prefix
    if '_' in lcp:
        candidate_prefix = lcp.rsplit('_', 1)[0] + '_'
        if all(s.startswith(candidate_prefix) for s in str_list):
            return candidate_prefix
            
    # Attempt 2: Check if LCP can be trimmed to end with '-' and still be a common prefix
    if '-' in lcp:
        candidate_prefix = lcp.rsplit('-', 1)[0] + '-'
        if all(s.startswith(candidate_prefix) for s in str_list):
            return candidate_prefix
            
    # Attempt 3: If the LCP itself is what all strings start with (e.g. "chr" for "chr1", "chr2")
    # and it's of reasonable length (e.g. >=2 chars). This is already handled by os.path.commonprefix.
    # Avoid overly short prefixes unless they are the only option.
    if len(lcp) >= 2: # Arbitrary minimum length for a "sensible" raw prefix
        return lcp
        
    return "" # Return empty if no suitable prefix found

def parse_synteny_for_ordering_and_genomes(input_synteny_file):
    """
    Parses the synteny file for chromosome ordering, genome sets, and link rules.
    Returns: (
        final_ordered_chromosome_list,
        col1_chrs_for_link_rules (Genome A for links),
        all_genome_A_chrs_set (all unique names from col1),
        all_genome_B_chrs_set (all unique names from col4),
        error_message
    )
    """
    all_col1_chrs_set = set()
    all_col4_chrs_set = set()
    link_counts = defaultdict(lambda: defaultdict(int))

    try:
        with open(input_synteny_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) < 4:
                    return None, None, None, None, f"Error: Line {line_num} in synteny file '{input_synteny_file}' has fewer than 4 columns."
                
                a_chr, b_chr = parts[0], parts[3]
                all_col1_chrs_set.add(a_chr)
                all_col4_chrs_set.add(b_chr)
                link_counts[a_chr][b_chr] += 1
                
    except FileNotFoundError:
        return None, None, None, None, f"Error: Input synteny file '{input_synteny_file}' not found."
    except Exception as e:
        return None, None, None, None, f"Error reading input synteny file '{input_synteny_file}': {e}"

    if not all_col1_chrs_set: # Genome A must have chromosomes for link rules
        return None, None, None, None, f"No data found in column 1 of '{input_synteny_file}' or file is unparsable."

    col1_chrs_for_link_rules = sorted(list(all_col1_chrs_set))

    genome_A_initial_order = list(col1_chrs_for_link_rules)
    genome_B_partners_for_A = []
    
    for a_chr in genome_A_initial_order:
        best_b_partner = None
        max_links = -1
        for b_chr_candidate in sorted(list(all_col4_chrs_set)): 
            current_links = link_counts[a_chr].get(b_chr_candidate, 0)
            if current_links > max_links:
                max_links = current_links
                best_b_partner = b_chr_candidate
        
        if best_b_partner is not None and max_links > 0:
            genome_B_partners_for_A.append(best_b_partner)
        else:
            genome_B_partners_for_A.append(None)

    genome_B_partner_segment_ordered = [p for p in genome_B_partners_for_A if p is not None]
    seen_b_partners = set()
    unique_genome_B_partner_segment = []
    for b_partner in genome_B_partner_segment_ordered:
        if b_partner not in seen_b_partners:
            unique_genome_B_partner_segment.append(b_partner)
            seen_b_partners.add(b_partner)
    genome_B_partner_segment_reversed = unique_genome_B_partner_segment[::-1]
    
    remaining_genome_B_chrs = sorted(list(all_col4_chrs_set - seen_b_partners))
    
    final_ordered_list_raw = genome_A_initial_order + genome_B_partner_segment_reversed + remaining_genome_B_chrs
    seen_final = set()
    final_ordered_chromosome_list = []
    for item in final_ordered_list_raw:
        if item not in seen_final:
            final_ordered_chromosome_list.append(item)
            seen_final.add(item)
            
    return final_ordered_chromosome_list, col1_chrs_for_link_rules, all_col1_chrs_set, all_col4_chrs_set, None

def generate_links_conf_content(col1_chrs_for_link_rules, synteny_file_basename):
    output_lines = ["<links>", ""]
    output_lines.append("<link>")
    output_lines.append(f"file    = {synteny_file_basename}")
    output_lines.append("radius  = 0.99r")
    output_lines.append("ribbon  = yes")
    output_lines.append("color   = grey")
    output_lines.append("thickness = 2p")
    output_lines.append("<rules>")
    for i, entry in enumerate(col1_chrs_for_link_rules):
        output_lines.append("<rule>")
        output_lines.append(f'condition = var(chr1) eq "{entry}"')
        output_lines.append(f"color     = rule_color_{i+1}")
        output_lines.append("</rule>")
    output_lines.append("</rules>")
    output_lines.append("</link>")
    output_lines.append("</links>")
    return "\n".join(output_lines)

def generate_colors_conf_content(col1_chrs_for_link_rules):
    output_lines = ["# Colors for synteny links"]
    for i, entry in enumerate(col1_chrs_for_link_rules):
        color_to_use = RULE_COLORS[i % len(RULE_COLORS)]
        output_lines.append(f"rule_color_{i+1} = {color_to_use} # For links from {entry}")
    output_lines.append("# <<include etc/colors.conf>>")
    output_lines.append("# <<include etc/brewer.conf>>")
    return "\n".join(output_lines)

def generate_circos_conf_content(final_ordered_chromosome_list, genome_A_chrs_set, genome_B_chrs_set, current_user, current_datetime_utc):
    chrom_order_str = ",".join(final_ordered_chromosome_list) if final_ordered_chromosome_list else ""
    
    order_explanation = """
# --- chromosomes_order Explanation ---
# Generated based on synteny: Genome A (synteny col1) first, then primary Genome B partners
# (synteny col4) in reverse order, then remaining Genome B.
# --- IMPORTANT: Review and MANUALLY EDIT 'chromosomes_order' if this doesn't match your desired layout! ---
"""
    chrom_order_conf_str = f"chromosomes_order = {chrom_order_str}" if chrom_order_str else ""

    prefix_A = find_sensible_prefix(genome_A_chrs_set)
    # For prefix_B, consider only chromosomes from col4 not already in col1 to find a distinct prefix for "pure B"
    # However, user wants to scale based on all chromosomes from col4 as "Genome B"
    prefix_B = find_sensible_prefix(genome_B_chrs_set) 

    scale_config_lines = []
    scale_comment = "# Chromosome scaling by auto-detected prefixes for Genome A (synteny col1) and Genome B (synteny col4)."
    scale_value = "0.47rn" # Target proportion for each genome half

    scale_parts = []
    if prefix_A:
        scale_parts.append(f"/{prefix_A}/={scale_value}")
    else:
        scale_comment += "\n# Could not determine a common prefix for Genome A chromosomes (from synteny col1)."
        scale_comment += "\n# You may need to define 'chromosomes_scale' manually for Genome A."

    if prefix_B:
        if prefix_A and prefix_A == prefix_B:
            scale_comment += f"\n# WARNING: Auto-detected prefix for Genome A ('{prefix_A}') and Genome B ('{prefix_B}') are IDENTICAL."
            scale_comment += "\n# This prefix-based scaling will not separate them effectively. Consider manual 'chromosomes_scale',"
            scale_comment += "\n# or ensure your karyotype files use distinct names/prefixes for the two genomes."
            # Do not add prefix_B to scale_parts if it's identical and problematic
        elif prefix_A and prefix_A != prefix_B :
             scale_parts.append(f"/{prefix_B}/={scale_value}")
        elif not prefix_A: # prefix_A is empty, prefix_B might be valid
             scale_parts.append(f"/{prefix_B}/={scale_value}")

    else: # prefix_B is empty
        scale_comment += "\n# Could not determine a common prefix for Genome B chromosomes (from synteny col4)."
        scale_comment += "\n# You may need to define 'chromosomes_scale' manually for Genome B."

    if scale_parts:
        scale_config_lines.append(scale_comment)
        scale_config_lines.append(f"chromosomes_scale = {','.join(scale_parts)}")
        scale_config_lines.append(f"# Each rn value (e.g., {scale_value}) is a proportion of the ideogram circle.")
    else:
        scale_config_lines.append(scale_comment)
        scale_config_lines.append("# 'chromosomes_scale' was not automatically generated due to missing or conflicting prefixes.")
        scale_config_lines.append("# Please define it manually if needed, e.g.: /prefix1_/=0.47rn,/prefix2_/=0.47rn")
    
    chromosomes_scale_str = "\n".join(scale_config_lines)

    content = f"""
# Circos Configuration File (Generated by script)
# Generated on: {current_datetime_utc} by user: {current_user}

karyotype = karyotype_sp1.txt,karyotype_sp2.txt

chromosomes_units = 1000000

{order_explanation.strip()}
{chrom_order_conf_str}

{chromosomes_scale_str}

angle_offset* = -80 

<ideogram>
<<include ideogram.conf>>
</ideogram>

<colors>
<<include colors.conf>>
</colors>

<<include etc/colors_fonts_patterns.conf>> 
<<include ticks.conf>>
<<include links.conf>>

<image>
<<include etc/image.conf>> 
</image>

<<include etc/housekeeping.conf>> 
anti_aliasing* = no 
"""
    return content.strip()

def generate_ideogram_conf_content(final_ordered_chromosome_list, genome_A_chrs_set):
    pairwise_spacing_str = ""
    if len(final_ordered_chromosome_list) > 1:
        last_A_in_order = None
        idx_last_A = -1
        for i, chrom in enumerate(final_ordered_chromosome_list):
            if chrom in genome_A_chrs_set: # Check if current chrom is considered Genome A
                idx_last_A = i
            else: 
                if idx_last_A != -1: 
                    break 
        
        if idx_last_A != -1 and idx_last_A < len(final_ordered_chromosome_list) -1 :
            last_A_in_order = final_ordered_chromosome_list[idx_last_A]
            first_B_after_A = final_ordered_chromosome_list[idx_last_A+1]
            pairwise_spacing_str += f"\n    <pairwise {last_A_in_order};{first_B_after_A}>\n    spacing = 25u\n    </pairwise>"

        first_A_in_order = final_ordered_chromosome_list[0]
        last_chr_in_order = final_ordered_chromosome_list[-1]
        if last_chr_in_order not in genome_A_chrs_set and last_chr_in_order != first_A_in_order:
             pairwise_spacing_str += f"\n    <pairwise {last_chr_in_order};{first_A_in_order}>\n    spacing = 25u\n    </pairwise>"

    return f"""
# Ideogram Configuration
<spacing>
default = 5u{pairwise_spacing_str}
</spacing>
radius           = 0.90r 
thickness        = 20p
fill             = yes 
stroke_color     = dgrey
stroke_thickness = 2p
show_label       = yes
label_font       = default 
label_radius     = 1r+75p 
label_size       = 15     
label_parallel   = yes
"""

def generate_ticks_conf_content():
    return """
# Tick Mark Configuration
show_ticks          = yes
show_tick_labels    = yes
<ticks>
radius           = 1r     
color            = black
thickness        = 2p
font             = default 
size             = 20p     
multiplier       = 1e-6
format           = %d
<tick>
spacing        = 5u 
size           = 10p 
color          = grey
</tick>
<tick>
spacing        = 25u 
size           = 15p
color          = black
show_label     = yes
label_size     = 12p 
label_offset   = 10p 
format         = %dMb 
</tick>
</ticks>
"""

def write_file(filename, content, output_dir="."):
    filepath = os.path.join(output_dir, filename)
    try:
        os.makedirs(output_dir, exist_ok=True)
        with open(filepath, "w") as outfile:
            outfile.write(content)
        print(f"Successfully generated '{filepath}'")
    except IOError as e:
        print(f"Error writing to file '{filepath}': {e}", file=sys.stderr)
        return False
    return True

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python generate_circos_configs.py <path_to_your_synteny_data_file>", file=sys.stderr)
        sys.exit(1)

    input_synteny_file = sys.argv[1]
    synteny_file_basename = os.path.basename(input_synteny_file)
    
    current_datetime_utc = datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC")
    try:
        current_user = os.getlogin()
    except Exception: 
        current_user = "N/A" # Fallback if os.getlogin() fails

    ordered_chrom_list, col1_chrs_links, genome_A_chrs, genome_B_chrs, error_msg = \
        parse_synteny_for_ordering_and_genomes(input_synteny_file)

    if error_msg:
        print(error_msg, file=sys.stderr)
        sys.exit(1)

    output_directory = "circos_config"

    links_content = generate_links_conf_content(col1_chrs_links, synteny_file_basename)
    colors_content = generate_colors_conf_content(col1_chrs_links)
    circos_main_content = generate_circos_conf_content(ordered_chrom_list, genome_A_chrs, genome_B_chrs, current_user, current_datetime_utc) 
    ideogram_content = generate_ideogram_conf_content(ordered_chrom_list, genome_A_chrs) 
    ticks_content = generate_ticks_conf_content()

    all_files_written = True
    if not write_file("links.conf", links_content, output_directory): all_files_written = False
    if not write_file("colors.conf", colors_content, output_directory): all_files_written = False
    if not write_file("circos.conf", circos_main_content, output_directory): all_files_written = False
    if not write_file("ideogram.conf", ideogram_content, output_directory): all_files_written = False
    if not write_file("ticks.conf", ticks_content, output_directory): all_files_written = False

    if all_files_written:
        print(f"\n--- Configuration files generated and backup in '{output_directory}/' directory ---")
        print(f"--- User: {current_user}, Timestamp: {current_datetime_utc} ---")
        print("\n--- IMPORTANT NEXT STEPS ---")
        print(f"1. CHROMOSOMES ORDER: Review 'chromosomes_order' in 'circos.conf' ")
        print(f"2. CHROMOSOMES SCALE: Open 'circos.conf'. **Review the auto-generated 'chromosomes_scale' section.**")
        print("   - Verify if the detected prefixes (e.g., /prefixA/) are correct for your genomes.")
        print(f"   - Adjust the '0.47rn' values if needed to balance visual proportions.")
        print("   - If prefixes were not detected correctly or are identical, you'll need to define 'chromosomes_scale' manually.")
        print(f"3. GENOME GAPS (SPACING): Review <pairwise> spacing in 'ideogram.conf'. Adjust '25u' if needed.")
        print(f"4. SYNTENY DATA FILE: update path in links.conf.")
        print(f"5. ETC FILES: Ensure Circos can find its 'etc' directory.")
        print(f"6. CUSTOMIZE: Review all '.conf' files for further adjustments.")
        print("\nTo run Circos (navigate into the output directory):")
        print("  circos -conf circos.conf")
    else:
        print("\nErrors occurred while writing one or more configuration files.", file=sys.stderr)
        sys.exit(1)
