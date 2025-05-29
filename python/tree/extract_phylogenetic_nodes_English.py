import sys
import io
import os
import argparse
import re
from Bio import Phylo # Bio.Phylo is now only used if --replace is NOT specified

def get_names_from_phylotree( # Renamed to clarify it uses Bio.Phylo tree object
    start_clade, 
    node_type="all",
    include_unnamed_prefix=None,
    min_depth=None,
    max_depth=None,
    min_branch_length=None,
    max_branch_length=None,
    quiet=False
):
    """
    (Used when --replace is NOT active)
    From a Bio.Phylo.Clade object, extracts node names based on criteria.
    Depth and traversal are relative to start_clade.
    """
    if not quiet:
        print(f"INFO: (Bio.Phylo mode) Starting to extract node list from subtree rooted at '{start_clade.name or 'UnnamedSubtreeRoot'}'...", file=sys.stderr)

    node_names_extracted = []
    unnamed_counter = 0
    
    clade_depths = None
    if min_depth is not None or max_depth is not None:
        try:
            clade_depths = start_clade.depths(unit_branch_lengths=False) # Depths relative to start_clade
            if not quiet:
                start_clade_relative_depth = clade_depths.get(start_clade, 'N/A')
                print(f"INFO: (Bio.Phylo mode) Calculated node depths relative to the subtree root. Subtree root '{start_clade.name or 'UnnamedSubtreeRoot'}' relative depth is {start_clade_relative_depth}", file=sys.stderr)
        except Exception as e:
            if not quiet:
                print(f"WARNING: (Bio.Phylo mode) Error calculating node depths: {e}. Depth filtering will be ignored.", file=sys.stderr)
            clade_depths = None

    for clade in start_clade.find_clades(order='preorder'): # Iterate from start_clade
        current_name = clade.name
        
        if current_name is None and include_unnamed_prefix:
            current_name = f"{include_unnamed_prefix}{unnamed_counter}"
            unnamed_counter += 1
        
        if current_name is None: # If still no name (e.g., --include-unnamed not specified or it was already named None)
            continue

        # 1. Type filtering
        is_correct_type = False
        if node_type == "all": is_correct_type = True
        elif node_type == "internal":
            if not clade.is_terminal(): is_correct_type = True
        elif node_type == "leaf":
            if clade.is_terminal(): is_correct_type = True
        
        if not is_correct_type: continue

        # 2. Depth filtering
        if clade_depths:
            depth = clade_depths.get(clade) 
            if depth is not None: # Ensure clade is in depths dictionary
                if min_depth is not None and depth < min_depth: continue
                if max_depth is not None and depth > max_depth: continue

        # 3. Branch length filtering
        branch_len = clade.branch_length
        if min_branch_length is not None or max_branch_length is not None:
            if branch_len is None: continue # If filtering by branch length, skip nodes without it
            if min_branch_length is not None and branch_len < min_branch_length: continue
            if max_branch_length is not None and branch_len > max_branch_length: continue
        
        node_names_extracted.append(current_name)
        
    if not quiet:
        print(f"INFO: (Bio.Phylo mode) Initially extracted {len(node_names_extracted)} node names from the subtree based on basic criteria.", file=sys.stderr)
    return node_names_extracted

def main():
    parser = argparse.ArgumentParser(
        description="Extracts node names from a phylogenetic tree file or performs text replacement and outputs the modified file content.",
        formatter_class=argparse.RawTextHelpFormatter # Better display for multiline help
    )
    
    # --- Input Parameters ---
    parser.add_argument(
        "tree_input", 
        help="Required: File path containing tree information, or a direct tree string (if not found as a file path)."
    )
    parser.add_argument(
        "--encoding", 
        default="utf-8", 
        help="Encoding for reading and writing files.\nDefault: utf-8"
    )

    # --- Text Replacement Mode Parameter ---
    parser.add_argument(
        "--replace",
        nargs=2,
        action='append', # Allow multiple uses
        metavar=('OLD_NAME', 'NEW_NAME'),
        dest='replacements',
        default=[], # Ensure default is an empty list
        help=("Perform plain text replacement: replaces OLD_NAME with NEW_NAME in the file.\n"
              "This option can be used multiple times for several replacements.\n"
              "If this option is used, the script will output the modified full file content,\n"
              "and most Bio.Phylo related parameters below (like --format, --subtree-from, --node-type, etc.) will be ignored.\n"
              "Replacement uses the regex '\\bOLD_NAME\\b' to match whole words.")
    )
    
    # --- Parameters for Bio.Phylo mode (when --replace is NOT used) ---
    biophylo_group = parser.add_argument_group(
        'Bio.Phylo Mode Options (when --replace is NOT used)'
    )
    biophylo_group.add_argument(
        "--format", 
        default="newick", 
        help="Input tree format (e.g., newick, nexus, phyloxml).\nDefault: newick"
    )
    biophylo_group.add_argument(
        "--subtree-from",
        metavar="NODE_NAME",
        dest="subtree_from_node_name",
        help="Extract name list only from the subtree rooted at this internal node NODE_NAME."
    )
    biophylo_group.add_argument(
        "--node-type", 
        choices=["all", "internal", "leaf"], 
        default="all", 
        dest="node_type",
        help="Type of nodes to extract (Default: all)."
    )
    biophylo_group.add_argument(
        "--include-unnamed", 
        dest="unnamed_prefix", 
        nargs="?", # Optional argument for the flag itself
        const="UnnamedNode_", # Value if flag is present without an argument
        default=None, # Value if flag is not present
        help="Include unnamed nodes and assign them a prefix.\nE.g.: --include-unnamed MyPrefix_\nIf only the flag is provided, uses default prefix 'UnnamedNode_'."
    )
    biophylo_group.add_argument(
        "--min-depth", type=int, metavar="INT", help="Minimum depth of nodes to extract (subtree root is 0)."
    )
    biophylo_group.add_argument(
        "--max-depth", type=int, metavar="INT", help="Maximum depth of nodes to extract."
    )
    biophylo_group.add_argument(
        "--min-branch", type=float, dest="min_branch_length", metavar="FLOAT", help="Minimum branch length."
    )
    biophylo_group.add_argument(
        "--max-branch", type=float, dest="max_branch_length", metavar="FLOAT", help="Maximum branch length."
    )
    biophylo_group.add_argument(
        "--filter-name-regex", dest="name_regex", metavar="PATTERN", help="Regular expression to filter the node name list."
    )
    biophylo_group.add_argument(
        "--sort", choices=["asc", "desc", "none"], default="none", help="Sort the output node name list (Default: none)."
    )
    biophylo_group.add_argument(
        "--separator", default="\\n", help="Separator used when outputting the name list (Default: '\\n')."
    )

    # --- General Output Parameters ---
    parser.add_argument(
        "-o", "--output", 
        dest="output_filepath", 
        metavar="FILEPATH",
        help="File path to save the results (modified text or name list).\nIf not provided, prints to standard output."
    )
    parser.add_argument(
        "-q", "--quiet", 
        action="store_true", 
        help="Quiet mode, suppresses all informational messages to stderr."
    )
    args = parser.parse_args()

    # --- Logging functions ---
    def log_info(message):
        if not args.quiet: print(f"INFO: {message}", file=sys.stderr)
    def log_error(message):
        # Error messages are always printed to stderr
        print(f"ERROR: {message}", file=sys.stderr)

    # --- 1. Read input content (as string first for both modes) ---
    input_content_as_string = ""
    log_info(f"Attempting to read input content from '{args.tree_input}'...")
    if os.path.exists(args.tree_input):
        try:
            with open(args.tree_input, "r", encoding=args.encoding) as f:
                input_content_as_string = f.read() # Read entire file
            log_info(f"Successfully read content from file '{args.tree_input}' (encoding: {args.encoding}).")
        except Exception as e:
            log_error(f"Error reading file '{args.tree_input}': {e}"); sys.exit(1)
    else:
        log_info(f"Path '{args.tree_input}' is not an existing file, treating it as a direct input string.")
        input_content_as_string = args.tree_input
    
    if not input_content_as_string.strip(): # Check if string is empty or only whitespace
        log_error("Input content is empty or could not be read successfully."); sys.exit(1)

    # --- 2. Decide mode based on --replace ---
    if args.replacements:
        # --- MODE: Pure text replacement and output ---
        log_info("Detected --replace parameter. Entering plain text replacement mode.")
        log_info("Will directly modify the text of the input content and output it. Bio.Phylo related parameters will be ignored.")
        if args.subtree_from_node_name and not args.quiet : # Warn if --subtree-from is used in this mode
            log_info("WARNING: --subtree-from parameter is ineffective for output content in plain text replacement mode. The entire modified text will be output.")
        # Warn if other Bio.Phylo specific args are used
        if any([args.node_type != "all", args.unnamed_prefix, args.min_depth, args.max_depth, 
                args.min_branch_length, args.max_branch_length, args.name_regex, 
                args.sort != "none", args.separator != "\\n", args.format != "newick"]) and not args.quiet:
            log_info("WARNING: One or more Bio.Phylo mode parameters were specified but will be ignored in plain text replacement mode.")


        current_text_content = input_content_as_string
        replacement_happened = False # Flag to track if any replacement actually occurred
        for old_name, new_name in args.replacements:
            try:
                # Using \b for word boundaries. This helps to match 'NODE1' but not 'NODE10'.
                # Node names with special regex characters should be handled by re.escape.
                pattern = r'\b' + re.escape(old_name) + r'\b'
                modified_content, num_subs = re.subn(pattern, new_name, current_text_content)
                
                if num_subs > 0:
                    current_text_content = modified_content
                    replacement_happened = True
                    log_info(f"  Text replacement: '{old_name}' -> '{new_name}' ({num_subs} occurrences, using regex '\\b{re.escape(old_name)}\\b').")
                elif not args.quiet: # Log if a rule didn't match anything
                    log_info(f"  Text rule '{old_name}' -> '{new_name}': No matches found in text (using regex '\\b{re.escape(old_name)}\\b').")
            except re.error as e:
                log_error(f"Error during regex replacement (rule: '{old_name}' -> '{new_name}'): {e}"); sys.exit(1)
        
        if replacement_happened:
            log_info("Plain text replacement completed.")
        else:
            log_info("No text replacements were performed (no rules matched the content).")


        # Output the (potentially) modified text content
        if args.output_filepath:
            try:
                with open(args.output_filepath, "w", encoding=args.encoding) as outfile:
                    outfile.write(current_text_content)
                log_info(f"Modified text content saved to: '{args.output_filepath}'")
            except Exception as e:
                log_error(f"Error writing modified text to file '{args.output_filepath}': {e}"); sys.exit(1)
        else:
            sys.stdout.write(current_text_content) # Output the modified string to stdout
        
        sys.exit(0) # Exit after text replacement mode

    else:
        # --- MODE: Bio.Phylo based operations (name list extraction) ---
        log_info("No --replace parameter detected. Entering Bio.Phylo mode to extract node name list.")
        
        full_tree_obj = None
        try:
            tree_file_like = io.StringIO(input_content_as_string.strip()) # Use stripped content for Phylo.read
            full_tree_obj = Phylo.read(tree_file_like, args.format)
            log_info(f"(Bio.Phylo mode) Full tree (format: {args.format}) parsed successfully.")
        except Exception as e:
            log_error(f"(Bio.Phylo mode) Error parsing tree (format: {args.format}): {e}"); sys.exit(1)

        # Determine the starting clade for name extraction
        clade_for_operation = full_tree_obj.root 
        if args.subtree_from_node_name:
            log_info(f"(Bio.Phylo mode) Attempting to find specified subtree starting node: '{args.subtree_from_node_name}'...")
            found_target_node = None
            for clade_search in full_tree_obj.find_clades(): # Search in the full tree
                if clade_search.name == args.subtree_from_node_name:
                    found_target_node = clade_search; break
            if found_target_node:
                if not found_target_node.is_terminal(): # Must be an internal node
                    clade_for_operation = found_target_node
                    log_info(f"(Bio.Phylo mode) Found internal node '{args.subtree_from_node_name}'. Operations will start from this node.")
                else:
                    log_error(f"(Bio.Phylo mode) Node '{args.subtree_from_node_name}' was found, but it is a leaf node. --subtree-from requires an internal node name."); sys.exit(1)
            else:
                log_error(f"(Bio.Phylo mode) Could not find the specified --subtree-from node '{args.subtree_from_node_name}' in the tree."); sys.exit(1)
        
        # Extract names using Bio.Phylo logic
        names = get_names_from_phylotree( # Use the renamed function
            clade_for_operation, # Start extraction from this clade
            node_type=args.node_type.lower(),
            include_unnamed_prefix=args.unnamed_prefix,
            min_depth=args.min_depth,
            max_depth=args.max_depth,
            min_branch_length=args.min_branch_length,
            max_branch_length=args.max_branch_length,
            quiet=args.quiet 
        )
        
        if names is None: # Should not happen with current get_names_from_phylotree logic
            log_error("(Bio.Phylo mode) Internal error occurred during node name extraction."); sys.exit(1)
        if not names and not args.quiet: # Log if no names were extracted after filtering
            log_info(f"(Bio.Phylo mode) No node names were extracted from the target tree/subtree based on the filter criteria.")

        # Apply regex filtering to the list of names
        if args.name_regex:
            log_info(f"(Bio.Phylo mode) Applying regex filter to name list: '{args.name_regex}'")
            try:
                regex = re.compile(args.name_regex)
                original_count = len(names)
                names = [name for name in names if regex.search(name)]
                log_info(f"(Bio.Phylo mode) Regex filtering changed node count from {original_count} to {len(names)}.")
            except re.error as e:
                log_error(f"(Bio.Phylo mode) Invalid regex pattern '{args.name_regex}': {e}"); sys.exit(1)
                
        # Sort the list of names
        if args.sort == "asc": names.sort(); log_info("(Bio.Phylo mode) Node name list sorted in ascending order.")
        elif args.sort == "desc": names.sort(reverse=True); log_info("(Bio.Phylo mode) Node name list sorted in descending order.")
        elif args.sort == "none": log_info("(Bio.Phylo mode) Node name list not sorted (maintaining tree traversal order).")

        # Prepare output content (list of names)
        try:
            # Handle escaped characters like \n, \t in separator string
            separator = args.separator.encode('latin-1', 'backslashreplace').decode('unicode-escape')
        except Exception:
            log_info(f"(Bio.Phylo mode) Could not parse escape sequences in separator '{args.separator}', using it literally."); separator = args.separator # Fallback
        output_content = separator.join(names)

        # Output the list of names
        if args.output_filepath:
            log_info(f"(Bio.Phylo mode) Saving name list to file: '{args.output_filepath}' (encoding: {args.encoding})")
            try:
                with open(args.output_filepath, "w", encoding=args.encoding) as outfile:
                    outfile.write(output_content)
                    if separator == "\n" and output_content: outfile.write("\n") # Ensure newline at EOF if using newlines
                    elif not output_content and not args.quiet: # Log if output is empty
                         log_info(f"Final output content is empty, file '{args.output_filepath}' might be empty.")
                log_info("(Bio.Phylo mode) Name list saved successfully.")
            except Exception as e:
                log_error(f"(Bio.Phylo mode) Error writing name list to file '{args.output_filepath}': {e}"); sys.exit(1)
        else:
            # Print to stdout
            if names or output_content: # Print even if names is empty but separator might make output_content non-empty
                sys.stdout.write(output_content)
                if separator == "\n" and output_content : sys.stdout.write("\n") # Ensure newline at EOF
                sys.stdout.flush()
            elif not args.quiet : # Log if output is empty and not in quiet mode
                 log_info("(Bio.Phylo mode) Final extracted node name list is empty, nothing to print to standard output.")

if __name__ == "__main__":
    main()
