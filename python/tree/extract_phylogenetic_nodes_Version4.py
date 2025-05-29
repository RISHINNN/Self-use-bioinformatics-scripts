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
    log_func=None # Pass a logging function for verbose output
):
    """
    (Used when --replace is NOT active)
    From a Bio.Phylo.Clade object, extracts node names based on criteria.
    Depth and traversal are relative to start_clade.
    """
    if log_func:
        log_func(f"(Bio.Phylo mode) Starting to extract node list from subtree rooted at '{start_clade.name or 'UnnamedSubtreeRoot'}'...")

    node_names_extracted = []
    unnamed_counter = 0
    
    clade_depths = None
    if min_depth is not None or max_depth is not None:
        try:
            clade_depths = start_clade.depths(unit_branch_lengths=False) # Depths relative to start_clade
            if log_func:
                start_clade_relative_depth = clade_depths.get(start_clade, 'N/A')
                log_func(f"(Bio.Phylo mode) Calculated node depths relative to the subtree root. Subtree root '{start_clade.name or 'UnnamedSubtreeRoot'}' relative depth is {start_clade_relative_depth}")
        except Exception as e:
            # Errors from here should probably go to stderr directly or a dedicated error logger
            print(f"WARNING: (Bio.Phylo mode) Error calculating node depths: {e}. Depth filtering will be ignored.", file=sys.stderr)
            clade_depths = None

    for clade in start_clade.find_clades(order='preorder'): # Iterate from start_clade
        current_name = clade.name
        
        if current_name is None and include_unnamed_prefix:
            current_name = f"{include_unnamed_prefix}{unnamed_counter}"
            unnamed_counter += 1
        
        if current_name is None: 
            continue

        is_correct_type = False
        if node_type == "all": is_correct_type = True
        elif node_type == "internal":
            if not clade.is_terminal(): is_correct_type = True
        elif node_type == "leaf":
            if clade.is_terminal(): is_correct_type = True
        
        if not is_correct_type: continue

        if clade_depths:
            depth = clade_depths.get(clade) 
            if depth is not None: 
                if min_depth is not None and depth < min_depth: continue
                if max_depth is not None and depth > max_depth: continue

        branch_len = clade.branch_length
        if min_branch_length is not None or max_branch_length is not None:
            if branch_len is None: continue 
            if min_branch_length is not None and branch_len < min_branch_length: continue
            if max_branch_length is not None and branch_len > max_branch_length: continue
        
        node_names_extracted.append(current_name)
        
    if log_func:
        log_func(f"(Bio.Phylo mode) Initially extracted {len(node_names_extracted)} node names from the subtree based on basic criteria.")
    return node_names_extracted

def main():
    parser = argparse.ArgumentParser(
        description="Extracts node names from a phylogenetic tree file or performs text replacement and outputs the modified file content.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument(
        "tree_input", 
        help="Required: File path containing tree information, or a direct tree string (if not found as a file path)."
    )
    parser.add_argument(
        "--encoding", 
        default="utf-8", 
        help="Encoding for reading and writing files.\nDefault: utf-8"
    )
    parser.add_argument(
        "--replace",
        nargs=2,
        action='append',
        metavar=('OLD_NAME', 'NEW_NAME'),
        dest='replacements',
        default=[],
        help=("Perform plain text replacement: replaces OLD_NAME with NEW_NAME in the file.\n"
              "This option can be used multiple times for several replacements.\n"
              "If this option is used, the script will output the modified full file content,\n"
              "and most Bio.Phylo related parameters below will be ignored.\n"
              "Replacement uses the regex '\\bOLD_NAME\\b' to match whole words.")
    )
    
    biophylo_group = parser.add_argument_group(
        'Bio.Phylo Mode Options (when --replace is NOT used)'
    )
    biophylo_group.add_argument(
        "--format", default="newick", help="Input tree format (e.g., newick, nexus, phyloxml).\nDefault: newick"
    )
    biophylo_group.add_argument(
        "--subtree-from", metavar="NODE_NAME", dest="subtree_from_node_name",
        help="Extract name list only from the subtree rooted at this internal node NODE_NAME."
    )
    biophylo_group.add_argument(
        "--node-type", choices=["all", "internal", "leaf"], default="all", dest="node_type",
        help="Type of nodes to extract (Default: all)."
    )
    biophylo_group.add_argument(
        "--include-unnamed", dest="unnamed_prefix", nargs="?", const="UnnamedNode_", default=None,
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
        help="Quiet mode, suppresses all informational and warning (from log_info) messages to stderr. Errors (from log_error) will still be shown."
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Verbose mode, shows detailed informational messages to stderr. Overridden by --quiet."
    )
    args = parser.parse_args()

    # --- Logging functions ---
    def log_info(message):
        if not args.quiet and args.verbose: # Only print if not quiet AND verbose
            print(f"INFO: {message}", file=sys.stderr)
    
    def log_error(message): # Errors are always printed
        print(f"ERROR: {message}", file=sys.stderr)

    # --- 1. Read input content (as string first for both modes) ---
    input_content_as_string = ""
    log_info(f"Attempting to read input content from '{args.tree_input}'...") # This log_info call will now respect verbose/quiet
    if os.path.exists(args.tree_input):
        try:
            with open(args.tree_input, "r", encoding=args.encoding) as f:
                input_content_as_string = f.read()
            log_info(f"Successfully read content from file '{args.tree_input}' (encoding: {args.encoding}).")
        except Exception as e:
            log_error(f"Error reading file '{args.tree_input}': {e}"); sys.exit(1)
    else:
        log_info(f"Path '{args.tree_input}' is not an existing file, treating it as a direct input string.")
        input_content_as_string = args.tree_input
    
    if not input_content_as_string.strip():
        log_error("Input content is empty or could not be read successfully."); sys.exit(1)

    # --- 2. Decide mode based on --replace ---
    if args.replacements:
        log_info("Detected --replace parameter. Entering plain text replacement mode.")
        log_info("Will directly modify the text of the input content and output it. Bio.Phylo related parameters will be ignored.")
        if args.subtree_from_node_name: # This is a specific warning if a likely incompatible flag is used
             if not args.quiet and args.verbose: # Only show this warning if verbose
                print("WARNING: --subtree-from parameter is ineffective for output content in plain text replacement mode. The entire modified text will be output.", file=sys.stderr)
        
        current_text_content = input_content_as_string
        replacement_happened = False
        for old_name, new_name in args.replacements:
            try:
                pattern = r'\b' + re.escape(old_name) + r'\b'
                modified_content, num_subs = re.subn(pattern, new_name, current_text_content)
                
                if num_subs > 0:
                    current_text_content = modified_content
                    replacement_happened = True
                    log_info(f"  Text replacement: '{old_name}' -> '{new_name}' ({num_subs} occurrences, using regex '\\b{re.escape(old_name)}\\b').")
                else: # Log only if verbose and no replacement
                    log_info(f"  Text rule '{old_name}' -> '{new_name}': No matches found in text (using regex '\\b{re.escape(old_name)}\\b').")
            except re.error as e:
                log_error(f"Error during regex replacement (rule: '{old_name}' -> '{new_name}'): {e}"); sys.exit(1)
        
        if replacement_happened:
            log_info("Plain text replacement completed.")
        else:
            log_info("No text replacements were performed (no rules matched the content).")

        if args.output_filepath:
            try:
                with open(args.output_filepath, "w", encoding=args.encoding) as outfile:
                    outfile.write(current_text_content)
                log_info(f"Modified text content saved to: '{args.output_filepath}'")
            except Exception as e:
                log_error(f"Error writing modified text to file '{args.output_filepath}': {e}"); sys.exit(1)
        else:
            sys.stdout.write(current_text_content)
        sys.exit(0)

    else:
        log_info("No --replace parameter detected. Entering Bio.Phylo mode to extract node name list.")
        
        full_tree_obj = None
        try:
            tree_file_like = io.StringIO(input_content_as_string.strip())
            full_tree_obj = Phylo.read(tree_file_like, args.format)
            log_info(f"(Bio.Phylo mode) Full tree (format: {args.format}) parsed successfully.")
        except Exception as e:
            log_error(f"(Bio.Phylo mode) Error parsing tree (format: {args.format}): {e}"); sys.exit(1)

        clade_for_operation = full_tree_obj.root 
        if args.subtree_from_node_name:
            log_info(f"(Bio.Phylo mode) Attempting to find specified subtree starting node: '{args.subtree_from_node_name}'...")
            found_target_node = None
            for clade_search in full_tree_obj.find_clades():
                if clade_search.name == args.subtree_from_node_name:
                    found_target_node = clade_search; break
            if found_target_node:
                if not found_target_node.is_terminal():
                    clade_for_operation = found_target_node
                    log_info(f"(Bio.Phylo mode) Found internal node '{args.subtree_from_node_name}'. Operations will start from this node.")
                else:
                    log_error(f"(Bio.Phylo mode) Node '{args.subtree_from_node_name}' was found, but it is a leaf node. --subtree-from requires an internal node name."); sys.exit(1)
            else:
                log_error(f"(Bio.Phylo mode) Could not find the specified --subtree-from node '{args.subtree_from_node_name}' in the tree."); sys.exit(1)
        
        names = get_names_from_phylotree(
            clade_for_operation,
            node_type=args.node_type.lower(),
            include_unnamed_prefix=args.unnamed_prefix,
            min_depth=args.min_depth,
            max_depth=args.max_depth,
            min_branch_length=args.min_branch_length,
            max_branch_length=args.max_branch_length,
            log_func=log_info # Pass the log_info function for internal verbose logging
        )
        
        if names is None: log_error("(Bio.Phylo mode) Internal error occurred during node name extraction."); sys.exit(1)
        log_info(f"(Bio.Phylo mode) Total names extracted after all filters: {len(names)}")


        if args.name_regex:
            log_info(f"(Bio.Phylo mode) Applying regex filter to name list: '{args.name_regex}'")
            try:
                regex = re.compile(args.name_regex)
                original_count = len(names)
                names = [name for name in names if regex.search(name)]
                log_info(f"(Bio.Phylo mode) Regex filtering changed node count from {original_count} to {len(names)}.")
            except re.error as e:
                log_error(f"(Bio.Phylo mode) Invalid regex pattern '{args.name_regex}': {e}"); sys.exit(1)
                
        if args.sort == "asc": names.sort(); log_info("(Bio.Phylo mode) Node name list sorted in ascending order.")
        elif args.sort == "desc": names.sort(reverse=True); log_info("(Bio.Phylo mode) Node name list sorted in descending order.")
        elif args.sort == "none": log_info("(Bio.Phylo mode) Node name list not sorted.")

        try:
            separator = args.separator.encode('latin-1', 'backslashreplace').decode('unicode-escape')
        except Exception:
            log_info(f"(Bio.Phylo mode) Could not parse escape sequences in separator '{args.separator}', using it literally."); separator = args.separator
        output_content = separator.join(names)

        if args.output_filepath:
            log_info(f"(Bio.Phylo mode) Saving name list to file: '{args.output_filepath}' (encoding: {args.encoding})")
            try:
                with open(args.output_filepath, "w", encoding=args.encoding) as outfile:
                    outfile.write(output_content)
                    if separator == "\n" and output_content: outfile.write("\n")
                log_info("(Bio.Phylo mode) Name list saved successfully.")
            except Exception as e:
                log_error(f"(Bio.Phylo mode) Error writing name list to file '{args.output_filepath}': {e}"); sys.exit(1)
        else:
            if names or output_content: 
                sys.stdout.write(output_content)
                if separator == "\n" and output_content : sys.stdout.write("\n")
                sys.stdout.flush()
            else: # Only log this if not quiet and no output was produced
                 log_info("(Bio.Phylo mode) Final extracted node name list is empty, nothing to print to standard output.")

if __name__ == "__main__":
    main()