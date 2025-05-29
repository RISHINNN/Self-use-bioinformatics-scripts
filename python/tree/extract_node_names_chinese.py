
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
    """
    if not quiet:
        print(f"信息: (Bio.Phylo模式) 开始从节点 '{start_clade.name or 'UnnamedSubtreeRoot'}' 提取子树节点列表...", file=sys.stderr)

    node_names_extracted = []
    unnamed_counter = 0
    
    clade_depths = None
    if min_depth is not None or max_depth is not None:
        try:
            clade_depths = start_clade.depths(unit_branch_lengths=False)
            if not quiet:
                start_clade_relative_depth = clade_depths.get(start_clade, 'N/A')
                print(f"信息: (Bio.Phylo模式) 已计算相对于子树根的节点深度。子树根 '{start_clade.name or 'UnnamedSubtreeRoot'}' 的相对深度为 {start_clade_relative_depth}", file=sys.stderr)
        except Exception as e:
            if not quiet:
                print(f"警告: (Bio.Phylo模式) 计算节点深度时出错: {e}。将忽略深度筛选。", file=sys.stderr)
            clade_depths = None

    for clade in start_clade.find_clades(order='preorder'):
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
        
    if not quiet:
        print(f"信息: (Bio.Phylo模式) 从子树中初步提取到 {len(node_names_extracted)} 个符合基本条件的节点名称。", file=sys.stderr)
    return node_names_extracted

def main():
    parser = argparse.ArgumentParser(
        description="从系统发育树文件提取节点名称，或进行文本替换并输出修改后的文件内容。",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument(
        "tree_input", 
        help="必需：包含树信息的文件路径，或直接的树字符串 (当不作为文件路径找到时)。"
    )
    parser.add_argument(
        "--encoding", 
        default="utf-8", 
        help="用于读写文件的编码。\n默认: utf-8"
    )
    parser.add_argument(
        "--replace",
        nargs=2,
        action='append',
        metavar=('OLD_NAME', 'NEW_NAME'),
        dest='replacements',
        default=[],
        help=("进行纯文本替换：将文件中的 OLD_NAME 替换为 NEW_NAME。\n"
              "可以多次使用此选项进行多次替换。\n"
              "如果使用此选项，脚本将输出修改后的完整文件内容，\n"
              "并忽略下方大多数 Bio.Phylo 相关参数 (如 --format, --subtree-from, --node-type 等)。\n"
              "替换使用正则表达式 '\\bOLD_NAME\\b' 以匹配完整单词。")
    )
    
    # --- Parameters for Bio.Phylo mode (when --replace is NOT used) ---
    biophylo_group = parser.add_argument_group(
        'Bio.Phylo 模式选项 (当不使用 --replace 时)'
    )
    biophylo_group.add_argument(
        "--format", 
        default="newick", 
        help="输入树的格式 (例如: newick, nexus, phyloxml)。\n默认: newick"
    )
    biophylo_group.add_argument(
        "--subtree-from",
        metavar="NODE_NAME",
        dest="subtree_from_node_name",
        help="仅从此内部节点 NODE_NAME 为根的子树中提取名称列表。"
    )
    biophylo_group.add_argument(
        "--node-type", 
        choices=["all", "internal", "leaf"], 
        default="all", 
        dest="node_type",
        help="要提取的节点类型 (默认: all)。"
    )
    biophylo_group.add_argument(
        "--include-unnamed", 
        dest="unnamed_prefix", 
        nargs="?", 
        const="UnnamedNode_", 
        default=None,
        help="包含未命名节点，并为其指定前缀 (默认: UnnamedNode_ 无参数时)。"
    )
    biophylo_group.add_argument(
        "--min-depth", type=int, metavar="INT", help="提取节点的最小深度 (子树根为0)。"
    )
    biophylo_group.add_argument(
        "--max-depth", type=int, metavar="INT", help="提取节点的最大深度。"
    )
    biophylo_group.add_argument(
        "--min-branch", type=float, dest="min_branch_length", metavar="FLOAT", help="最小枝长。"
    )
    biophylo_group.add_argument(
        "--max-branch", type=float, dest="max_branch_length", metavar="FLOAT", help="最大枝长。"
    )
    biophylo_group.add_argument(
        "--filter-name-regex", dest="name_regex", metavar="PATTERN", help="用于筛选节点名称列表的正则表达式。"
    )
    biophylo_group.add_argument(
        "--sort", choices=["asc", "desc", "none"], default="none", help="对输出的节点名称列表进行排序 (默认: none)。"
    )
    biophylo_group.add_argument(
        "--separator", default="\\n", help="输出名称列表时使用的分隔符 (默认: '\\n')."
    )

    # --- General Output Parameters ---
    parser.add_argument(
        "-o", "--output", 
        dest="output_filepath", 
        metavar="FILEPATH",
        help="将结果 (修改后的文本 或 名称列表) 保存到的文件路径。\n如果未提供，则打印到标准输出。"
    )
    parser.add_argument(
        "-q", "--quiet", 
        action="store_true", 
        help="安静模式，禁止所有信息性消息输出到 stderr。"
    )
    args = parser.parse_args()

    def log_info(message):
        if not args.quiet: print(f"信息: {message}", file=sys.stderr)
    def log_error(message):
        print(f"错误: {message}", file=sys.stderr)

    # --- 1. Read input content (as string first for both modes) ---
    input_content_as_string = ""
    log_info(f"尝试从 '{args.tree_input}' 读取输入内容...")
    if os.path.exists(args.tree_input):
        try:
            with open(args.tree_input, "r", encoding=args.encoding) as f:
                input_content_as_string = f.read() # Read entire file
            log_info(f"成功从文件 '{args.tree_input}' 读取内容 (编码: {args.encoding})。")
        except Exception as e:
            log_error(f"读取文件 '{args.tree_input}' 时出错: {e}"); sys.exit(1)
    else:
        log_info(f"路径 '{args.tree_input}' 不是一个存在的文件，将其作为直接的输入字符串处理。")
        input_content_as_string = args.tree_input
    
    if not input_content_as_string.strip(): # Check if string is empty or only whitespace
        log_error("输入内容为空或未能成功读取。"); sys.exit(1)

    # --- 2. Decide mode based on --replace ---
    if args.replacements:
        # --- MODE: Pure text replacement and output ---
        log_info("检测到 --replace 参数。进入纯文本替换模式。")
        log_info("将直接修改输入内容的文本并输出。Bio.Phylo 相关参数将被忽略。")
        if args.subtree_from_node_name and not args.quiet :
            log_info("警告: --subtree-from 参数在纯文本替换模式下对输出内容无效。将输出整个修改后的文本。")
        if any([args.node_type != "all", args.unnamed_prefix, args.min_depth, args.max_depth, 
                args.min_branch_length, args.max_branch_length, args.name_regex, 
                args.sort != "none", args.separator != "\\n", args.format != "newick"]) and not args.quiet:
            log_info("警告: 一个或多个 Bio.Phylo 模式的参数被指定，但在纯文本替换模式下它们将被忽略。")


        current_text_content = input_content_as_string
        replacement_happened = False
        for old_name, new_name in args.replacements:
            try:
                # Using \b for word boundaries. Assumes node names are "word-like".
                # For names with non-alphanumeric characters (e.g., hyphens), this might need adjustment
                # or users should use more specific regex if this is too broad/narrow.
                pattern = r'\b' + re.escape(old_name) + r'\b'
                modified_content, num_subs = re.subn(pattern, new_name, current_text_content)
                
                if num_subs > 0:
                    current_text_content = modified_content
                    replacement_happened = True
                    log_info(f"  文本替换: '{old_name}' -> '{new_name}' ({num_subs} 次，使用正则 '\\b{re.escape(old_name)}\\b')。")
                elif not args.quiet:
                    log_info(f"  文本规则 '{old_name}' -> '{new_name}': 未在文本中找到匹配项 (使用正则 '\\b{re.escape(old_name)}\\b')。")
            except re.error as e:
                log_error(f"正则表达式替换时出错 (规则: '{old_name}' -> '{new_name}'): {e}"); sys.exit(1)
        
        if replacement_happened:
            log_info("纯文本替换完成。")
        else:
            log_info("未执行任何文本替换（没有规则匹配到内容）。")


        if args.output_filepath:
            try:
                with open(args.output_filepath, "w", encoding=args.encoding) as outfile:
                    outfile.write(current_text_content)
                log_info(f"修改后的文本内容已保存到: '{args.output_filepath}'")
            except Exception as e:
                log_error(f"写入修改后的文本到文件 '{args.output_filepath}' 时出错: {e}"); sys.exit(1)
        else:
            sys.stdout.write(current_text_content) # Output the modified string
        
        sys.exit(0)

    else:
        # --- MODE: Bio.Phylo based operations (name list extraction) ---
        log_info("未检测到 --replace 参数。进入 Bio.Phylo 模式以提取节点名称列表。")
        
        full_tree_obj = None
        try:
            tree_file_like = io.StringIO(input_content_as_string.strip()) # Use stripped content
            full_tree_obj = Phylo.read(tree_file_like, args.format)
            log_info(f"(Bio.Phylo模式) 完整树 (格式: {args.format}) 解析成功。")
        except Exception as e:
            log_error(f"(Bio.Phylo模式) 解析树 (格式: {args.format}) 时出错: {e}"); sys.exit(1)

        clade_for_operation = full_tree_obj.root 
        if args.subtree_from_node_name:
            log_info(f"(Bio.Phylo模式) 尝试查找指定的子树起始节点: '{args.subtree_from_node_name}'...")
            found_target_node = None
            for clade_search in full_tree_obj.find_clades():
                if clade_search.name == args.subtree_from_node_name:
                    found_target_node = clade_search; break
            if found_target_node:
                if not found_target_node.is_terminal():
                    clade_for_operation = found_target_node
                    log_info(f"(Bio.Phylo模式) 已找到内部节点 '{args.subtree_from_node_name}'。操作将从此节点开始。")
                else:
                    log_error(f"(Bio.Phylo模式) 节点 '{args.subtree_from_node_name}' 被找到，但它是一个叶节点。--subtree-from 需要一个内部节点名称。"); sys.exit(1)
            else:
                log_error(f"(Bio.Phylo模式) 未能在树中找到指定的 --subtree-from 节点 '{args.subtree_from_node_name}'。"); sys.exit(1)
        
        names = get_names_from_phylotree( # Use the renamed function
            clade_for_operation,
            node_type=args.node_type.lower(),
            include_unnamed_prefix=args.unnamed_prefix,
            min_depth=args.min_depth,
            max_depth=args.max_depth,
            min_branch_length=args.min_branch_length,
            max_branch_length=args.max_branch_length,
            quiet=args.quiet 
        )
        
        if names is None: log_error("(Bio.Phylo模式) 提取节点名称时发生内部错误。"); sys.exit(1)
        if not names and not args.quiet:
            log_info(f"(Bio.Phylo模式) 根据筛选条件，未能从目标树/子树中提取到任何节点名称。")

        if args.name_regex:
            log_info(f"(Bio.Phylo模式) 应用正则表达式筛选名称列表: '{args.name_regex}'")
            try:
                regex = re.compile(args.name_regex)
                original_count = len(names)
                names = [name for name in names if regex.search(name)]
                log_info(f"(Bio.Phylo模式) 正则筛选后，节点数量从 {original_count} 变为 {len(names)}。")
            except re.error as e:
                log_error(f"(Bio.Phylo模式) 无效的正则表达式模式 '{args.name_regex}': {e}"); sys.exit(1)
                
        if args.sort == "asc": names.sort(); log_info("(Bio.Phylo模式) 已按升序对节点名称列表进行排序。")
        elif args.sort == "desc": names.sort(reverse=True); log_info("(Bio.Phylo模式) 已按降序对节点名称列表进行排序。")
        elif args.sort == "none": log_info("(Bio.Phylo模式) 节点名称列表未排序。")

        try:
            separator = args.separator.encode('latin-1', 'backslashreplace').decode('unicode-escape')
        except Exception:
            log_info(f"(Bio.Phylo模式) 无法解析分隔符 '{args.separator}'，将按原样使用。"); separator = args.separator
        output_content = separator.join(names)

        if args.output_filepath:
            log_info(f"(Bio.Phylo模式) 将名称列表保存到文件: '{args.output_filepath}' (编码: {args.encoding})")
            try:
                with open(args.output_filepath, "w", encoding=args.encoding) as outfile:
                    outfile.write(output_content)
                    if separator == "\n" and output_content: outfile.write("\n")
                    elif not output_content and not args.quiet:
                         log_info(f"最终输出内容为空，文件 '{args.output_filepath}' 可能为空。")
                log_info("(Bio.Phylo模式) 名称列表已成功保存。")
            except Exception as e:
                log_error(f"(Bio.Phylo模式) 写入名称列表到文件 '{args.output_filepath}' 时出错: {e}"); sys.exit(1)
        else:
            if names or output_content: 
                sys.stdout.write(output_content)
                if separator == "\n" and output_content : sys.stdout.write("\n")
                sys.stdout.flush()
            elif not args.quiet :
                 log_info("(Bio.Phylo模式) 最终提取到的节点名称列表为空，无内容打印到标准输出。")

if __name__ == "__main__":
    main()
