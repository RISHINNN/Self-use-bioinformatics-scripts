#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# 默认参数
my @fields = ();
my $delimiter = ';';
my $input_file = '';
my $output_file = '';
my $help = 0;

# 获取命令行参数
GetOptions(
    'field|f=s@'     => \@fields,      # 可以指定多个字段，如 -f product= -f gene=
    'delimiter|d=s'  => \$delimiter,   # 终止分隔符，默认是 ;
    'input|i=s'      => \$input_file,  # 输入文件
    'output|o=s'     => \$output_file, # 输出文件（可选）
    'help|h'         => \$help,        # 帮助信息
) or die "Error in command line arguments\n";

# 显示帮助信息
if ($help || scalar(@fields) == 0) {
    print_help();
    exit;
}

# 确定输入源：管道���文件或STDIN
my $fh_in;
if ($input_file) {
    # 从指定文件读取
    open($fh_in, '<', $input_file) or die "Cannot open input file '$input_file': $!\n";
} elsif (! -t STDIN) {
    # 从管道读取（STDIN不是终端）
    $fh_in = *STDIN;
} else {
    die "No input provided.  Use -i <file> or pipe input\n";
}

# 打开输出文件（如果指定）
my $fh_out;
if ($output_file) {
    open($fh_out, '>', $output_file) or die "Cannot open output file '$output_file': $!\n";
} else {
    $fh_out = *STDOUT;
}

# 处理每一行
while (my $line = <$fh_in>) {
    chomp $line;
    
    # 存储提取的结果
    my @results = ();
    
    # 对每个指定的字段进行提取
    foreach my $field (@fields) {
        my $value = extract_value($line, $field, $delimiter);
        push @results, $value if defined $value;
    }
    
    # 输出结果（如果有提取到内容）
    if (@results) {
        print $fh_out join("\t", @results) . "\n";
    }
}

close($fh_in) if $input_file;
close($fh_out) if $output_file;

# 提取函数
sub extract_value {
    my ($text, $field_name, $end_delimiter) = @_;
    
    # 转义特殊字符
    my $escaped_field = quotemeta($field_name);
    my $escaped_delim = quotemeta($end_delimiter);
    
    # 匹配模式：字段名后面的内容，直到遇到分隔符或行尾
    if ($text =~ /${escaped_field}([^${escaped_delim}]+)/) {
        my $value = $1;
        # 去除前后空格
        $value =~ s/^\s+|\s+$//g;
        return $value;
    }
    
    return undef;
}

# 帮助信息
sub print_help {
    print <<'HELP';
用法: perl zhua.pl -f <field1> [-f <field2> ...  ] [options]
      或通过管道:  cat file | perl zhua.pl -f <field1>

必需参数:
    -f, --field <string>     要提取的字段标识符（可以指定多次）
                             例如: -f product= -f gene= -f protein_id=

可选参数:
    -i, --input <file>       输入文件路径（不指定则从管道或STDIN读取）
    -d, --delimiter <char>   字段终止符（默认: ;）
    -o, --output <file>      输出文件路径（默认: 标准输出）
    -h, --help               显示此帮助信息

示例:  
    # 从文件读取
    perl zhua.pl -i input.txt -f product=
    
    # 从管道读取
    cat input.txt | perl zhua.pl -f product=
    
    # 从 grep 结果读取
    grep "CDS" input.txt | perl zhua.pl -f product= -f gene=
    
    # 提取多个字段并保存
    cat input.txt | perl zhua.pl -f product= -f gene= -o output.txt
    
    # 链式管道操作
    cat input.txt | grep "isoform" | perl zhua.pl -f product= | sort | uniq

HELP
}
