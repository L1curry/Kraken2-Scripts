#!/usr/bin/env python3
"""
在一个目录下的所有 FASTA 文件中检查指定 ID 存在于哪个文件，并生成可扩展的输出结果。

fasta_id_mapper.py:
    扫描目录，读取 FASTA 文件中的序列 ID，将查询列表中的 ID 映射到出现它们的文件。

Usage:
    python fasta_id_mapper.py \
        --input-dir <INPUT_DIR> \
        --query-file <QUERY_FILE> \
        [--segment] [--recursive] [--output <OUT_FILE>] [--format {text,json}]

Options:
    -i, --input-dir      输入目录，包含要扫描的 FASTA 文件
    -q, --query-file     包含查询序列 ID (每行一个) 的文本文件
    -s, --segment        启用分段匹配：将 FASTA ID 按非字母数字分隔，匹配独立部分或前缀
    -r, --recursive      递归扫描子目录
    -o, --output         输出映射结果到指定文件 (默认 stdout)
    -f, --format         输出格式：text (默认) 或 json
    -h, --help           显示帮助并退出
"""

import sys
import argparse
import re
import json
from pathlib import Path
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(description="Map sequence IDs in FASTA files to query IDs.")
    parser.add_argument('-i', '--input-dir', required=True, type=Path,
                        help='输入目录，包含 FASTA 文件')
    parser.add_argument('-q', '--query-file', required=True, type=Path,
                        help='查询序列 ID 列表文件')
    parser.add_argument('-s', '--segment', action='store_true',
                        help='使用分段匹配模式')
    parser.add_argument('-r', '--recursive', action='store_true',
                        help='递归扫描子目录')
    parser.add_argument('-o', '--output', type=Path,
                        help='输出文件路径 (默认写到 stdout)')
    parser.add_argument('-f', '--format', choices=('text', 'json'), default='text',
                        help='输出格式: text or json')
    return parser.parse_args()


def load_query_ids(path):
    try:
        return [line.strip() for line in path.read_text(encoding='utf-8').splitlines() if line.strip()]
    except Exception as e:
        sys.exit(f"Error loading query file '{path}': {e}")


def find_fasta_files(directory, recursive=False):
    pattern = '**/*' if recursive else '*'
    return [p for p in directory.glob(pattern) if p.is_file()]


def is_fasta_file(path, max_lines=10):
    try:
        with path.open('r', encoding='utf-8', errors='ignore') as f:
            for _ in range(max_lines):
                if f.readline().startswith('>'):
                    return True
    except Exception:
        return False
    return False


def extract_ids(path):
    ids = set()
    with path.open('r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('>'):
                header = line[1:].strip()
                ids.add(header.split()[0])
    return ids


def build_pattern(qid):
    # 分段或前缀匹配：^qid($|非字母数字)
    return re.compile(rf'^{re.escape(qid)}($|[^0-9A-Za-z])')


def map_ids_to_files(query_ids, fasta_ids_map, segment=False):
    result = defaultdict(list)
    patterns = {qid: build_pattern(qid) for qid in query_ids} if segment else {}
    for qid in query_ids:
        for fname, ids in fasta_ids_map.items():
            for sid in ids:
                if sid == qid or (segment and patterns[qid].match(sid)):
                    result[qid].append(str(fname))
                    break
    return result


def summarize(mapping):
    total = len(mapping)
    found = sum(bool(files) for files in mapping.values())
    missing = [qid for qid, files in mapping.items() if not files]
    multi = {qid: files for qid, files in mapping.items() if len(files) > 1}
    return {
        'total': total,
        'found': found,
        'missing_count': len(missing),
        'missing_ids': missing,
        'multi_files': multi
    }


def output_text(mapping, summary, out_stream):
    out_stream.write("Sequence ID mapping results:\n")
    for qid, files in mapping.items():
        line = f"{qid}: {', '.join(files) if files else 'NOT FOUND'}"
        out_stream.write(line + "\n")
    out_stream.write("\nSummary:\n")
    out_stream.write(f"Total queries: {summary['total']}\n")
    out_stream.write(f"Found: {summary['found']}\n")
    out_stream.write(f"Missing: {summary['missing_count']}\n")
    if summary['missing_count']:
        out_stream.write(f"Missing IDs: {', '.join(summary['missing_ids'])}\n")
    out_stream.write(f"IDs in multiple files: {len(summary['multi_files'])}\n")
    for qid, files in summary['multi_files'].items():
        out_stream.write(f"  {qid}: {', '.join(files)}\n")


def output_json(mapping, summary, out_stream):
    json.dump({'mapping': mapping, 'summary': summary}, out_stream, ensure_ascii=False, indent=2)


def main():
    args = parse_args()
    queries = load_query_ids(args.query_file)

    fasta_files = find_fasta_files(args.input_dir, args.recursive)
    fasta_ids_map = {}
    for path in fasta_files:
        if is_fasta_file(path):
            fasta_ids_map[path] = extract_ids(path)

    mapping = map_ids_to_files(queries, fasta_ids_map, segment=args.segment)
    summary = summarize(mapping)

    # 输出
    out_stream = args.output.open('w', encoding='utf-8') if args.output else sys.stdout
    if args.format == 'json':
        output_json(mapping, summary, out_stream)
    else:
        output_text(mapping, summary, out_stream)
    if args.output:
        out_stream.close()

if __name__ == '__main__':
    main()
