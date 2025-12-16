#!/usr/bin/env python3
"""
脚本功能：根据映射文件，为指定的 FASTA 文件中所有序列的 header 行（以 `>` 开头）添加或替换后缀字符串。
用法示例：
    python AddSuffix_toFastaTitle.py -newtitle mapping.txt -p before
    python AddSuffix_toFastaTitle.py -newtitle mapping.txt -p _

参数说明：
  -newtitle    映射文件路径，制表符或空白分隔的两列：
               第一列为 FASTA 文件名，第二列为要添加的字符串
  -p, --position  指定添加位置：
               before  在 header（`>`）后插入字符串；
               after   在 header 行末尾插入字符串；
               其他   视作匹配字符或模式，将首次匹配到的部分替换为该字符串。
"""
import argparse
import os
import sys
from datetime import datetime

def parse_args():
    parser = argparse.ArgumentParser(
        description='根据映射文件，为FASTA header添加或替换后缀字符串',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '-newtitle', required=True, metavar='FILE',
        help='映射文件：每行两列，第一列FASTA文件名，第二列要添加的字符串'
    )
    parser.add_argument(
        '-p', '--position', required=True,
        help='添加位置参数：\n'
             '  before：在`>`后插入字符串；\n'
             '  after：在行尾插入字符串；\n'
             '  其他：将首次匹配到的该参数字符串替换为新字符串'
    )
    return parser.parse_args()

def read_mapping(path):
    """
    读取映射文件，返回 (filename, suffix) 列表。
    """
    mapping = []
    with open(path, 'r', encoding='utf-8') as f:
        for lineno, line in enumerate(f, start=1):
            parts = line.strip().split()
            if len(parts) != 2:
                print(f"[警告] 映射文件 '{path}' 第{lineno}行格式不符：{line.strip()}")
                continue
            mapping.append((parts[0], parts[1]))
    return mapping


def process_fasta(fasta_file, suffix, position):
    """
    根据 position 参数，修改 FASTA 文件 header：
      - before：在 ">" 后添加 suffix
      - after：在行尾添加 suffix
      - 其他：将 header 中首次匹配到 position 的部分替换为 suffix
    返回修改的序列数。
    """
    count = 0
    lines = []
    with open(fasta_file, 'r', encoding='utf-8') as fr:
        for line in fr:
            if line.startswith('>'):
                header = line[1:].rstrip('\n')
                if position == 'before':
                    new_header = suffix + header
                elif position == 'after':
                    new_header = header + suffix
                else:
                    if position in header:
                        new_header = header.replace(position, suffix, 1)
                    else:
                        # 未找到匹配，保留原 header 并直接加后缀
                        new_header = header + suffix
                lines.append('>' + new_header + '\n')
                count += 1
            else:
                lines.append(line)
    # 写回文件
    with open(fasta_file, 'w', encoding='utf-8') as fw:
        fw.writelines(lines)
    return count


def main():
    args = parse_args()

    # 检查映射文件
    if not os.path.isfile(args.newtitle):
        print(f"[错误] 映射文件不存在: {args.newtitle}")
        sys.exit(1)

    # 开始计时
    start = datetime.now()
    print(f"[信息] 脚本开始: {start.strftime('%Y-%m-%d %H:%M:%S')} \n")

    mapping = read_mapping(args.newtitle)
    if not mapping:
        print("[错误] 未读取到有效的映射关系。")
        sys.exit(1)

    for fasta_file, suffix in mapping:
        if not os.path.isfile(fasta_file):
            print(f"[警告] FASTA 文件不存在，跳过: {fasta_file}")
            continue
        print(f"[信息] 处理 '{fasta_file}'，添加字符串 '{suffix}' (position={args.position})...")
        modified = process_fasta(fasta_file, suffix, args.position)
        print(f"[信息] 共修改 {modified} 条序列 header。\n")

    # 结束计时
    end = datetime.now()
    print(f"[信息] 脚本完成: {end.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"[信息] 总耗时: {end - start}")

if __name__ == '__main__':
    main()
